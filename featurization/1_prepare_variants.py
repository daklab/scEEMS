#!/usr/bin/env python
"""
Read a table of variants, check them against the reference genome and pair each with the genes it is
scored for. The later steps read the two files written here.

Input: a tab-separated file with a header row and one variant per row (GRCh38, 1-based positions), either
    chrom  pos  ref  alt  [name]  [gene_id]
or a variant_id column (chrom:pos:ref:alt) instead of the first four. chrom may be written with or without
"chr"; only autosomes are scored (the models were trained on chromosomes 1-22). SNVs, insertions,
deletions and multi-base substitutions are accepted. Optional columns: name (e.g. an rsID, carried to the
outputs) and gene_id (Ensembl gene ID): a row with a gene_id pairs its variant with that gene; a variant
without one is paired with every gene whose TSS lies within --window bp of it.

Alleles are oriented to the genome: a variant whose REF does not match GRCh38 but whose ALT does has REF
and ALT swapped (and is reported); a variant matching neither is dropped. The sequence-model features
need this orientation: they compare the genome sequence with the sequence carrying ALT, so a variant
written with the genome allele as ALT would get no effect.

usage:   python 1_prepare_variants.py VARIANTS.tsv OUT_DIR --fasta GRCh38.fa [--window 1000000]
outputs: OUT_DIR/variants.tsv   variant_id, CHR, BP, REF, ALT, SNP (name, else variant_id), input_variant
         OUT_DIR/pairs.tsv      variant_id, gene_id, gene_name
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("variants", help="tab-separated table of variants (see above)")
ap.add_argument("out_dir")
ap.add_argument("--fasta", required=True, help="GRCh38 FASTA (featurization/README.md)")
ap.add_argument("--window", type=int, default=1_000_000, help="largest variant-TSS distance when pairing by window")
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

raw = pd.read_csv(a.variants, sep="\t", dtype=str, keep_default_na=False)
col = {c.strip().lower(): c for c in raw.columns}
if "variant_id" in col:
    parts = raw[col["variant_id"]].str.split(":", expand=True)
    if parts.shape[1] != 4:
        sys.exit("variant_id must be chrom:pos:ref:alt")
    v = pd.DataFrame({"chrom": parts[0], "pos": parts[1], "ref": parts[2], "alt": parts[3]})
elif all(c in col for c in ("chrom", "pos", "ref", "alt")):
    v = pd.DataFrame({c: raw[col[c]] for c in ("chrom", "pos", "ref", "alt")})
else:
    sys.exit("the variant table needs columns chrom, pos, ref, alt (or variant_id)")
v = v.apply(lambda s: s.str.strip())
v["input_variant"] = v["chrom"] + ":" + v["pos"] + ":" + v["ref"] + ":" + v["alt"]
v["name"] = raw[col["name"]].str.strip() if "name" in col else ""
v["gene_id"] = raw[col["gene_id"]].str.strip().str.replace(r"\.\d+$", "", regex=True) if "gene_id" in col else ""
v["chrom"] = "chr" + v["chrom"].str.replace(r"^chr", "", regex=True)
v["ref"], v["alt"] = v["ref"].str.upper(), v["alt"].str.upper()
n_in = len(v)

problems = {}
ok = v["chrom"].isin(AUTOSOMES)
problems["not on an autosome"] = (~ok).sum()
valid = v["pos"].str.fullmatch(r"\d+") & v["ref"].str.fullmatch("[ACGT]+") & v["alt"].str.fullmatch("[ACGT]+") \
    & (v["ref"] != v["alt"])
problems["position or alleles not valid (alleles must be A/C/G/T, REF != ALT)"] = (ok & ~valid).sum()
v = v[ok & valid].copy()
v["pos"] = v["pos"].astype(int)

fa = pysam.FastaFile(a.fasta)


def orientation(chrom, pos, ref, alt):
    """'ref' if REF is the genome sequence at pos, 'swap' if ALT is, else 'none'."""
    if fa.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() == ref:
        return "ref"
    return "swap" if fa.fetch(chrom, pos - 1, pos - 1 + len(alt)).upper() == alt else "none"


v["orient"] = [orientation(r.chrom, r.pos, r.ref, r.alt) for r in v.itertuples()]
problems["neither allele matches GRCh38"] = (v["orient"] == "none").sum()
swap = v["orient"] == "swap"
v.loc[swap, ["ref", "alt"]] = v.loc[swap, ["alt", "ref"]].to_numpy()
v = v[v["orient"] != "none"]

v["variant_id"] = v["chrom"] + ":" + v["pos"].astype(str) + ":" + v["ref"] + ":" + v["alt"]
v["SNP"] = np.where(v["name"] != "", v["name"], v["variant_id"])
variants = v.drop_duplicates("variant_id")[["variant_id", "chrom", "pos", "ref", "alt", "SNP", "input_variant"]]
variants = variants.rename(columns={"chrom": "CHR", "pos": "BP", "ref": "REF", "alt": "ALT"})
variants.to_csv(f"{a.out_dir}/variants.tsv", sep="\t", index=False)

# ---- variant-gene pairs: the given genes, else the genes with a TSS within --window bp
genes = pd.read_csv(path("genes_file"), sep="\t")
given = v.loc[v["gene_id"] != "", ["variant_id", "chrom", "gene_id"]].drop_duplicates()
given = given.merge(genes, on="gene_id", how="left", suffixes=("", "_gene"))
unknown = given["gene_TSS"].isna()
other_chrom = given["chrom_gene"].notna() & (given["chrom_gene"] != given["chrom"])
problems["gene_id not in genes.tsv (pair dropped)"] = unknown.sum()
problems["gene_id on another chromosome (pair dropped)"] = other_chrom.sum()
pairs = [given.loc[~unknown & ~other_chrom, ["variant_id", "gene_id", "gene_name"]]]
by_chrom = {c: g.sort_values("gene_TSS") for c, g in genes.dropna(subset=["chrom"]).groupby("chrom")}
for r in variants[~variants["variant_id"].isin(given["variant_id"])].itertuples():
    g = by_chrom[r.CHR]
    tss = g["gene_TSS"].to_numpy()
    lo, hi = np.searchsorted(tss, r.BP - a.window, "left"), np.searchsorted(tss, r.BP + a.window, "right")
    pairs.append(g.iloc[lo:hi][["gene_id", "gene_name"]].assign(variant_id=r.variant_id))
pairs = pd.concat(pairs)[["variant_id", "gene_id", "gene_name"]].drop_duplicates()
pairs.to_csv(f"{a.out_dir}/pairs.tsv", sep="\t", index=False)

print(f"{n_in} input rows -> {len(variants)} variants ({int(swap.sum())} with REF/ALT swapped to GRCh38), "
      f"{len(pairs)} variant-gene pairs, {pairs['gene_id'].nunique()} genes")
for what, n in problems.items():
    if n:
        print(f"  dropped: {n} {what}")
print(f"wrote {a.out_dir}/variants.tsv, {a.out_dir}/pairs.tsv")
