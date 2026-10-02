#!/usr/bin/env python
"""
Optional, step 1 of 3: GPN-STAR scores for the non-European MAGMA panel. The GPN-STAR scores of the data
release cover the main variant set only, so the non-European panel's SNVs are scored separately before
score_noneur.py. This writes the per-chromosome SNV inputs of the GPN-STAR scoring (step 2).

GPN-STAR scores a single-position log-likelihood ratio, so only SNVs (one A/C/G/T base to another) are
scored; insertions, deletions, MNVs and spanning deletions (ALT "*") are listed in
unscored_variants.parquet and get no GPN-STAR feature. Each unique variant is scored once; variant_id
(chrN:pos:ref:alt) matches the GPN-STAR file of the data release.

usage:   python aggregate_noneur_variants.py
outputs: {magma_dir}/noneur_gpn_star/snv_input/chr{N}.parquet (variant_id, chrom, pos, ref, alt),
         unscored_variants.parquet, MANIFEST.txt
"""
import os
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared"))
from config import path

SRC = path("noneur_variants_file")                 # CHR ('chr10'), BP, REF, ALT
OUT = f"{path('magma_dir')}/noneur_gpn_star"
ACGT = {"A", "C", "G", "T"}
CHRS = list(range(1, 23))


def classify(df):
    """Add a vtype column; SNV = single ACGT ref -> single ACGT alt (scoreable by GPN-STAR)."""
    ref, alt = df["ref"], df["alt"]
    rl, al = ref.str.len(), alt.str.len()
    is_star = alt.eq("*")
    is_snv = ref.isin(ACGT) & alt.isin(ACGT)
    df["vtype"] = np.select(
        [is_snv, is_star, rl.eq(al), al.gt(rl)],
        ["SNV", "spanning_del", "MNV", "insertion"], default="deletion")
    return df


def main():
    os.makedirs(f"{OUT}/snv_input", exist_ok=True)
    src = pd.read_csv(SRC, sep="\t", usecols=["CHR", "BP", "REF", "ALT"], dtype={"CHR": str})
    src = src.rename(columns={"BP": "pos", "REF": "ref", "ALT": "alt"})
    src["chrom"] = src["CHR"].str.replace("chr", "", regex=False)
    src["variant_id"] = src["CHR"] + ":" + src["pos"].astype(str) + ":" + src["ref"] + ":" + src["alt"]
    src = src.drop_duplicates("variant_id").reset_index(drop=True)

    counts = {v: 0 for v in ("SNV", "MNV", "insertion", "deletion", "spanning_del")}
    per_chrom_snv, total = {}, 0
    uw = None
    for i in CHRS:
        df = classify(src[src["chrom"] == str(i)].copy())
        total += len(df)
        for v, c in df["vtype"].value_counts().items():
            counts[v] = counts.get(v, 0) + int(c)
        snv = (df.loc[df["vtype"] == "SNV", ["variant_id", "chrom", "pos", "ref", "alt"]]
                 .sort_values("pos").reset_index(drop=True))
        snv.to_parquet(f"{OUT}/snv_input/chr{i}.parquet", index=False)
        per_chrom_snv[i] = len(snv)
        un = df.loc[df["vtype"] != "SNV", ["variant_id", "chrom", "pos", "ref", "alt", "vtype"]]
        if len(un):
            tbl = pa.Table.from_pandas(un, preserve_index=False)
            if uw is None:
                uw = pq.ParquetWriter(f"{OUT}/unscored_variants.parquet", tbl.schema)
            uw.write_table(tbl)
        print(f"chr{i:<2}  total={len(df):>8,}  SNV={len(snv):>8,}  unscored={len(df)-len(snv):>7,}", flush=True)
    if uw is not None:
        uw.close()

    n_snv = sum(per_chrom_snv.values())
    with open(f"{OUT}/MANIFEST.txt", "w") as f:
        f.write("GPN-STAR inputs of the non-European MAGMA panel\n")
        f.write(f"total unique variants: {total:,}\n")
        for v in ("SNV", "MNV", "insertion", "deletion", "spanning_del"):
            f.write(f"  {v:13} {counts[v]:>10,} ({100 * counts[v] / total:.1f}%)\n")
        f.write(f"\nSNVs per chromosome in snv_input/:\n")
        for i in CHRS:
            f.write(f"  chr{i}: {per_chrom_snv[i]:,}\n")
    print(f"\ntotal {total:,} | SNV {n_snv:,} | unscored {total - n_snv:,} -> {OUT}")


if __name__ == "__main__":
    main()
