#!/usr/bin/env python
"""
Export the eQTL credible sets in the format of the fine_mapping/ folder of the data release: one
bgzipped, tabix-indexed TSV per prior x cell type x chromosome, laid out like the prediction files
(CHROM without "chr", POS, ID, REF, ALT first; header line starts with "#"). Only variants inside a 95%
credible set are exported. All 22 chromosomes are written; a chromosome without a credible set gets a
header-only file.

Columns: CHROM POS ID REF ALT GENE_ID CS_INDEX CS_SIZE CS_PURITY CS_COVERAGE CS_LBF PIP ALPHA VARIANT_LBF
PRIOR_WEIGHT (see the README of the data release). GENE_ID + CS_INDEX (the SuSiE effect index) identify a
credible set.

usage:   python extract_credible_sets_tsv.py [--priors uniform EMS ...] [--cells Mic Ast ...]
input:   {aggregate_dir}/finemapping_comparison/{cell}_{prior}_cs_variants.tsv and _cs.tsv (aggregate_finemap.R)
output:  {release_export_dir}/fine_mapping/{prior}/{cell}_mega_eQTL/credible_sets_{cell}_mega_eQTL_{prior}_{chrom}.tsv.gz
         (+ .tbi) and {release_export_dir}/fine_mapping/credible_set_counts.tsv
"""
import argparse
import os
import sys

import pandas as pd
import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

AGG = os.path.join(path("aggregate_dir"), "finemapping_comparison")
OUT_ROOT = os.path.join(path("release_export_dir"), "fine_mapping")
PRIORS = ["uniform", "EMS", "scEEMS_Weighted_Full", "scEEMS_Unweighted_Full", "scEEMS_Weighted_Restricted"]
CELLS = ["Ast", "Exc", "Inh", "Mic", "Oli", "OPC"]
CHROMS = [str(c) for c in range(1, 23)]
COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "GENE_ID", "CS_INDEX", "CS_SIZE", "CS_PURITY",
           "CS_COVERAGE", "CS_LBF", "PIP", "ALPHA", "VARIANT_LBF", "PRIOR_WEIGHT"]


def build_table(cell, prior):
    """Join each variant row to its credible-set row and return the release columns."""
    var = pd.read_csv(f"{AGG}/{cell}_{prior}_cs_variants.tsv", sep="\t",
                      dtype={"gene_id": str, "chr": str, "variant_id": str})
    cs = pd.read_csv(f"{AGG}/{cell}_{prior}_cs.tsv", sep="\t", dtype={"gene_id": str, "chr": str})
    cs = cs[["gene_id", "cs_idx", "cs_size", "min_abs_corr", "coverage", "lbf"]]
    assert not cs.duplicated(["gene_id", "cs_idx"]).any(), "credible set keys are not unique"
    df = var.merge(cs, on=["gene_id", "cs_idx"], how="left", validate="many_to_one")
    assert df["cs_size"].notna().all(), "variant rows without a matching credible set"

    parts = df["variant_id"].str.split(":", expand=True)
    assert parts.shape[1] == 4 and parts.notna().all().all(), "variant IDs are not chr:pos:ref:alt"
    assert (parts[0] == df["chr"]).all(), "variant chromosome differs from gene chromosome"
    out = pd.DataFrame({
        "CHROM": parts[0].str.replace("chr", "", regex=False),
        "POS": parts[1].astype(int),
        "ID": df["variant_id"],
        "REF": parts[2],
        "ALT": parts[3],
        "GENE_ID": df["gene_id"],
        "CS_INDEX": df["cs_idx"].astype(int),
        "CS_SIZE": df["cs_size"].astype(int),
        "CS_PURITY": df["min_abs_corr"],
        "CS_COVERAGE": df["coverage"],
        "CS_LBF": df["lbf"],
        "PIP": df["pip"],
        "ALPHA": df["alpha"],
        "VARIANT_LBF": df["lbf_variable"],
        "PRIOR_WEIGHT": df["prior_weight"],
    })[COLUMNS]
    sizes = out.groupby(["GENE_ID", "CS_INDEX"]).agg(n=("ID", "size"), size=("CS_SIZE", "first"))
    assert (sizes["n"] == sizes["size"]).all(), "credible sets with missing variants"
    assert len(sizes) == len(cs), "credible sets lost in the join"
    return out.sort_values(["POS", "GENE_ID", "CS_INDEX", "ID"], kind="mergesort")


def write_tabix(df_chr, path_tsv):
    with open(path_tsv, "w") as f:
        f.write("#" + "\t".join(COLUMNS) + "\n")
        df_chr.to_csv(f, sep="\t", index=False, header=False, float_format="%.6g")
    return pysam.tabix_index(path_tsv, seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)


def main():
    ap = argparse.ArgumentParser(description="Export credible sets as tabix-indexed TSVs")
    ap.add_argument("--priors", nargs="+", default=PRIORS, choices=PRIORS)
    ap.add_argument("--cells", nargs="+", default=CELLS, choices=CELLS)
    args = ap.parse_args()

    counts = []
    for prior in args.priors:
        for cell in args.cells:
            cohort = f"{cell}_mega_eQTL"
            out_dir = f"{OUT_ROOT}/{prior}/{cohort}"
            os.makedirs(out_dir, exist_ok=True)
            df = build_table(cell, prior)
            assert set(df["CHROM"]) <= set(CHROMS), "unexpected chromosome"
            written = 0
            for chrom in CHROMS:
                df_chr = df[df["CHROM"] == chrom]
                p = write_tabix(df_chr, f"{out_dir}/credible_sets_{cohort}_{prior}_{chrom}.tsv")
                written += sum(1 for _ in pysam.TabixFile(p).fetch(chrom)) if len(df_chr) else 0
            assert written == len(df), f"{prior}/{cell}: wrote {written} rows, expected {len(df)}"
            n_cs = df[["GENE_ID", "CS_INDEX"]].drop_duplicates().shape[0]
            counts.append({"prior": prior, "cell_type": cohort, "genes_with_cs": df["GENE_ID"].nunique(),
                           "credible_sets": n_cs, "rows": len(df), "unique_variants": df["ID"].nunique()})
            print(f"{prior:28s} {cohort:15s} genes={counts[-1]['genes_with_cs']:6d} "
                  f"credible_sets={n_cs:6d} rows={len(df):8d}", flush=True)

    summary = pd.DataFrame(counts)
    summary_path = f"{OUT_ROOT}/credible_set_counts.tsv"
    if os.path.exists(summary_path) and (set(args.priors) != set(PRIORS) or set(args.cells) != set(CELLS)):
        old = pd.read_csv(summary_path, sep="\t")         # partial run: replace only the regenerated rows
        keep = ~(old["prior"].isin(args.priors) & old["cell_type"].isin([f"{c}_mega_eQTL" for c in args.cells]))
        summary = pd.concat([old[keep], summary], ignore_index=True)
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
