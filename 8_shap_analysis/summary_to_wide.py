#!/usr/bin/env python
"""
Reshape the SHAP summary (aggregate_shap.py) to one row per cell type and region with one
{category}_shap_sum column per feature category (the mean summed |SHAP|), the table behind the
manuscript's SHAP figure.

usage:   python summary_to_wide.py [SUMMARY_TSV]     (default: {aggregate_dir}/shap/all_cohorts_shap_summary.tsv)
output:  the same path with _wide.tsv
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

RENAME = {"Enformer": "enformer_shap_sum", "chromBPNet": "chrombpnet_shap_sum",
          "ABC": "abc_score_shap_sum", "CRE": "brain_CRE_shap_sum",
          "distance": "distance_shap_sum", "variant_type": "variant_shap_sum",
          "gene_lof": "gene_shap_sum", "baseline": "baseline_shap_sum",
          "TF": "all_TF_shap_sum", "gnomad_MAF": "gnomad_MAF_shap_sum",
          "abs_gpn": "abs_gpn_shap_sum"}
WIDE_COLS = ["genomic_region", "enformer_shap_sum", "chrombpnet_shap_sum", "abc_score_shap_sum",
             "brain_CRE_shap_sum", "distance_shap_sum", "variant_shap_sum", "gene_shap_sum",
             "baseline_shap_sum", "all_TF_shap_sum", "gnomad_MAF_shap_sum", "abs_gpn_shap_sum", "cohort"]

src = sys.argv[1] if len(sys.argv) > 1 else f"{path('aggregate_dir')}/shap/all_cohorts_shap_summary.tsv"
long = pd.read_csv(src, sep="\t")
wide = (long.pivot(index=["cohort", "region"], columns="category", values="mean_abs_shap")
        .rename(columns=RENAME).reset_index().rename(columns={"region": "genomic_region"}))
wide.columns.name = None
for c in WIDE_COLS:
    if c not in wide.columns:
        wide[c] = 0.0
wide = wide[WIDE_COLS]
dest = src.replace(".tsv", "_wide.tsv")
wide.to_csv(dest, sep="\t", index=False)
print(f"{src} -> {dest}: {len(wide)} rows")
