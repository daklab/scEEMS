#!/usr/bin/env python
"""
Summarize the per-chromosome SHAP tables (shap_analysis.py) of all cell types: mean summed |SHAP| per
feature category over the predicted eQTLs, separately for promoter-like and enhancer-like variants, and
each category's share of the total.

usage:   python aggregate_shap.py
output:  {aggregate_dir}/shap/all_cohorts_shap_summary.tsv with columns
         cohort, region, category, mean_abs_shap, share, n_variants
"""
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import COHORTS, path

CAT_ORDER = ["ABC", "CRE", "chromBPNet", "Enformer", "TF", "abs_gpn",
             "gene_lof", "variant_type", "baseline", "distance", "gnomad_MAF"]
OUT_DIR = f"{path('aggregate_dir')}/shap"
os.makedirs(OUT_DIR, exist_ok=True)

rows = []
for coh in COHORTS:
    files = sorted(glob.glob(f"{path('shap_dir', cohort=coh)}/shap_category_chr*.parquet"))
    if not files:
        print(f"  {coh}: no shap_category files -- skipped")
        continue
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    n = df.groupby("region").size().rename("n_variants")
    m = df.groupby("region")[CAT_ORDER].mean()
    long = m.reset_index().melt(id_vars="region", var_name="category", value_name="mean_abs_shap")
    long["share"] = long["mean_abs_shap"] / long.groupby("region")["mean_abs_shap"].transform("sum")
    long = long.merge(n.reset_index(), on="region")
    long["cohort"] = coh
    rows.append(long)
    print(f"  {coh}: {len(files)}/22 chromosomes, {len(df):,} predicted eQTLs "
          f"(enhancer-like {int((df.region == 'enhancer').sum()):,}, promoter-like {int((df.region == 'promoter').sum()):,})")

allr = pd.concat(rows, ignore_index=True)
allr["category"] = pd.Categorical(allr["category"], CAT_ORDER, ordered=True)
allr = allr.sort_values(["cohort", "region", "category"])[
    ["cohort", "region", "category", "mean_abs_shap", "share", "n_variants"]]
allr.to_csv(f"{OUT_DIR}/all_cohorts_shap_summary.tsv", sep="\t", index=False)
print(f"wrote {OUT_DIR}/all_cohorts_shap_summary.tsv: {len(allr)} rows")
