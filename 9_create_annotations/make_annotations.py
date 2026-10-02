#!/usr/bin/env python
"""
Predicted-eQTL annotations for the S-LDSC threshold sweep, for one chromosome, cell type and model.

A variant of the baseline model gets the value 1 at threshold p if the model predicts it to be a causal
eQTL with pred_prob > p for any gene, and 0 otherwise. Twenty thresholds, p = 0.80, 0.81, ..., 0.99, are
written as the twenty columns {cohort}_{model}_pred_prob_{p} of one annotation file. Each column is then
scored against the baseline on its own (ldscore_regression.py), and tau_star.py picks the threshold.

Predictions are joined to the baseline variants on position and alleles (BP, A1 = alternative allele,
A2 = reference allele), as in every annotation of this step.

usage:   python make_annotations.py CHR COHORT MODEL
         MODEL: weighted_full, unweighted_full or weighted_restricted
output:  {sldsc_dir}/pareto/{cohort}_{model}/MLxQTL_chr{CHR}.annot.gz and MLxQTL_chr{CHR}.l2.M
"""
import os
import sys

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

THRESHOLDS = [round(0.80 + 0.01 * i, 2) for i in range(20)]     # 0.80, 0.81, ..., 0.99
MODELS = ["weighted_full", "unweighted_full", "weighted_restricted"]

chrom = sys.argv[1]
cohort = cohort_name(sys.argv[2])
model = sys.argv[3]
assert model in MODELS, f"unknown model {model}"
dask.config.set({"temporary_directory": path("scratch_dir")})
prefix = f"{cohort}_{model}"
out_dir = f"{path('sldsc_dir')}/pareto/{cohort}_{model}"
os.makedirs(out_dir, exist_ok=True)

# baseline variants (CHR SNP BP A1 A2) of this chromosome
baseline = pd.read_csv(f"{path('baseline_annot_dir')}/baseline_chr{chrom}.annot.gz", sep="\t",
                       usecols=["CHR", "SNP", "BP", "A1", "A2"])
annotation_df = baseline.copy(deep=True)

# this chromosome's predictions, max pred_prob over genes per variant (BP, A1 = alt, A2 = ref)
pred = dd.read_parquet(f"{path('predictions_parquet_dir', cohort=cohort)}/{model}/predictions.parquet")
pred = pred[pred["chr"] == f"chr{chrom}"].compute()
pred = pred.rename(columns={"pos": "BP", "ref": "A2", "alt": "A1"})[["BP", "A1", "A2", "pred_prob"]]
pred = pred.groupby(["BP", "A1", "A2"], as_index=False)["pred_prob"].max()

for t in THRESHOLDS:
    col = f"{prefix}_pred_prob_{t}"
    top = pred.loc[pred["pred_prob"] > t, ["BP", "A1", "A2"]].copy()
    top[col] = 1
    m = baseline.merge(top, on=["BP", "A1", "A2"], how="left")
    m[col] = m[col].fillna(0).astype(int)
    annotation_df = annotation_df.merge(m, on=["CHR", "SNP", "BP", "A1", "A2"], how="left")

annotation_df.to_csv(f"{out_dir}/MLxQTL_chr{chrom}.annot.gz", sep="\t", index=False, compression="gzip")
sums = annotation_df.drop(columns=["CHR", "SNP", "BP", "A1", "A2"]).to_numpy().sum(axis=0).reshape(1, -1)
np.savetxt(f"{out_dir}/MLxQTL_chr{chrom}.l2.M", sums, fmt="%f", delimiter=" ")
print(f"[{cohort} {model} chr{chrom}] {len(THRESHOLDS)} thresholds -> {out_dir}/MLxQTL_chr{chrom}.annot.gz "
      f"(annotated variants per threshold {int(sums.min())}..{int(sums.max())})")
