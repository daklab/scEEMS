#!/usr/bin/env python
"""
Fine-mapped eQTL annotation for one chromosome and cell type: the comparator of the predicted-eQTL
annotation in the S-LDSC comparison "fine-mapped eQTLs (PIP > 0.10)" vs "predicted eQTLs (pred_prob > tau*)".

A variant of the baseline model gets the value 1 if any gene fine-maps it with SuSiE PIP > 0.10, and 0
otherwise. The PIP is the FunGen-xQTL fine-mapping PIP carried in the prediction dataset, the same in every
model's dataset; the weighted_full dataset is read. Like make_annotations.py, this takes the maximum over
the genes of the prediction dataset (eQTL and other genes), so both annotations cover the same variants.

usage:   python make_annotations_pip.py CHR COHORT
output:  {sldsc_dir}/pip/{cohort}/MLxQTL_chr{CHR}.annot.gz and MLxQTL_chr{CHR}.l2.M (column {cohort}_pip_binary)
"""
import os
import sys

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

PIP_THRESHOLD = 0.10

chrom = sys.argv[1]
cohort = cohort_name(sys.argv[2])
dask.config.set({"temporary_directory": path("scratch_dir")})
col = f"{cohort}_pip_binary"
out_dir = f"{path('sldsc_dir')}/pip/{cohort}"
os.makedirs(out_dir, exist_ok=True)

# baseline variants (CHR SNP BP A1 A2) of this chromosome
baseline = pd.read_csv(f"{path('baseline_annot_dir')}/baseline_chr{chrom}.annot.gz", sep="\t",
                       usecols=["CHR", "SNP", "BP", "A1", "A2"])

# this chromosome's scored variants, max PIP over genes per variant (BP, A1 = alt, A2 = ref)
pred = dd.read_parquet(f"{path('predictions_parquet_dir', cohort=cohort)}/weighted_full/predictions.parquet")
pred = pred[pred["chr"] == f"chr{chrom}"].compute()
pred = pred.rename(columns={"pos": "BP", "ref": "A2", "alt": "A1"})[["BP", "A1", "A2", "pip"]]
pred = pred.groupby(["BP", "A1", "A2"], as_index=False)["pip"].max()

top = pred.loc[pred["pip"] > PIP_THRESHOLD, ["BP", "A1", "A2"]].copy()
top[col] = 1
annotation_df = baseline.merge(top, on=["BP", "A1", "A2"], how="left")
annotation_df[col] = annotation_df[col].fillna(0).astype(int)

annotation_df.to_csv(f"{out_dir}/MLxQTL_chr{chrom}.annot.gz", sep="\t", index=False, compression="gzip")
sums = annotation_df.drop(columns=["CHR", "SNP", "BP", "A1", "A2"]).to_numpy().sum(axis=0).reshape(1, -1)
np.savetxt(f"{out_dir}/MLxQTL_chr{chrom}.l2.M", sums, fmt="%f", delimiter=" ")
print(f"[{cohort} PIP > {PIP_THRESHOLD} chr{chrom}] {int(sums.sum())} of {len(annotation_df)} baseline "
      f"variants annotated -> {out_dir}/MLxQTL_chr{chrom}.annot.gz")
