#!/usr/bin/env python
"""
Rebuild the scEEMS (weighted_full) prediction parquet dataset of one cell type from the tabix-indexed
TSVs of the data release (predictions/{cohort}/predictions_{cohort}_{N}.tsv.gz), so that the scEEMS
analyses of step 9, which read the parquet dataset, can be run without rerunning steps 5-7.
The values are identical to those of step 7: the TSVs store every float at full precision and are
parsed back exactly.

usage:   python import_release_predictions.py COHORT
output:  {predictions_parquet_dir}/weighted_full/predictions.parquet/chr=chr{N}/part.0.parquet
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

cohort = cohort_name(sys.argv[1])
src = path("release_predictions_dir", cohort=cohort)
dest = f"{path('predictions_parquet_dir', cohort=cohort)}/weighted_full/predictions.parquet"

total = 0
for chrom in range(1, 23):
    f = f"{src}/predictions_{cohort}_{chrom}.tsv.gz"
    if not os.path.exists(f):
        print(f"  missing {f}", flush=True)
        continue
    # round_trip: the default fast float parser can be off by one unit in the last place
    t = pd.read_csv(f, sep="\t", dtype={"#CHROM": str}, float_precision="round_trip")
    t = t.rename(columns={"POS": "pos", "ID": "variant_id", "REF": "ref", "ALT": "alt", "GENE_ID": "gene_id",
                          "PIP": "pip", "PRED_PROBABILITY": "pred_prob"})
    t = t[["variant_id", "pos", "ref", "alt", "pip", "gene_id", "pred_prob"]]
    os.makedirs(f"{dest}/chr=chr{chrom}", exist_ok=True)
    t.to_parquet(f"{dest}/chr=chr{chrom}/part.0.parquet", index=False)
    total += len(t)
    print(f"  chr{chrom}: {len(t):,}", flush=True)
print(f"[{cohort}] {total:,} predictions -> {dest}")
