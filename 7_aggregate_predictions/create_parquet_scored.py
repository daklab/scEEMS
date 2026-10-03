#!/usr/bin/env python
"""
Collect one model's per-gene prediction tables (step 6) into a single parquet dataset partitioned by
chromosome, the input of steps 8-11 and of the TSV export (extract_predictions_tsv.py).

usage:   python create_parquet_scored.py COHORT [MODEL]
         MODEL: weighted_full (default), unweighted_full or weighted_restricted
output:  {predictions_parquet_dir}/{model}/predictions.parquet/chr=chr{N}/ with columns
         variant_id, pos, ref, alt, pip, gene_id, pred_prob (chr is the partition key)
"""
import os
import sys

import dask
from dask import dataframe as dd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

cohort = cohort_name(sys.argv[1])
model = sys.argv[2] if len(sys.argv) > 2 else "weighted_full"
assert model in ("weighted_full", "unweighted_full", "weighted_restricted"), f"unknown model {model}"
dask.config.set({"temporary_directory": path("scratch_dir")})

src = f"{path('predictions_dir', cohort=cohort)}/{model}"
dest = f"{path('predictions_parquet_dir', cohort=cohort)}/{model}/predictions.parquet"
os.makedirs(os.path.dirname(dest), exist_ok=True)

df = dd.read_csv(f"{src}/*_predictions.tsv", sep="\t")
df = df.repartition(partition_size="100MB")
df.to_parquet(dest, engine="pyarrow", compression="snappy", write_index=False, partition_on=["chr"],
              overwrite=True)
print(f"[{cohort} {model}] {src} -> {dest}")
