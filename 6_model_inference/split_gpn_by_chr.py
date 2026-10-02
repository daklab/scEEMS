#!/usr/bin/env python
"""
One-time preparation for inference: split the GPN-STAR scores (gpn_star_scores_all.parquet, 15.3 million
SNVs) into one file per chromosome, so each per-gene job loads only its chromosome's scores.
shared/featurize.load_gpn_map uses these files when they exist and otherwise filters the full file.

usage:   python split_gpn_by_chr.py
output:  {gpn_star_by_chr_dir}/chr{N}.parquet (variant_id, gpn_star_llr)
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
import featurize as F
from config import path

out_dir = path("gpn_star_by_chr_dir")
os.makedirs(out_dir, exist_ok=True)
g = pd.read_parquet(path("gpn_star_file"), columns=["variant_id", F.GPN_COL])
g["chrom"] = g["variant_id"].str.split(":", n=1).str[0]
n = 0
for chrom, sub in g.groupby("chrom"):
    sub[["variant_id", F.GPN_COL]].to_parquet(f"{out_dir}/{chrom}.parquet", index=False)
    n += len(sub)
    print(f"{chrom}: {len(sub):,}", flush=True)
print(f"done: {n:,} scores -> {out_dir}")
