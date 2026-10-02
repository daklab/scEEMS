#!/usr/bin/env python
"""
Optional, step 3 of 3: combine the per-chromosome GPN-STAR scores of the non-European panel's SNVs
(run_score_gpn_star.sh) into the file score_noneur.py reads (config noneur_gpn_star_file).

usage:   python merge_scores_noneur.py
output:  {noneur_gpn_star_file}: variant_id, gpn_star_llr (SNVs only)
"""
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "shared"))
from config import path

SC = f"{path('magma_dir')}/noneur_gpn_star/snv_scored"
rows = []
for sc in sorted(glob.glob(f"{SC}/chr*.parquet"), key=lambda p: int(os.path.basename(p)[3:-8])):
    part = pd.read_parquet(sc, columns=["variant_id", "gpn_star_llr"])
    rows.append(part)
    print(f"{os.path.basename(sc)}: {len(part):,} variants")

out = pd.concat(rows, ignore_index=True).drop_duplicates("variant_id")
dest = path("noneur_gpn_star_file")
os.makedirs(os.path.dirname(dest), exist_ok=True)
out.to_parquet(dest, index=False)
print(f"\n{dest}: {len(out):,} SNVs scored")
