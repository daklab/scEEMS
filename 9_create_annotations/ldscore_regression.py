#!/usr/bin/env python
"""
Stratified LD score regression (S-LDSC) of one annotation column against the GWAS, with PolyFun's ldsc.py.

Each run fits the baseline model plus ONE focal column of an annotation set. The columns of a set are
nested thresholds or overlapping variant sets, so fitting them jointly would split one signal across
near-collinear regressors; comparing columns is therefore a ranking over separate models, not a variance
decomposition. Every set is fitted with the same summary statistics, baseline, regression weights and
flags, so all results are comparable. --print-delete-vals writes the block-jackknife values tau_star.py
needs.

usage:   python ldscore_regression.py SET COLUMN
         SET: an annotation set under {sldsc_dir} (pareto/{cohort}_{model}, pip/{cohort}, top{N}/{cohort},
              cs/{cohort}); COLUMN: a column name, or its 1-based position in the set (SLURM array index)
output:  {sldsc_dir}/results/{SET}/{column}.results (+ .log, .part_delete, .delete)
env:     the PolyFun environment (paths.ldsc_dir is the PolyFun directory)
"""
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

META = ["CHR", "SNP", "BP", "A1", "A2"]

annot_set, column = sys.argv[1], sys.argv[2]
BL = f"{path('baseline_annot_dir')}/baseline_chr"
WEIGHT = f"{path('weights_dir')}/weights_chr"
MLX = f"{path('sldsc_dir')}/{annot_set}/MLxQTL_chr"
OUT = f"{path('sldsc_dir')}/results/{annot_set}"
os.makedirs(OUT, exist_ok=True)

bl_cols = [c for c in pd.read_csv(f"{BL}22.annot.gz", sep="\t", nrows=0).columns if c not in META]
mlx_cols = [c for c in pd.read_csv(f"{MLX}22.annot.gz", sep="\t", nrows=0).columns if c not in META]
if column.isdigit():
    column = mlx_cols[int(column) - 1]
assert column in mlx_cols, f"{column} is not one of the {len(mlx_cols)} columns of {MLX}22.annot.gz"

cmd = [sys.executable, "ldsc.py",
       "--h2", path("sumstats_file"),
       "--ref-ld-chr", f"{BL},{MLX}",
       "--w-ld-chr", WEIGHT,
       "--anno", ",".join(bl_cols + [column]),
       "--out", f"{OUT}/{column}",
       "--overlap-annot", "--not-M-5-50", "--print-coefficients", "--print-delete-vals"]
print(" ".join(cmd), flush=True)
sys.exit(subprocess.call(cmd, cwd=path("ldsc_dir")))
