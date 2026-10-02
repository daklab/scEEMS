#!/usr/bin/env python
"""
LD scores of one annotation set for one chromosome, with PolyFun's compute_ldscores.py (1 Mb window, on
the PLINK reference panel). LD scores depend on the annotation values, so they have to be recomputed
whenever an annotation changes. All columns of a set are computed together, once; ldscore_regression.py
then scores each column separately.

Annotation sets (directories under {sldsc_dir}): pareto/{cohort}_{model}, pip/{cohort}, top{N}/{cohort},
cs/{cohort}.

usage:   python compute_ldscores.py SET CHR          e.g. python compute_ldscores.py pip/Mic_mega_eQTL 22
output:  {sldsc_dir}/{SET}/MLxQTL_chr{CHR}.l2.ldscore.parquet
env:     the PolyFun environment (paths.ldsc_dir is the PolyFun directory)
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

annot_set, chrom = sys.argv[1], sys.argv[2]
ann = f"{path('sldsc_dir')}/{annot_set}"

cmd = [sys.executable, "compute_ldscores.py",
       "--annot", f"{ann}/MLxQTL_chr{chrom}.annot.gz",
       "--bfile", f"{path('ldsc_bfile_prefix')}{chrom}",
       "--out", f"{ann}/MLxQTL_chr{chrom}.l2.ldscore.parquet",
       "--allow-missing",
       "--ld-wind-kb", "1000"]
print(" ".join(cmd), flush=True)
sys.exit(subprocess.call(cmd, cwd=path("ldsc_dir")))
