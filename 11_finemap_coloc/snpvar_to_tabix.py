#!/usr/bin/env python
"""
Collapse PolyFun's per-chromosome SNPVAR (compute_snpvar.py) into one sorted, bgzipped, tabix-indexed
table for region queries by precompute_gwas.R. Uses the constrained (non-negative) estimates.

usage:   python snpvar_to_tabix.py
input:   {snpvar_dir}/bellenguez_sldsc83.{1-22}.snpvar_ridge_constrained.gz
output:  {snpvar_dir}/bellenguez_sldsc83.snpvar.hg38.tsv.gz (+ .tbi): CHR POS A1 A2 SNPVAR, no header;
         needs bgzip and tabix (htslib) on PATH
"""
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

d = path("snpvar_dir")
out_tsv = os.path.join(d, "bellenguez_sldsc83.snpvar.hg38.tsv")

parts = []
for chrom in range(1, 23):
    f = os.path.join(d, f"bellenguez_sldsc83.{chrom}.snpvar_ridge_constrained.gz")
    df = pd.read_csv(f, sep="\t", usecols=["CHR", "BP", "A1", "A2", "SNPVAR"])
    if (df["SNPVAR"] < 0).any():
        sys.exit(f"{f}: negative SNPVAR in the constrained estimates")
    parts.append(df.sort_values("BP", kind="mergesort")[["CHR", "BP", "A1", "A2", "SNPVAR"]])
    print(f"  chr{chrom}: {len(df):,}", flush=True)
out = pd.concat(parts, ignore_index=True)
out.to_csv(out_tsv, sep="\t", index=False, header=False)
subprocess.run(["bgzip", "-f", out_tsv], check=True)
subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", out_tsv + ".gz"], check=True)
print(f"done: {out_tsv}.gz (+ .tbi), {len(out):,} variants")
