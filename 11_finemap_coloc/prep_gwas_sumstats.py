#!/usr/bin/env python
"""
Convert the AD GWAS of Bellenguez et al. (2022), stage 1 (GWAS Catalog GCST90027158, GRCh38), into the
sorted, bgzipped, tabix-indexed table that precompute_ld.R, build_coloc_table.R and build_credset_table.R
query by region.

Output columns: CHR POS REF ALT Z N MAF SNP, with Z = BETA / standard_error and rows without a Z removed.
CHR has no "chr" prefix. Alleles are kept as in the source; precompute_ld.R orients Z to the LD panel.

usage:   python prep_gwas_sumstats.py
input:   {gwas_sumstats_raw_file}   (GCST90027158_buildGRCh38.tsv.gz)
output:  {gwas_sumstats_file} (+ .tbi); needs bgzip and tabix (htslib) on PATH
"""
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

src = path("gwas_sumstats_raw_file")
out_gz = path("gwas_sumstats_file")
assert out_gz.endswith(".gz"), "gwas_sumstats_file must end in .gz"
out_tsv = out_gz[:-3]
os.makedirs(os.path.dirname(out_gz), exist_ok=True)

parts = []
for chunk in pd.read_csv(src, sep="\t", chunksize=2_000_000, low_memory=False,
                         usecols=["SNP", "CHR", "POS", "ALT", "REF", "MAF", "BETA", "standard_error", "N"]):
    chunk["Z"] = chunk["BETA"] / chunk["standard_error"]
    parts.append(chunk[["CHR", "POS", "REF", "ALT", "Z", "N", "MAF", "SNP"]])
df = pd.concat(parts, ignore_index=True)
del parts
n = len(df)
df = df.dropna(subset=["Z"])
print(f"{n:,} rows, {n - len(df):,} without Z removed", flush=True)

df = df.sort_values(["CHR", "POS"], kind="mergesort").reset_index(drop=True)
df.to_csv(out_tsv, sep="\t", index=False)
subprocess.run(["bgzip", "-f", out_tsv], check=True)
subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "2", "-S", "1", "-f", out_gz], check=True)
print(f"done: {out_gz} (+ .tbi)")
