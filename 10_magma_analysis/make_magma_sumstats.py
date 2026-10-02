#!/usr/bin/env python
"""
Write the GWAS summary statistics in MAGMA's --pval format (space-separated SNP CHR BP P N):

  bellenguez_MAGMA_sumstats.txt   European AD GWAS (Bellenguez et al. 2022), from the munged parquet
                                  (config sumstats_file, the same file S-LDSC uses in step 9)
  ADGC_{AFR,AMR,EAS}_MAGMA_sumstats.txt
                                  ADGC African American, Hispanic/Latino and East Asian AD GWAS
                                  (APOE-adjusted, common variants), matched to the ADSP reference panel
                                  of the same population on chr:pos:ref:alt so that SNP carries the
                                  panel's variant ID; N is the study sample size

The European file is also the SNP location file of the positional gene annotation
(make_positional_annotation.py).

usage:   python make_magma_sumstats.py
outputs: {magma_dir}/sumstats/{bellenguez,ADGC_AFR,ADGC_AMR,ADGC_EAS}_MAGMA_sumstats.txt
"""
import os
import sys

import dask.dataframe as dd
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

# ADGC summary statistics file and sample size per population
ADGC = {"AFR": ("AFA_common_apoe_adj_p-valueOnly.txt", 6728),
        "AMR": ("HISP_common_apoe_adj_p-valueOnly.txt", 8899),
        "EAS": ("EAS_common_apoe_adj_p-valueOnly.txt", 3232)}

out_dir = f"{path('magma_dir')}/sumstats"
os.makedirs(out_dir, exist_ok=True)

# ---- European: Bellenguez ----
eur = pd.read_parquet(path("sumstats_file"), engine="pyarrow")
eur[["SNP", "CHR", "BP", "P", "N"]].to_csv(f"{out_dir}/bellenguez_MAGMA_sumstats.txt", sep=" ", header=True,
                                          index=False)
print(f"bellenguez: {len(eur):,} variants", flush=True)

# ---- ADGC AFR / AMR / EAS ----
for pop, (fname, n) in ADGC.items():
    bim = dd.read_csv(f"{path('plink_ref_dir')}/{pop}/plink/ADSP_{pop}_chr*.bim", sep="\t", header=None,
                      names=["CHR", "SNP", "CM", "BP", "ALT", "REF"]).drop(columns=["CM"])
    df = pd.read_csv(f"{path('adgc_sumstats_dir')}/{fname}", sep="\t")
    df[["CHR", "BP", "REF", "ALT"]] = df["MarkerName"].str.split(":", expand=True)
    df = df.rename(columns={"P-value": "P"})
    df["CHR"] = df["CHR"].apply(lambda x: x.replace("chr", ""))
    df = df.astype({"CHR": int, "BP": int})
    df = df[["CHR", "BP", "P", "REF", "ALT"]].sort_values(by=["CHR", "BP"]).reset_index(drop=True)
    df = df.merge(bim.compute(), how="inner", on=["CHR", "BP", "REF", "ALT"])
    df["N"] = n
    df = df[["SNP", "CHR", "BP", "P", "N"]].drop_duplicates(subset=["SNP"], keep="first").reset_index(drop=True)
    df.to_csv(f"{out_dir}/ADGC_{pop}_MAGMA_sumstats.txt", sep=" ", header=True, index=False)
    print(f"ADGC {pop}: {len(df):,} variants matched to the ADSP {pop} panel", flush=True)
