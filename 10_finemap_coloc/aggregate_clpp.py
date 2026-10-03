#!/usr/bin/env python
"""
Combine the per-chromosome CLPP tables of compute_clpp.R and compute the gene-level CLPP of the EMS paper
(Wang et al. 2021): the maximum over cell types and variants of PIP_eQTL x PIP_GWAS; a gene is prioritized
at CLPP > 0.1.

usage:   python aggregate_clpp.py
outputs: {aggregate_dir}/clpp/clpp_all.tsv and clpp_ems_gene_level.tsv
"""
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

d = os.path.join(path("aggregate_dir"), "clpp")
files = sorted(glob.glob(f"{d}/clpp_*_chr*.tsv"))
if not files:
    sys.exit(f"no per-chromosome CLPP tables in {d}: run compute_clpp.R first")
df = pd.concat((pd.read_csv(f, sep="\t") for f in files), ignore_index=True)
df.to_csv(f"{d}/clpp_all.tsv", sep="\t", index=False, float_format="%.6g")
gene = (df.groupby(["gwas_id", "gwas_prior", "eqtl_prior", "gene_id"])["max_product"]
          .max().rename("clpp_ems").reset_index())
gene.to_csv(f"{d}/clpp_ems_gene_level.tsv", sep="\t", index=False, float_format="%.6g")
print(f"{len(df):,} rows from {len(files)} chromosomes; genes with CLPP > 0.1 by prior:")
print(pd.crosstab(gene.loc[gene.clpp_ems > 0.1, "gwas_prior"], gene.loc[gene.clpp_ems > 0.1, "eqtl_prior"]))
