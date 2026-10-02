#!/usr/bin/env python
"""
Build the positional MAGMA gene annotation, the baseline the eMAGMA results are compared with: a SNP is
assigned to a gene only if it lies within the gene body (MAGMA's default window of 0 kb up- and
downstream). The annotation uses only gene coordinates (NCBI build 38 gene locations distributed with
MAGMA) and the positions of the European summary-statistics SNPs, so it does not depend on the models.

usage:   python make_positional_annotation.py           (after make_magma_sumstats.py)
output:  {magma_dir}/positional/bellenguez_MAGMA.genes.annot
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

out_dir = f"{path('magma_dir')}/positional"
os.makedirs(out_dir, exist_ok=True)
cmd = [path("magma_binary"), "--annotate",
       "--snp-loc", f"{path('magma_dir')}/sumstats/bellenguez_MAGMA_sumstats.txt",
       "--gene-loc", path("magma_gene_loc_file"),
       "--out", f"{out_dir}/bellenguez_MAGMA"]                 # MAGMA appends .genes.annot
print(" ".join(cmd), flush=True)
subprocess.run(cmd, check=True)
