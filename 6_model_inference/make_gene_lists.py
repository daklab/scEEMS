#!/usr/bin/env python
"""
List every gene to score for each cell type: the genes of list_genes.csv (MEGA eQTL genes) and
list_genes_other.csv (other genes), both written by step 3 and both featurized in all_variants/.
Prints the number of genes, which is the size of the SLURM array for run_inference.sh.

usage:   python make_gene_lists.py [COHORT ...]          (default: all six cell types)
output:  {gene_list_file}, tab-separated gene_id, chr
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import COHORTS, cohort_name, path

print(f"{'cohort':>16}  {'genes':>6}   array")
for coh in [cohort_name(c) for c in sys.argv[1:]] or COHORTS:
    src = path("gene_list_dir", cohort=coh)
    parts = [pd.read_csv(f"{src}/{fn}", sep="\t")[["gene_id", "chr"]]
             for fn in ("list_genes.csv", "list_genes_other.csv") if os.path.exists(f"{src}/{fn}")]
    if not parts:
        print(f"{coh:>16}  no list_genes.csv / list_genes_other.csv in {src}")
        continue
    genes = pd.concat(parts, ignore_index=True).drop_duplicates("gene_id").reset_index(drop=True)
    dest = path("gene_list_file", cohort=coh)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    genes.to_csv(dest, sep="\t", index=False)
    print(f"{coh:>16}  {len(genes):>6}   --array=1-{len(genes)}%200")
