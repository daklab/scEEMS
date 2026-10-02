#!/usr/bin/env python
"""
List the genes model_inference.py has not scored yet, as 1-based rows of the gene list, so they can be
resubmitted (for example with more memory):
    sbatch --export=ALL,cohort=Mic_mega_eQTL --array=$(paste -sd, missing_genes_Mic_mega_eQTL.txt) run_inference.sh
A gene counts as scored when its weighted_full and unweighted_full prediction tables both exist.

usage:   python find_missing_genes.py COHORT
output:  missing_genes_{cohort}.txt in the current directory
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

cohort = cohort_name(sys.argv[1])
genes = pd.read_csv(path("gene_list_file", cohort=cohort), sep="\t")
pred_dir = path("predictions_dir", cohort=cohort)
missing = [str(i) for i, g in enumerate(genes["gene_id"], start=1)
           if not all(os.path.exists(f"{pred_dir}/{m}/{g}_predictions.tsv") for m in ("weighted_full", "unweighted_full"))]
with open(f"missing_genes_{cohort}.txt", "w") as fh:
    fh.write("\n".join(missing) + ("\n" if missing else ""))
print(f"{cohort}: {len(missing)} of {len(genes)} genes not scored -> missing_genes_{cohort}.txt")
