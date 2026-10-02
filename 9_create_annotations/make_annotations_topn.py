#!/usr/bin/env python
"""
Size-matched S-LDSC, step 2: annotation columns for one chromosome and cell type, one column per ranking
of select_topn.py, each marking that ranking's genome-wide top N variants. Every column holds exactly N
baseline variants genome-wide, so Prop._SNPs is the same for all of them.

The ranking is done genome-wide by select_topn.py, not per chromosome, or "top N" would mean top N per
chromosome.

Columns {cell}_top{N}_{ranking}:
  pred_weighted_full, pred_unweighted_full, pred_weighted_restricted   the three models, eQTL + other genes
  pip              the FunGen-xQTL fine-mapping (uniform prior), eQTL + other genes: the fine-mapped
                   comparator of pred_weighted_full
  pip_{prior}      the five fine-mapping priors of step 11, eQTL genes fine-mapped under all five priors

usage:   python make_annotations_topn.py CHR COHORT [N=5000]
output:  {sldsc_dir}/top{N}/{cohort}/MLxQTL_chr{CHR}.annot.gz and MLxQTL_chr{CHR}.l2.M
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

chrom = sys.argv[1]
cohort = cohort_name(sys.argv[2])
cell = cohort.replace("_mega_eQTL", "")
N = int(sys.argv[3]) if len(sys.argv) > 3 else 5000

SEL = f"{path('sldsc_dir')}/top{N}_selection/{cohort}_top{N}.tsv.gz"
out_dir = f"{path('sldsc_dir')}/top{N}/{cohort}"
os.makedirs(out_dir, exist_ok=True)

baseline = pd.read_csv(f"{path('baseline_annot_dir')}/baseline_chr{chrom}.annot.gz", sep="\t",
                       usecols=["CHR", "SNP", "BP", "A1", "A2"])
sel = pd.read_csv(SEL, sep="\t")
sel = sel[sel["CHR"].astype(str) == str(chrom)]

annot = baseline.copy()
for ranking in sorted(pd.read_csv(SEL, sep="\t", usecols=["ranking"])["ranking"].unique()):
    name = f"{cell}_top{N}_{ranking}"
    keys = sel.loc[sel["ranking"] == ranking, ["BP", "A1", "A2"]].drop_duplicates()
    if len(keys):
        m = baseline.merge(keys.assign(**{name: 1}), on=["BP", "A1", "A2"], how="left")
        m = m[["CHR", "SNP", "BP", "A1", "A2", name]]
    else:
        m = baseline.assign(**{name: np.nan})
    annot = annot.merge(m, on=["CHR", "SNP", "BP", "A1", "A2"], how="left")
    annot[name] = annot[name].fillna(0).astype(int)

annot.to_csv(f"{out_dir}/MLxQTL_chr{chrom}.annot.gz", sep="\t", index=False, compression="gzip")
sums = annot.drop(columns=["CHR", "SNP", "BP", "A1", "A2"]).to_numpy().sum(axis=0).reshape(1, -1)
np.savetxt(f"{out_dir}/MLxQTL_chr{chrom}.l2.M", sums, fmt="%f", delimiter=" ")
cols = [c for c in annot.columns if c not in ("CHR", "SNP", "BP", "A1", "A2")]
print(f"[{cohort} chr{chrom}] {len(annot):,} baseline variants, {len(cols)} columns", flush=True)
for c, s in zip(cols, sums[0]):
    print(f"    {c:<48} {int(s):>7,}", flush=True)
