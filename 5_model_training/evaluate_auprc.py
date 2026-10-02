#!/usr/bin/env python
"""
Held-out AUPRC of the three models for one cell type, from the LOCO test predictions of train_loco.py.

The test predictions of all held-out chromosomes are pooled. Reported per model: pooled AUPRC (all
variants and SNVs only) and the mean of the per-chromosome AUPRCs; and, for Weighted (Full) against each
comparison model, the paired bootstrap 95% confidence interval of the AUPRC difference (rows resampled
with replacement, both models scored on the same resample).

usage:   python evaluate_auprc.py COHORT [N_BOOT=2000]
output:  {test_predictions_dir}/auprc.json
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

SEED = 9448
MODELS = ["weighted_full", "unweighted_full", "weighted_restricted"]
COMPARISONS = [("weighted_full", "unweighted_full"), ("weighted_full", "weighted_restricted")]

cohort = cohort_name(sys.argv[1])
NB = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
TESTPRED = path("test_predictions_dir", cohort=cohort)
files = sorted(glob.glob(f"{TESTPRED}/test_pred_chr*.parquet"))
if not files:
    sys.exit(f"no test predictions under {TESTPRED}: run train_loco.py first")
df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
y = df["y_true"].to_numpy().astype(int)
snv = df["is_snv"].to_numpy().astype(bool)


def ap(y_, p_, mask=None):
    """Average precision over rows with a prediction (and in mask); NaN if only one class is present."""
    m = ~np.isnan(p_)
    if mask is not None:
        m = m & mask
    yy, pp = y_[m], p_[m]
    return float(average_precision_score(yy, pp)) if len(np.unique(yy)) > 1 else float("nan")


pooled = {m: ap(y, df[f"pred_{m}"].to_numpy()) for m in MODELS}
pooled_snv = {m: ap(y, df[f"pred_{m}"].to_numpy(), snv) for m in MODELS}
mean_per_chr = {m: float(np.nanmean([ap(g["y_true"].to_numpy().astype(int), g[f"pred_{m}"].to_numpy())
                                     for _, g in df.groupby("chrom")])) for m in MODELS}

rng = np.random.default_rng(SEED)


def boot_ci(a, b, snv_only=False):
    """Paired bootstrap 95% CI of AUPRC(a) - AUPRC(b) over rows both models scored."""
    pa, pb = df[f"pred_{a}"].to_numpy(), df[f"pred_{b}"].to_numpy()
    m = ~np.isnan(pa) & ~np.isnan(pb)
    if snv_only:
        m = m & snv
    idx = np.flatnonzero(m)
    deltas = []
    for _ in range(NB):
        s = rng.choice(idx, size=idx.size, replace=True)
        ys = y[s]
        if len(np.unique(ys)) < 2:
            continue
        deltas.append(average_precision_score(ys, pa[s]) - average_precision_score(ys, pb[s]))
    deltas = np.asarray(deltas)
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return dict(delta=float(deltas.mean()), ci=[float(lo), float(hi)], significant=bool(lo > 0 or hi < 0),
                n=int(idx.size))


cis = {f"{a}_vs_{b}": dict(overall=boot_ci(a, b, False), snv=boot_ci(a, b, True)) for a, b in COMPARISONS}
out = dict(cohort=cohort, n_chromosomes=len(files), n_test=int(len(df)), n_pos=int(y.sum()),
           n_snv=int(snv.sum()), n_boot=NB, pooled_auprc=pooled, pooled_auprc_snv=pooled_snv,
           mean_per_chr_auprc=mean_per_chr, bootstrap_ci=cis)
with open(f"{TESTPRED}/auprc.json", "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2))
