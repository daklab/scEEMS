#!/usr/bin/env python
"""
Standardized effect size tau* with block-jackknife significance for every (cell type, model, pred_prob
threshold) of the S-LDSC threshold sweep, and the selected threshold tau* of each cell type.

tau* is the per-SNP heritability contributed by a one-standard-deviation increase in the annotation, as a
fraction of the total heritability (Gazal et al. 2017):

    tau*_c = tau_c * sd_c * M / h2_g          sd_c = sqrt(P (1 - P)) for a binary annotation (P = Prop._SNPs)

M (the number of regression SNPs) and h2_g are read from each run's ldsc .log. The standard error comes
from LDSC's block jackknife (--print-delete-vals: the coefficient with each of 200 blocks left out),
rescaled to tau* units:

    SE = sd(tau*_blocks) * sqrt(n_blocks - 1)       z = tau* / SE       p = 2 * P(Z > |z|)

tau* is a constant rescaling of the LDSC coefficient, so z and p equal LDSC's own (checked below); M, h2
and sd only set the magnitude, which makes tau* comparable across annotations and cell types.

Selection: tau* of a cell type is the threshold with the largest tau* for the weighted_full model.
Enrichment rises monotonically with the threshold, so it would always pick 0.99; tau* has an interior
maximum. The p-value of the selected threshold is reported, not used as a filter.

usage:   python tau_star.py [COHORT ...]              (default: all six cell types)
outputs: {aggregate_dir}/tau_star_jackknife.tsv       one row per cell type, model and threshold
         {aggregate_dir}/tau_star.json                {cohort: selected threshold}, for every cell type
                                                      whose 20-threshold sweep is complete; read by
                                                      step 8 and by aggregate_prediction_vs_pip.py
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import COHORTS, cohort_name, path

RES = f"{path('sldsc_dir')}/results/pareto"
OUT_DIR = path("aggregate_dir")
COHORT_LIST = [cohort_name(c) for c in sys.argv[1:]] or COHORTS
MODELS = ["weighted_full", "unweighted_full", "weighted_restricted"]
SEL_MODEL = "weighted_full"
N_THRESH = 20

_h2_re = re.compile(r"Total Observed scale h2:\s*([-\d.eE+]+)")
_m_re = re.compile(r"After merging with reference panel LD,\s*([\d]+)\s*SNPs remain")


def h2_and_M(log_path):
    """Total observed-scale h2 and the number of regression SNPs M, from a run's ldsc .log."""
    txt = open(log_path).read()
    h2, m = _h2_re.search(txt), _m_re.search(txt)
    return (float(h2.group(1)) if h2 else np.nan), (int(m.group(1)) if m else np.nan)


rows = []
for coh in COHORT_LIST:
    for model in MODELS:
        for rf in sorted(glob.glob(f"{RES}/{coh}_{model}/*.results")):
            stem = rf[:-len(".results")]
            pdel, log = f"{stem}.part_delete", f"{stem}.log"
            if not (os.path.exists(pdel) and os.path.exists(log)):
                print(f"  {coh} {model} {os.path.basename(stem)}: missing .part_delete or .log -- skipped")
                continue
            res = pd.read_csv(rf, sep="\t")
            cats = res["Category"].tolist()
            target = cats[-1]                                   # the focal annotation comes last
            if not target.startswith(f"{coh}_{model}_pred_prob_"):
                print(f"  {coh} {model}: last Category is {target!r}, not a pred_prob annotation -- skipped")
                continue
            h2, M = h2_and_M(log)
            r = res.iloc[-1]
            sd_annot = np.sqrt(r["Prop._SNPs"] * (1 - r["Prop._SNPs"]))
            k = sd_annot * M / h2                               # tau_c -> tau* rescaling

            dl = pd.read_csv(pdel, sep=" ", header=None, names=cats)
            blocks = dl[target].to_numpy() * k                  # per-block tau*
            n = len(blocks)
            tau_star = blocks.mean()
            se = blocks.std(ddof=1) * np.sqrt(n - 1)            # block-jackknife SE
            z = tau_star / se if se else np.nan
            rows.append(dict(
                cohort=coh, model=model,
                pred_prob=float(re.search(r"_pred_prob_([0-9.]+)", target).group(1)),
                category=target,
                prop_SNPs=r["Prop._SNPs"], prop_h2=r["Prop._h2"], prop_h2_SE=r["Prop._h2_std_error"],
                enrichment=r["Enrichment"], enrichment_SE=r["Enrichment_std_error"],
                enrichment_p=r["Enrichment_p"],
                coefficient=r["Coefficient"], coefficient_SE=r["Coefficient_std_error"],
                coef_z_ldsc=r["Coefficient_z-score"],
                tau_star=tau_star, tau_star_SE=se, tau_star_z=z,
                tau_star_p=2 * norm.sf(abs(z)), h2=h2, M=M, n_blocks=n))

if not rows:
    sys.exit(f"no S-LDSC results under {RES}: run ldscore_regression.py first")
df = pd.DataFrame(rows).sort_values(["cohort", "model", "pred_prob"])
os.makedirs(OUT_DIR, exist_ok=True)
df.to_csv(f"{OUT_DIR}/tau_star_jackknife.tsv", sep="\t", index=False)
df["sig"] = df["tau_star_p"] < 0.05
print(f"wrote {OUT_DIR}/tau_star_jackknife.tsv: {len(df)} rows "
      f"({df.cohort.nunique()} cell types x {df.model.nunique()} models x {df.pred_prob.nunique()} thresholds)\n")

# tau* is a constant rescaling of the LDSC coefficient, so its jackknife z must equal LDSC's z
dz = (df["tau_star_z"] - df["coef_z_ldsc"]).abs()
print(f"check |tau*_z - ldsc_coef_z|: max={dz.max():.4f} median={dz.median():.4f} (should be ~0)\n")

# ------------------------------------------------------------------ tau* selection
sel = {}
print(f"=== tau* selection (argmax tau*, model {SEL_MODEL}) ===")
for coh in COHORT_LIST:
    g = df[(df.model == SEL_MODEL) & (df.cohort == coh)]
    if len(g) < N_THRESH:
        print(f"  {coh:22} {len(g)}/{N_THRESH} thresholds -- incomplete, not selected")
        continue
    r = g.loc[g["tau_star"].idxmax()]
    sel[coh] = round(float(r["pred_prob"]), 2)
    print(f"  {coh:22} tau* = {sel[coh]}   (tau*={r['tau_star']:.3f}, p={r['tau_star_p']:.2g}, "
          f"prop_h2={100 * r['prop_h2']:.1f}%, enrichment={r['enrichment']:.1f}; "
          f"{int(g.sig.sum())}/{len(g)} thresholds significant)"
          f"{'' if r['tau_star_p'] < 0.05 else '   ** not significant: no heritability beyond the baseline **'}")

if sel:
    with open(f"{OUT_DIR}/tau_star.json", "w") as fh:
        json.dump(sel, fh, indent=2)
    print(f"\nwrote {OUT_DIR}/tau_star.json ({len(sel)} cell type{'s' if len(sel) != 1 else ''})")
