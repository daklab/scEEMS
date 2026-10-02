#!/usr/bin/env python
"""
Collect the two arms of the heritability comparison "predicted eQTLs vs fine-mapped eQTLs" into one table:

  prediction arm   weighted_full predictions with pred_prob > tau* (tau_star.json), taken from the
                   threshold sweep (tau_star_jackknife.tsv); no extra S-LDSC run is needed
  pip arm          fine-mapped eQTLs, PIP > 0.10 (make_annotations_pip.py)

Both arms come from the same S-LDSC setup (summary statistics, baseline, weights, flags) and the same
variants (maximum over the genes of the prediction dataset), so they are directly comparable.

usage:   python aggregate_prediction_vs_pip.py
output:  {aggregate_dir}/prediction_vs_pip_sldsc.tsv with columns
         Category ({cohort}_pred_prob_binary_1 or {cohort}_pip_binary_1), Prop._SNPs, Prop._h2,
         Prop._h2_std_error, Enrichment, Enrichment_std_error, Enrichment_p, Coefficient,
         Coefficient_std_error, Coefficient_z-score, cohort, arm (prediction | pip), model, threshold
"""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

AGG = path("aggregate_dir")
PIP_RES = f"{path('sldsc_dir')}/results/pip"
OUT = f"{AGG}/prediction_vs_pip_sldsc.tsv"
SEL_MODEL = "weighted_full"
PUB = ["Category", "Prop._SNPs", "Prop._h2", "Prop._h2_std_error", "Enrichment",
       "Enrichment_std_error", "Enrichment_p", "Coefficient", "Coefficient_std_error",
       "Coefficient_z-score"]

with open(f"{AGG}/tau_star.json") as fh:
    tau = {k: v for k, v in json.load(fh).items() if not k.startswith("__")}
rows = []

# ---- prediction arm: the tau* row of the threshold sweep ----
jk = pd.read_csv(f"{AGG}/tau_star_jackknife.tsv", sep="\t", keep_default_na=False)
jk = jk[jk["model"] == SEL_MODEL]
for coh, t in tau.items():
    r = jk[(jk["cohort"] == coh) & (jk["pred_prob"] == t)]
    if r.empty:
        print(f"  {coh}: no {SEL_MODEL} row at tau*={t} -- skipped")
        continue
    r = r.iloc[0]
    rows.append({"Category": f"{coh}_pred_prob_binary_1",
                 "Prop._SNPs": r["prop_SNPs"], "Prop._h2": r["prop_h2"],
                 "Prop._h2_std_error": r["prop_h2_SE"],
                 "Enrichment": r["enrichment"], "Enrichment_std_error": r["enrichment_SE"],
                 "Enrichment_p": r["enrichment_p"],
                 "Coefficient": r["coefficient"], "Coefficient_std_error": r["coefficient_SE"],
                 "Coefficient_z-score": r["coef_z_ldsc"],
                 "cohort": coh, "arm": "prediction", "model": SEL_MODEL, "threshold": t})
    print(f"  {coh}: prediction arm at tau*={t}")

# ---- pip arm: the focal annotation is the last Category of the .results ----
for coh in tau:
    rf = f"{PIP_RES}/{coh}/{coh}_pip_binary.results"
    if not os.path.exists(rf):
        print(f"  {coh}: missing {rf}")
        continue
    res = pd.read_csv(rf, sep="\t")
    target = res["Category"].tolist()[-1]
    if "pip" not in target:
        print(f"  {coh}: last Category is {target!r}, not the PIP annotation -- skipped")
        continue
    r = res.iloc[-1]
    rows.append({"Category": f"{coh}_pip_binary_1",
                 **{c: r[c] for c in PUB[1:]},
                 "cohort": coh, "arm": "pip", "model": "NA", "threshold": 0.10})
    print(f"  {coh}: pip arm ({target!r})")

out = pd.DataFrame(rows)
if out.empty:
    raise SystemExit("nothing to write: neither arm produced rows")
out = out[PUB + ["cohort", "arm", "model", "threshold"]]
out.to_csv(OUT, sep="\t", index=False)
n_pred, n_pip = int((out["arm"] == "prediction").sum()), int((out["arm"] == "pip").sum())
print(f"\nwrote {OUT}: {len(out)} rows ({n_pred} prediction at tau*, {n_pip} PIP > 0.10)")
if n_pip < len(tau):
    print(f"WARNING: only {n_pip}/{len(tau)} PIP arms present")
