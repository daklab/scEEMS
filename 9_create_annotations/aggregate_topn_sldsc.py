#!/usr/bin/env python
"""
Size-matched S-LDSC, last step: collect the runs of the top-N annotations (one per ranking) into one table.

Every column holds exactly N baseline variants, so Prop._SNPs is the same for all of them and enrichment
differences cannot come from annotation size. Each run fitted the baseline plus ONE focal column, so
comparing rows is a ranking over separate models, not a variance decomposition.

usage:   python aggregate_topn_sldsc.py [N=5000]
output:  {aggregate_dir}/topn_sldsc.tsv: the focal row of every run, plus
         arm (Predicted | Fine-mapped | Fine-mapped (re-fit prior)), label, cell_type, ranking,
         gene_set (MEGA+Other | MEGA_in_all_priors), n_target, cohort, column
"""
import glob
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

N = sys.argv[1] if len(sys.argv) > 1 else "5000"
RES = f"{path('sldsc_dir')}/results/top{N}"
OUT = f"{path('aggregate_dir')}/topn_sldsc.tsv"

COL_RE = re.compile(rf"^(?P<cell>[A-Za-z]+)_top{N}_(?P<ranking>.+)$")
MODEL = {"pred_weighted_full": "Weighted (Full)", "pred_unweighted_full": "Unweighted (Full)",
         "pred_weighted_restricted": "Weighted (Restricted)"}
PRIOR = {"pip_uniform": "Uniform", "pip_EMS": "EMS", "pip_scEEMS_Weighted_Full": "Weighted (Full)",
         "pip_scEEMS_Weighted_Restricted": "Weighted (Restricted)",
         "pip_scEEMS_Unweighted_Full": "Unweighted (Full)"}

rows = []
for rf in sorted(glob.glob(f"{RES}/*/*.results")):
    column = os.path.basename(rf)[: -len(".results")]
    m = COL_RE.match(column)
    if not m:
        print(f"  SKIP {column}: unparsed", flush=True)
        continue
    df = pd.read_csv(rf, sep="\t")
    focal = df[df["Category"].str.startswith(column)]
    if len(focal) != 1:
        print(f"  SKIP {rf}: {len(focal)} rows match", flush=True)
        continue
    r = focal.iloc[0].to_dict()
    rk = m["ranking"]
    r.update(cohort=os.path.basename(os.path.dirname(rf)), column=column,
             cell_type=m["cell"], ranking=rk, n_target=int(N))
    if rk in MODEL:
        r.update(arm="Predicted", label=MODEL[rk], gene_set="MEGA+Other")
    elif rk == "pip":
        # the FunGen-xQTL fine-mapping (uniform prior) of all genes: the comparator of the predictions
        r.update(arm="Fine-mapped", label="FunGen fine-mapping", gene_set="MEGA+Other")
    else:
        r.update(arm="Fine-mapped (re-fit prior)", label=PRIOR.get(rk, rk), gene_set="MEGA_in_all_priors")
    rows.append(r)

if not rows:
    raise SystemExit(f"no .results under {RES}")
out = pd.DataFrame(rows)
lead = ["arm", "label", "cell_type", "ranking", "gene_set", "n_target", "cohort", "column"]
out = out[lead + [c for c in out.columns if c not in lead]].sort_values(["arm", "cell_type", "label"])
out.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}  ({len(out)} rows)")
# the size match: Prop._SNPs must be constant within a cell type
chk = out.groupby("cell_type")["Prop._SNPs"].agg(["min", "max", "nunique"])
print("\nProp._SNPs per cell type (min and max are equal if the size match held):")
print(chk.to_string())
print()
print(out[["arm", "label", "cell_type", "Prop._SNPs", "Prop._h2", "Enrichment", "Enrichment_p"]].to_string(index=False))
