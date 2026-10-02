#!/usr/bin/env python
"""
Credible-set S-LDSC, last step: collect the runs of the credible-set annotation columns into one table.

Each run fitted the baseline plus ONE focal column, so a row is that annotation's contribution over the
baseline alone, and comparing rows is a ranking over separate models, not a variance decomposition.

  arm == "finemapped"   {cell}_cs: every 95% credible-set member of the FunGen-xQTL fine-mapping, eQTL +
                        other genes, no PIP cut (the comparator of the predicted eQTLs at tau*)
  arm == "prior"        {cell}_{prior}_cs_pip{NN}: credible-set members with PIP > NN/100 under each
                        fine-mapping prior, eQTL genes fine-mapped under all priors; pip00 is PIP > 0,
                        the whole credible set

usage:   python aggregate_cs_sldsc.py
output:  {aggregate_dir}/finemap_cs_sldsc.tsv: the focal row of every run, plus arm, cell_type, prior,
         pip_threshold, pip_label, gene_set, cohort, column
"""
import glob
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

RES = f"{path('sldsc_dir')}/results/cs"
OUT = f"{path('aggregate_dir')}/finemap_cs_sldsc.tsv"

PRIOR_RE = re.compile(r"^(?P<cell>[A-Za-z]+)_(?P<prior>.+)_cs_pip(?P<thr>\d{2})$")
FM_RE = re.compile(r"^(?P<cell>[A-Za-z]+)_cs$")

rows = []
for rf in sorted(glob.glob(f"{RES}/*/*.results")):
    column = os.path.basename(rf)[: -len(".results")]
    df = pd.read_csv(rf, sep="\t")
    focal = df[df["Category"].str.startswith(column)]
    if len(focal) != 1:
        print(f"  SKIP {rf}: {len(focal)} rows match {column}", flush=True)
        continue
    r = focal.iloc[0].to_dict()
    r["cohort"] = os.path.basename(os.path.dirname(rf))
    r["column"] = column
    m = FM_RE.match(column)
    if m:
        r.update(arm="finemapped", cell_type=m["cell"], prior="uniform_original",
                 pip_threshold=0.0, pip_label="PIP > 0 (whole credible set)", gene_set="MEGA+Other")
    else:
        m = PRIOR_RE.match(column)
        if not m:
            print(f"  SKIP {column}: unparsed column name", flush=True)
            continue
        thr = int(m["thr"]) / 100
        r.update(arm="prior", cell_type=m["cell"], prior=m["prior"], pip_threshold=thr,
                 pip_label=("PIP > 0 (whole credible set)" if thr == 0 else f"PIP > {thr:.2f}"),
                 gene_set="MEGA_in_all_priors")
    rows.append(r)

if not rows:
    raise SystemExit(f"no .results files under {RES}")
out = pd.DataFrame(rows)
lead = ["arm", "cell_type", "prior", "pip_threshold", "pip_label", "gene_set", "cohort", "column"]
out = out[lead + [c for c in out.columns if c not in lead]]
out = out.sort_values(["arm", "cell_type", "prior", "pip_threshold"])
out.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}  ({len(out)} rows, {out.cohort.nunique()} cell types)")
print(out[["arm", "cell_type", "prior", "pip_label", "Prop._SNPs", "Prop._h2",
           "Enrichment", "Enrichment_p"]].to_string(index=False))
