#!/bin/bash
# Step 9 on one machine, for one cell type: the four annotation sets, their LD scores and S-LDSC runs, then
# the result tables. Needs step 7 (all three models) and step 10. The LD scores and S-LDSC runs use the
# PolyFun environment (POLYFUN_ENV, default polyfun) through `conda run`; the rest runs in the current
# environment. tau_star.py and the aggregation scripts summarize every cell type with results, so rerun
# them after the last cell type.
#   bash run_pipeline.sh Mic_mega_eQTL
set -euo pipefail
COHORT=${1:-Mic_mega_eQTL}
N=${N:-5000}
PF="conda run --no-capture-output -n ${POLYFUN_ENV:-polyfun} python -u"

# ---- annotations
sets=()
for m in weighted_full unweighted_full weighted_restricted; do
    for chr in $(seq 1 22); do python make_annotations.py "$chr" "$COHORT" "$m"; done
    sets+=("pareto/${COHORT}_${m}")
done
python select_topn.py "$COHORT" "$N"
for chr in $(seq 1 22); do
    python make_annotations_pip.py "$chr" "$COHORT"
    python make_annotations_topn.py "$chr" "$COHORT" "$N"
    python make_annotations_cs.py "$chr" "$COHORT"
done
sets+=("pip/${COHORT}" "top${N}/${COHORT}" "cs/${COHORT}")

# ---- LD scores, then one S-LDSC run per column
for s in "${sets[@]}"; do
    for chr in $(seq 1 22); do $PF compute_ldscores.py "$s" "$chr"; done
    ncol=$(python -c "import sys, pandas as pd; sys.path.insert(0, '../shared'); from config import path; \
print(len(pd.read_csv(path('sldsc_dir') + '/$s/MLxQTL_chr22.annot.gz', sep='\t', nrows=0).columns) - 5)")
    for i in $(seq 1 "$ncol"); do $PF ldscore_regression.py "$s" "$i"; done
done

# ---- tables
python tau_star.py
python aggregate_prediction_vs_pip.py
python aggregate_topn_sldsc.py "$N"
python aggregate_cs_sldsc.py
