#!/bin/bash
# Step 5 on one machine, for one cell type: train the three models for each held-out chromosome, then
# compute the held-out AUPRC. Uses the feature weights of the data release; with SEARCH=1 it first runs
# the feature-weight search (odd and even chromosomes) and trains with that selection instead.
#   bash run_pipeline.sh Mic_mega_eQTL
#   SEARCH=1 bash run_pipeline.sh Mic_mega_eQTL
set -euo pipefail
COHORT=${1:-Mic_mega_eQTL}
WEIGHTS=""

if [ "${SEARCH:-0}" = 1 ]; then
    python -u search_feature_weights.py "$COHORT" odd
    python -u search_feature_weights.py "$COHORT" even
    python -u select_feature_weights.py "$COHORT"
    WEIGHTS=$(python -c "import sys; sys.path.insert(0, '../shared'); from config import path, cohort_name; c = cohort_name('$COHORT'); print(path('feature_weight_search_dir', cohort=c) + f'/best_configs_{c}.json')")
fi

for chr in $(seq 1 22); do
    python -u train_loco.py "$COHORT" "$chr" $WEIGHTS
done
python -u evaluate_auprc.py "$COHORT"
