#!/bin/bash
# Step 4 on one machine, for one cell type: the train, train_restricted and test sets of every chromosome.
#   bash run_pipeline.sh Mic_mega_eQTL [NPR]
set -eo pipefail
cd "$(dirname "$0")"
COHORT=${1:-Mic_mega_eQTL}
NPR=${2:-10}

for chr in $(seq 1 22); do
    python create_training_datasets.py "$chr" train "$COHORT" "$NPR"
    python create_training_datasets.py "$chr" train_restricted "$COHORT" "$NPR"
    python create_training_datasets.py "$chr" test "$COHORT" "$NPR"
done
