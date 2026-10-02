#!/bin/bash
# Step 7 on one machine: collect the per-gene predictions of every cell type into parquet datasets, then
# export the scEEMS (weighted_full) predictions as tabix-indexed TSVs.
set -euo pipefail

for cohort in Ast_mega_eQTL Exc_mega_eQTL Inh_mega_eQTL Mic_mega_eQTL Oli_mega_eQTL OPC_mega_eQTL; do
    for m in weighted_full unweighted_full weighted_restricted; do
        python create_parquet_scored.py "$cohort" "$m"
    done
    for chr in $(seq 1 22); do
        python extract_predictions_tsv.py "$chr" "$cohort"
    done
done
