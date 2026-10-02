#!/bin/bash
# Step 10 on one machine, in order. On a cluster, submit the run_*.sh scripts instead (README.md).
# Needs the step 5 models (non-European scoring), the weighted_full predictions of step 7 and tau* of step 9.
set -euo pipefail
COHORTS="Ast_mega_eQTL Exc_mega_eQTL Inh_mega_eQTL Mic_mega_eQTL Oli_mega_eQTL OPC_mega_eQTL"

python make_magma_sumstats.py
python make_positional_annotation.py
for chr in $(seq 1 22); do python run_magma.py "$chr" positional; done

for cohort in $COHORTS; do
    for chr in $(seq 1 22); do
        python score_noneur.py "$cohort" "$chr"
        python make_magma_files.py "$chr" "$cohort"
        python make_magma_files_knn.py "$chr" "$cohort"
        for pop in EUR AFR AMR EAS; do python run_magma.py "$chr" emagma "$cohort" "$pop"; done
        python run_magma.py "$chr" knn "$cohort"
    done
done

python aggregate_magma.py
