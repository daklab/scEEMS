#!/bin/bash
# Step 8 on one machine: TreeSHAP for every chromosome of every cell type, then the summary tables.
# On a cluster, use run_shap.sh (one SLURM array task per chromosome) instead.
set -euo pipefail

for cohort in Ast_mega_eQTL Exc_mega_eQTL Inh_mega_eQTL Mic_mega_eQTL Oli_mega_eQTL OPC_mega_eQTL; do
    for chr in $(seq 1 22); do
        python shap_analysis.py "$cohort" "$chr"
    done
done
python aggregate_shap.py
python summary_to_wide.py
