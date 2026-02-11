#!/bin/bash
# Step 8: SHAP Analysis
# Compute SHAP values for model interpretability.

set -e

echo "=== Step 8: SHAP Analysis ==="

for cohort in Ast_mega_eQTL Exc_mega_eQTL Inh_mega_eQTL Mic_mega_eQTL Oli_mega_eQTL OPC_mega_eQTL; do
    echo "Computing SHAP values for ${cohort}..."
    for chr in $(seq 1 22); do
        python shap_analysis.py $cohort $chr
    done
    echo "Aggregating SHAP scores for ${cohort}..."
    python aggregate_shap_scores.py $cohort
done

echo "Computing cross-cell-type summary..."
python summary_cell_types_shap.py

echo "=== Step 8 Complete ==="
