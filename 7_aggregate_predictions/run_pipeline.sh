#!/bin/bash
# Step 7: Aggregate Predictions
# Combine per-gene predictions into chromosome-level Parquet files.

set -e

echo "=== Step 7: Aggregate Predictions ==="

for cohort in Ast_mega_eQTL Exc_mega_eQTL Inh_mega_eQTL Mic_mega_eQTL Oli_mega_eQTL OPC_mega_eQTL; do
    echo "Aggregating predictions for ${cohort}..."
    python aggregate_predictions.py $cohort
done

echo "=== Step 7 Complete ==="
