#!/bin/bash
# Step 5: Model Training
# Train CatBoost classifiers using LOCO cross-validation.

set -e

echo "=== Step 5: Model Training ==="

COHORT=${1:-Mic_mega_eQTL}
GENE_LOF_FILE=${2:-"../data/41588_2024_1820_MOESM4_ESM.xlsx"}

for chr in $(seq 1 22); do
    echo "Training model with chr${chr} held out..."
    python train_model.py $COHORT $chr \
        --gene_lof_file "$GENE_LOF_FILE" \
        --yaml_path data_params.yaml
done

echo "=== Step 5 Complete ==="
