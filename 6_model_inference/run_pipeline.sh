#!/bin/bash
# Step 6: Model Inference
# Score all variants genome-wide for each gene.

set -e

echo "=== Step 6: Model Inference ==="

COHORT=${1:-Mic_mega_eQTL}
GENE_LIST="../data/training_data/${COHORT}/list_genes.csv"
NUM_GENES=$(( $(wc -l < "$GENE_LIST") - 1 ))

echo "Scoring ${NUM_GENES} genes for cohort ${COHORT}..."

for i in $(seq 0 $(( NUM_GENES - 1 ))); do
    echo "Gene ${i}/${NUM_GENES}..."
    python model_inference.py $COHORT $i
done

echo "=== Step 6 Complete ==="
