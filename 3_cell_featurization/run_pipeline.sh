#!/bin/bash
# Step 3: Cell Featurization
# Create per-gene feature matrices for each cohort.
# This step is parallelized per gene.

set -e

echo "=== Step 3: Cell Featurization ==="

COHORT=${1:-Mic_mega_eQTL}
OTHER=${2:-F}

# Count number of genes
GENE_LIST="../data/training_data/${COHORT}/list_genes.csv"
NUM_GENES=$(tail -n +2 "$GENE_LIST" | wc -l)

echo "Processing ${NUM_GENES} genes for cohort ${COHORT}..."

for i in $(seq 1 $NUM_GENES); do
    echo "Gene ${i}/${NUM_GENES}..."
    python create_gene_datasets.py $i $COHORT $OTHER
done

echo "=== Step 3 Complete ==="
