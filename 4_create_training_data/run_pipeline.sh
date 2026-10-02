#!/bin/bash
# Step 4: Create Training Data
# Sample positive/negative training examples per chromosome.

set -e

echo "=== Step 4: Create Training Data ==="

COHORT=${1:-Mic_mega_eQTL}
NPR=${2:-10}

for chr in $(seq 1 22); do
    echo "Creating training data for chr${chr}..."
    python create_training_datasets.py $chr train $COHORT $NPR
    python create_training_datasets.py $chr test $COHORT $NPR
    python create_training_datasets.py $chr train_restricted $COHORT $NPR
done

echo "=== Step 4 Complete ==="
