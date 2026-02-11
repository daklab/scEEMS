#!/bin/bash
# Step 2: Annotate Variants
# Run variant annotation for each chromosome.

set -e

echo "=== Step 2: Annotate Variants ==="

for chr in $(seq 1 22); do
    echo "Annotating chromosome ${chr}..."
    python annotate_variants.py $chr
done

echo "=== Step 2 Complete ==="
