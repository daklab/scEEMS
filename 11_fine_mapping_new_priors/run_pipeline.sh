#!/bin/bash
# Step 11: Fine-Mapping with New Priors
# Re-run SuSiE fine-mapping with scEEMS-informed priors.

set -e

echo "=== Step 11: Fine-Mapping with New Priors ==="

# 1. Download reference panel
echo "Downloading reference panel..."
python download_data.py

# 2. Run fine-mapping with predicted priors
CELL_TYPE=${1:-Mic}
REGION_FILE="path/to/region_list.txt"
NUM_REGIONS=$(wc -l < "$REGION_FILE")

echo "Running fine-mapping with predicted priors for ${CELL_TYPE}..."
for i in $(seq 1 $NUM_REGIONS); do
    Rscript template_predicted.R $i $CELL_TYPE
done

echo "Running fine-mapping with uniform priors for ${CELL_TYPE}..."
for i in $(seq 1 $NUM_REGIONS); do
    Rscript template_uniform.R $i $CELL_TYPE
done

# 3. Compare results
echo "Comparing PIPs..."
python compare_pips.py

echo "=== Step 11 Complete ==="
