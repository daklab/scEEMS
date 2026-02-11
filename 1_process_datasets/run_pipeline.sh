#!/bin/bash
# Step 1: Process Datasets
# Internal provenance scripts for creating training inputs from raw fine-mapping RDS files.
# Raw inputs are not public; this step is documented but not expected to run for public users.

set -e

echo "=== Step 1: Process Datasets ==="

# 1. Raw fine-mapping RDS download is not public
echo "Step 1 raw inputs are not publicly available."
echo "If you have private internal RDS inputs under paths.data_dir/release_04_2024, continuing..."

# 2. Create gene list from downloaded files
echo "Creating gene list..."
python create_job_list.py

# 3. Extract variant PIPs from RDS files (run per gene)
echo "Extracting variant PIPs..."
NUM_GENES=$(wc -l < genes_list.txt)
for i in $(seq 1 $NUM_GENES); do
    Rscript get_vars_pips.R $i
done

# 4. Convert CSV results to Parquet
echo "Converting to Parquet..."
python create_parquet_files.py

# 5. Create unique variant lists per chromosome
echo "Creating unique variant lists..."
for chr in $(seq 1 22); do
    python create_unique_variant_list.py $chr
done

echo "=== Step 1 Complete ==="
