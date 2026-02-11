#!/bin/bash
# Step 9: Create LDSC annotations and run partitioned heritability analysis
# This script shows the full execution sequence. Adapt for your cluster.

set -e

CELL_TYPES=("Mic_mega_eQTL" "Ast_mega_eQTL" "Exc_mega_eQTL" "Inh_mega_eQTL" "Oli_mega_eQTL" "OPC_mega_eQTL")

echo "=== Step 9a: Create annotations ==="
for cell_type in "${CELL_TYPES[@]}"; do
    echo "Submitting annotation jobs for ${cell_type}..."
    # Revised model
    sbatch --export=cell_type=${cell_type} run_jobs.sh
    # Original weighted model
    sbatch --export=cell_type=${cell_type} run_jobs_original.sh
    # Unweighted model
    sbatch --export=cell_type=${cell_type} run_jobs_unweighted.sh
done

echo "Wait for all annotation jobs to complete before proceeding."
echo ""

echo "=== Step 9b: Compute LD scores ==="
for cell_type in "${CELL_TYPES[@]}"; do
    echo "Submitting LD score jobs for ${cell_type}..."
    # Set ANNOT_DIR for each variant before submitting
    # Revised
    sbatch --export=cell_type=${cell_type} template_ldscore.sh
    # Original
    sbatch --export=cell_type=${cell_type} template_ldscore_original.sh
    # Unweighted
    sbatch --export=cell_type=${cell_type} template_ldscore_unweighted.sh
done

echo "Wait for all LD score jobs to complete before proceeding."
echo ""

echo "=== Step 9c: Run LDSC regression ==="
cd ldscore_regression
for cell_type in "${CELL_TYPES[@]}"; do
    # Revised model
    sbatch run_jobs.sh
    # Original model
    sbatch --export=cell_type=${cell_type} run_jobs_original.sh
    # Unweighted model
    sbatch run_jobs_unweighted.sh
done

echo "Wait for all regression jobs to complete before proceeding."
echo ""

echo "=== Step 9d: Aggregate results ==="
python aggregate_data.py

echo "Step 9 complete."
