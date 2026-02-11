#!/bin/bash
# Step 10: MAGMA gene-set enrichment analysis
# This script shows the full execution sequence. Adapt for your cluster.

set -e

CELL_TYPES=("Mic_mega_eQTL" "Ast_mega_eQTL" "Exc_mega_eQTL" "Inh_mega_eQTL" "Oli_mega_eQTL" "OPC_mega_eQTL")

echo "=== Step 10a: Prepare summary statistics ==="
cd run_magma_analysis
python make_MAGMA_sumstats.py
python make_MAGMA_sumstats_v2.py
cd ..

echo "=== Step 10b: Create MAGMA annotation files ==="
for cell_type in "${CELL_TYPES[@]}"; do
    echo "Submitting MAGMA file generation for ${cell_type}..."
    # Revised model
    sbatch --export=cell_type=${cell_type} run_jobs.sh
    # KNN model
    sbatch --export=cell_type=${cell_type} run_jobs_knn.sh
done

echo "Wait for all annotation jobs to complete before proceeding."
echo ""

echo "=== Step 10c: Run MAGMA gene-based tests ==="
cd run_magma_analysis
for cell_type in "${CELL_TYPES[@]}"; do
    echo "Submitting MAGMA analysis for ${cell_type}..."
    # Bellenguez (EUR)
    sbatch --export=cell_type=${cell_type} run_magma.sh
    # Kunkle (EUR)
    sbatch --export=cell_type=${cell_type} run_magma_kunkle.sh
    # Kunkle (AFR)
    sbatch --export=cell_type=${cell_type} run_magma_kunkle_AFR.sh
    # ADSP (AFR)
    sbatch --export=cell_type=${cell_type} run_magma_ADSP_AFR.sh
    # ADSP (AMR)
    sbatch --export=cell_type=${cell_type} run_magma_ADSP_AMR.sh
    # ADGC (AFR, AMR, EAS)
    sbatch --export=cell_type=${cell_type} run_magma_ADGC_AFR.sh
    sbatch --export=cell_type=${cell_type} run_magma_ADGC_AMR.sh
    sbatch --export=cell_type=${cell_type} run_magma_ADGC_EAS.sh
    # KNN
    sbatch --export=cell_type=${cell_type} run_magma_knn.sh
done
cd ..

echo "Wait for all MAGMA jobs to complete before proceeding."
echo ""

echo "=== Step 10d: Extract significant genes ==="
python extract_significant_MAGMA_genes.py
python extract_significant_MAGMA_genes_knn.py

echo "Step 10 complete."
