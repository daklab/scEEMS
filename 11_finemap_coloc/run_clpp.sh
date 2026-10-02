#!/bin/bash
#SBATCH --job-name=sceems_clpp
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=24G
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

# CLPP for every cell type, eQTL prior and GWAS prior, one task per chromosome. Run after all
# fine-mapping is done; then python aggregate_clpp.py.
#   sbatch run_clpp.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS_R
cd "${SLURM_SUBMIT_DIR}"

Rscript --vanilla compute_clpp.R "chr${SLURM_ARRAY_TASK_ID}"
