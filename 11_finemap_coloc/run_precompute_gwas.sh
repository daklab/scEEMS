#!/bin/bash
#SBATCH --job-name=sceems_precompute_gwas
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1

# GWAS fine-mapping (uniform and PolyFun priors) of each gene's window, one array task per row of
# {finemap_dir}/genes.tsv. Run after run_precompute_ld.sh has finished:
#   sbatch --array=1-13206%200 run_precompute_gwas.sh
# A few percent of genes need more memory; rerun the array with --mem=150G (and, for genes left without
# a fit, --export=ALL,RESCUE_MAX_VAR=100000 --time=24:00:00); finished genes are skipped.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS_R
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1

Rscript --vanilla precompute_gwas.R "${SLURM_ARRAY_TASK_ID}"
