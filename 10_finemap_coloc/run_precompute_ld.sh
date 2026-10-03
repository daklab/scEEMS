#!/bin/bash
#SBATCH --job-name=sceems_precompute_ld
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=50G
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=5

# LD of each gene's GWAS window, one array task per row of {finemap_dir}/genes.tsv (make_gene_list.R
# prints the count):
#   sbatch --array=1-13206%150 run_precompute_ld.sh
# Genes with very wide windows (mostly pericentromeric chr9) need --mem=100G --cpus-per-task=10; rerun
# the array with those settings and finished genes are skipped.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS_R
cd "${SLURM_SUBMIT_DIR}"

Rscript --vanilla precompute_ld.R "${SLURM_ARRAY_TASK_ID}"
