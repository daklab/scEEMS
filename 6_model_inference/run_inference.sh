#!/bin/bash
#SBATCH --job-name=sceems_inference
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=30G
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=2
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Score every gene of one cell type, one array task per gene (row of the gene list). Run after step 5
# (all 22 held-out chromosomes) and after make_gene_lists.py, which prints the array size, e.g.
#   sbatch --export=ALL,cohort=Mic_mega_eQTL --array=1-7100%200 run_inference.sh
# Logs go to /dev/null because there is one task per gene; find_missing_genes.py lists genes to rerun.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u model_inference.py "${cohort}" "${SLURM_ARRAY_TASK_ID}"
