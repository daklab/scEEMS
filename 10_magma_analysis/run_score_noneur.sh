#!/bin/bash
#SBATCH --job-name=sceems_magma_noneur
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=50G
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=5
#SBATCH --array=1-22

# Score the non-European MAGMA panel with the scEEMS model of one cell type, one task per chromosome. Needs
# the step 5 models and the non-European GPN-STAR scores (noneur_gpn_star/).
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_score_noneur.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u score_noneur.py "${cohort}" "${SLURM_ARRAY_TASK_ID}"
