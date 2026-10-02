#!/bin/bash
#SBATCH --job-name=sceems_export
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

# Export the scEEMS predictions of one cell type as tabix-indexed TSVs, one task per chromosome:
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_export.sh
# Finished chromosomes are skipped, so resubmitting reruns only failed tasks.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u extract_predictions_tsv.py "${SLURM_ARRAY_TASK_ID}" "${cohort}"
