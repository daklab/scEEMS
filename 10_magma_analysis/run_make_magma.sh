#!/bin/bash
#SBATCH --job-name=sceems_magma_annot
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=5
#SBATCH --array=1-22

# eMAGMA gene annotations (prediction and pip) and their TSS-distance controls for one cell type, one task
# per chromosome. Needs tau* ({aggregate_dir}/tau_star.json, step 9) and the non-European predictions
# (run_score_noneur.sh).
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_make_magma.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u make_magma_files.py "${SLURM_ARRAY_TASK_ID}" "${cohort}"
python -u make_magma_files_knn.py "${SLURM_ARRAY_TASK_ID}" "${cohort}"
