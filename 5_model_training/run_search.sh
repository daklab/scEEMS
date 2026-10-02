#!/bin/bash
#SBATCH --job-name=sceems_search
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=100G
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --array=1-2
#SBATCH --requeue

# Feature-category weight search for one cell type: array task 1 searches the odd chromosomes, task 2 the
# even ones (83 trials each). The study is stored on disk, so a requeued or resubmitted task resumes.
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_search.sh
# When both tasks have finished:  python select_feature_weights.py Mic_mega_eQTL
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

PARITY=$([ "${SLURM_ARRAY_TASK_ID}" -eq 1 ] && echo odd || echo even)
python -u search_feature_weights.py "${cohort}" "${PARITY}"
