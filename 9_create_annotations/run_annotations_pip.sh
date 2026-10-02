#!/bin/bash
#SBATCH --job-name=sceems_annot_pip
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=10
#SBATCH --array=1-22

# Fine-mapped eQTL (PIP > 0.10) annotation for one cell type, one task per chromosome:
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_annotations_pip.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u make_annotations_pip.py "${SLURM_ARRAY_TASK_ID}" "${cohort}"
