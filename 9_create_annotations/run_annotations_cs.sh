#!/bin/bash
#SBATCH --job-name=sceems_annot_cs
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

# Credible-set annotations (31 columns) for one cell type, one task per chromosome. Needs the step 10
# credible sets and the step 1 PIP_top exports:
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_annotations_cs.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u make_annotations_cs.py "${SLURM_ARRAY_TASK_ID}" "${cohort}"
