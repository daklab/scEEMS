#!/bin/bash
#SBATCH --job-name=sceems_annot_topn
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

# Size-matched top-N annotations for one cell type, one task per chromosome. Run after run_select_topn.sh:
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_annotations_topn.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u make_annotations_topn.py "${SLURM_ARRAY_TASK_ID}" "${cohort}" ${N:-5000}
