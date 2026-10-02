#!/bin/bash
#SBATCH --job-name=sceems_annot
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=10
#SBATCH --array=1-22

# Threshold-sweep annotations (pred_prob > 0.80 ... 0.99) for one cell type and model, one task per
# chromosome. Run after step 7 for that model:
#   for m in weighted_full unweighted_full weighted_restricted; do
#       sbatch --export=ALL,cohort=Mic_mega_eQTL,model=$m run_annotations.sh; done
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u make_annotations.py "${SLURM_ARRAY_TASK_ID}" "${cohort}" "${model}"
