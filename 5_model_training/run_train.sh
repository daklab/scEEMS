#!/bin/bash
#SBATCH --job-name=sceems_train
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=80G
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=10
#SBATCH --array=1-22

# LOCO training of the three models for one cell type, one array task per held-out chromosome.
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_train.sh
# To train with your own feature-weight selection instead of the data release's, add
#   weights=/path/to/best_configs_Mic_mega_eQTL.json   to --export.
# When all 22 tasks have finished:  python evaluate_auprc.py Mic_mega_eQTL
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u train_loco.py "${cohort}" "${SLURM_ARRAY_TASK_ID}" ${weights:-}
