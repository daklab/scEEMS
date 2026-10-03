#!/bin/bash
#SBATCH --job-name=sceems_sldsc
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=72G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1

# S-LDSC for the columns of one annotation set, one array task per column (baseline + that column).
# Run after the set's LD scores (all 22 chromosomes). Set the array to the number of columns:
#   sbatch --export=ALL,set=pareto/Mic_mega_eQTL_weighted_full --array=1-20 run_ldscore_regression.sh
#   sbatch --export=ALL,set=pip/Mic_mega_eQTL --array=1 run_ldscore_regression.sh
#   sbatch --export=ALL,set=top5000/Mic_mega_eQTL --array=1-9 run_ldscore_regression.sh
#   sbatch --export=ALL,set=cs/Mic_mega_eQTL --array=1-31 run_ldscore_regression.sh
# Runs in the PolyFun environment (POLYFUN_ENV, default polyfun).
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${POLYFUN_ENV:-polyfun}"
cd "${SLURM_SUBMIT_DIR}"

python -u ldscore_regression.py "${set}" "${SLURM_ARRAY_TASK_ID}"
