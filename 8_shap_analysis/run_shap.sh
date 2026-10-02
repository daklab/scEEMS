#!/bin/bash
#SBATCH --job-name=sceems_shap
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=5
#SBATCH --array=1-22

# TreeSHAP of the predicted eQTLs of one cell type, one array task per chromosome. Needs the step 5
# models, the step 7 predictions, the step 3 all_variants tables and tau_star.json from step 9:
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_shap.sh
# When all cell types have finished:  python aggregate_shap.py && python summary_to_wide.py
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u shap_analysis.py "${cohort}" "${SLURM_ARRAY_TASK_ID}"
