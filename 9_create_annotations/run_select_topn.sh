#!/bin/bash
#SBATCH --job-name=sceems_select_topn
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1

# Genome-wide top-N selection for the size-matched comparison, one job per cell type (the ranking must be
# genome-wide, so this is not split by chromosome). Needs the step 7 predictions of all three models and
# the step 11 credible sets:
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_select_topn.sh            (N=5000 by default)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

python -u select_topn.py "${cohort}" ${N:-5000}
