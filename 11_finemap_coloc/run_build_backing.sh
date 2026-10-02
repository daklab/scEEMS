#!/bin/bash
#SBATCH --job-name=sceems_backing
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=30G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22

# LD reference panel (ADSP European samples) as bigsnpr backing files, one task per chromosome.
#   sbatch run_build_backing.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS_R
cd "${SLURM_SUBMIT_DIR}"

Rscript --vanilla build_backing.R "${SLURM_ARRAY_TASK_ID}"
