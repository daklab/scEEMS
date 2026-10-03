#!/bin/bash
#SBATCH --job-name=sceems_snpvar
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10

# PolyFun per-SNP priors (SNPVAR) for the GWAS, then one tabix-indexed table of them.
# Runs in PolyFun's own environment (set POLYFUN_ENV; default "polyfun"):
#   sbatch run_snpvar.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${POLYFUN_ENV:-polyfun}"
cd "${SLURM_SUBMIT_DIR}"

python -u compute_snpvar.py
python -u snpvar_to_tabix.py
