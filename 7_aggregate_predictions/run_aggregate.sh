#!/bin/bash
#SBATCH --job-name=sceems_aggregate
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=10

# Collect the per-gene predictions of all three models of one cell type into parquet datasets. Run after
# step 6 has scored every gene (find_missing_genes.py reports none missing):
#   sbatch --export=ALL,cohort=Mic_mega_eQTL run_aggregate.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

for m in weighted_full unweighted_full weighted_restricted; do
    python -u create_parquet_scored.py "${cohort}" "${m}"
done
