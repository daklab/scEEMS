#!/bin/bash
#SBATCH --job-name=sceems_magma
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=5
#SBATCH --array=1-22

# MAGMA gene analysis, one task per chromosome (run_magma.py). Submit, after run_make_magma.sh:
#   for c in Ast Exc Inh Mic Oli OPC; do
#     for pop in EUR AFR AMR EAS; do
#       sbatch --export=ALL,annotation=emagma,cohort=${c}_mega_eQTL,pop=${pop} run_magma.sh; done
#     sbatch --export=ALL,annotation=knn,cohort=${c}_mega_eQTL run_magma.sh; done
#   sbatch --export=ALL,annotation=positional run_magma.sh       # once, after make_positional_annotation.py
: "${annotation:?set annotation=emagma, knn or positional}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scEEMS
cd "${SLURM_SUBMIT_DIR}"

if [ "${annotation}" = positional ]; then
    python -u run_magma.py "${SLURM_ARRAY_TASK_ID}" positional
else
    : "${cohort:?set cohort=Mic_mega_eQTL}"
    python -u run_magma.py "${SLURM_ARRAY_TASK_ID}" "${annotation}" "${cohort}" "${pop:-EUR}"
fi
