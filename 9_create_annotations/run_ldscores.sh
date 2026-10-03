#!/bin/bash
#SBATCH --job-name=sceems_ldscores
#SBATCH --partition=cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

# LD scores of one annotation set, one task per chromosome. Run after the set's annotations:
#   sbatch --export=ALL,set=pareto/Mic_mega_eQTL_weighted_full run_ldscores.sh
#   sbatch --export=ALL,set=pip/Mic_mega_eQTL run_ldscores.sh
#   sbatch --export=ALL,set=top5000/Mic_mega_eQTL run_ldscores.sh
#   sbatch --export=ALL,set=cs/Mic_mega_eQTL run_ldscores.sh
# The MHC region makes chromosome 6 the slowest task. Runs in the PolyFun environment (POLYFUN_ENV, default
# polyfun).
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${POLYFUN_ENV:-polyfun}"
cd "${SLURM_SUBMIT_DIR}"

# compute_ldscores.py is multithreaded through BLAS: match the thread count to the allocation, or it
# starts one thread per core of the node
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${OMP_NUM_THREADS}
export OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS}

python -u compute_ldscores.py "${set}" "${SLURM_ARRAY_TASK_ID}"
