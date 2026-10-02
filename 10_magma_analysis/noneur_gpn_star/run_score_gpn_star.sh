#!/bin/bash
#SBATCH --job-name=sceems_gpn_star
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user@example.com
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

# Optional, step 2 of 3: GPN-STAR scores of the non-European panel's SNVs, one task per chromosome, on a
# GPU, in an environment with the gpn package (https://github.com/songlab-cal/gpn). msa is the GPN-STAR
# multiple sequence alignment directory and model the GPN-STAR model (the published scores used the p243
# alignment, a window of 256 and the gpn-star-p243-200m model); dir is {magma_dir}/noneur_gpn_star.
#   sbatch --export=ALL,dir=/path/to/magma/noneur_gpn_star,msa=/path/to/msa/p243,model=/path/to/gpn-star-p243-200m run_score_gpn_star.sh
# Then: python merge_scores_noneur.py
: "${dir:?set dir}" "${msa:?set msa}" "${model:?set model}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gpn                       # an environment with the gpn package
cd "${SLURM_SUBMIT_DIR}"

CHR=${SLURM_ARRAY_TASK_ID}
mkdir -p "${dir}/snv_scored"
python -u score_gpn_star_chunked.py "${dir}/snv_input/chr${CHR}.parquet" "${dir}/snv_scored/chr${CHR}.parquet" \
    "${msa}" 256 "${model}" --chunk 100000 --batch 16 --workers 8 --workdir "${dir}/work/chr${CHR}"
