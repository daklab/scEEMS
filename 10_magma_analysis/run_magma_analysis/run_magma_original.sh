#!/bin/bash
#SBATCH --job-name=magma           # Job name
#SBATCH --partition=cpu                      # Partition Name
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.com      # Where to send mail
#SBATCH --mem=100G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=4:20:00                       # Time limit 4 hours
#SBATCH --cpus-per-task=10               # Standard output and error log
#SBATCH --array=[2-2]%1

# ---------------------------------------------------------------
# Required shell variables (set these from config.yaml before running):
#   MAGMA_DIR       - paths.magma_dir       (directory containing MAGMA binary)
#   PLINK_REF_DIR   - paths.plink_ref_dir   (PLINK reference panel directory)
#   MAGMA_OUTPUT_DIR - paths.output_dir/aggregate_results/MAGMA
# ---------------------------------------------------------------

pop=EUR

# Activate the environment
conda activate scEEMS

chr=${SLURM_ARRAY_TASK_ID}

plink_file=${PLINK_REF_DIR}/${pop}/plink/ADSP_${pop}_chr${chr}

sumstats=${MAGMA_DIR}/bellenguez_MAGMA_sumstats.txt

annot_file_original=${MAGMA_OUTPUT_DIR}/annotations_gene/bellenguez_MAGMA.genes.annot.genes.annot

write_dir=${MAGMA_OUTPUT_DIR}/MAGMA_output_original

mkdir -p ${write_dir}

write_file=${write_dir}/original_chr${chr}_MAGMA_prediction

$MAGMA_DIR/magma --bfile ${plink_file} \
--gene-annot ${annot_file_original} \
--pval ${sumstats} ncol=N \
--gene-settings adap-permp=10000 \
--out ${write_file}
