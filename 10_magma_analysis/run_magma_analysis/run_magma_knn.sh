#!/bin/bash
#SBATCH --job-name=magma           # Job name
#SBATCH --partition=cpu                      # Partition Name
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.com      # Where to send mail
#SBATCH --mem=50G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=0:20:00                       # Time limit 4 hours
#SBATCH --cpus-per-task=5               # Standard output and error log
#SBATCH --array=[1-22]%22

# ---------------------------------------------------------------
# Required shell variables (set these from config.yaml before running):
#   MAGMA_DIR       - paths.magma_dir       (directory containing MAGMA binary)
#   PLINK_REF_DIR   - paths.plink_ref_dir   (PLINK reference panel directory)
#   MAGMA_OUTPUT_DIR - paths.output_dir/aggregate_results/MAGMA_knn
# ---------------------------------------------------------------

pop=EUR

# Activate the environment
conda activate scEEMS

chr=${SLURM_ARRAY_TASK_ID}

plink_file=${PLINK_REF_DIR}/${pop}/plink/ADSP_${pop}_chr${chr}

sumstats=${MAGMA_DIR}/bellenguez_MAGMA_sumstats.txt

annot_file_pred=${MAGMA_OUTPUT_DIR}/${cell_type}/prediction/${cell_type}_chr${chr}_MAGMA.genes.annot

annot_file_pip=${MAGMA_OUTPUT_DIR}/${cell_type}/pip/${cell_type}_chr${chr}_MAGMA.genes.annot

write_dir=${MAGMA_OUTPUT_DIR}/${cell_type}/MAGMA_output

mkdir -p ${write_dir}

write_file_prediction=${write_dir}/${cell_type}_chr${chr}_MAGMA_prediction

write_file_pip=${write_dir}/${cell_type}_chr${chr}_MAGMA_pip

$MAGMA_DIR/magma --bfile ${plink_file} \
--gene-annot ${annot_file_pred} \
--pval ${sumstats} ncol=N \
--gene-settings adap-permp=10000 \
--out ${write_file_prediction}

$MAGMA_DIR/magma --bfile ${plink_file} \
--gene-annot ${annot_file_pip} \
--pval ${sumstats} ncol=N \
--gene-settings adap-permp=10000 \
--out ${write_file_pip}
