#!/bin/bash
# Minimal wrapper for single held-out chromosome training demo.
set -euo pipefail

COHORT=${1:-Mic_mega_eQTL}
CHR=${2:-2}
GENE_LOF_FILE=${3:-"../sample_data/41588_2024_1820_MOESM4_ESM.xlsx"}
YAML_PATH=${4:-data_params.yaml}

python train_model.py "${COHORT}" "${CHR}" \
  --gene_lof_file "${GENE_LOF_FILE}" \
  --yaml_path "${YAML_PATH}" \
  --single_chromosome_demo
