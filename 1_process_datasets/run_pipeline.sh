#!/bin/bash
# Step 1 on one machine: list the genes of the fine-mapping exports, extract every gene's PIPs and top loci,
# convert them to parquet, collect each cell type's "other" genes and list the unique variants of each
# chromosome. Needs the FunGen-xQTL SuSiE exports (finemapping_rds_dir), which are not distributed.
# The Python scripts run in the scEEMS environment and get_vars_pips.R in scEEMS_R (environment_r.yml).
# With about 18,000 genes, run get_vars_pips.R as an array job, one task per gene, on a cluster.
set -eo pipefail
cd "$(dirname "$0")"
source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate scEEMS
python create_job_list.py
GENES=$(python -c "import sys; sys.path.insert(0, '../shared'); from config import path; print(path('susie_dir') + '/genes_list.txt')")

conda activate scEEMS_R
for i in $(seq 1 "$(wc -l < "$GENES")"); do
    Rscript --vanilla get_vars_pips.R "$i"
done

conda activate scEEMS
python create_parquet_files.py
for cell in Ast Exc Inh Mic Oli OPC; do
    python create_parquet_files_other.py "$cell"
done
for chr in $(seq 1 22); do
    python create_unique_variant_list.py "$chr"
done
