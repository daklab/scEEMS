#!/bin/bash
# Step 6 on one machine, for one cell type: score every gene with the three models of step 5.
#   bash run_pipeline.sh Mic_mega_eQTL
# On a cluster, use run_inference.sh (one SLURM array task per gene) instead.
set -euo pipefail
COHORT=${1:-Mic_mega_eQTL}

python split_gpn_by_chr.py
N=$(python make_gene_lists.py "$COHORT" | awk 'NR > 1 {print $2}')
for i in $(seq 1 "$N"); do
    python model_inference.py "$COHORT" "$i"
done
python find_missing_genes.py "$COHORT"
