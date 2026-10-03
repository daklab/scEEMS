#!/bin/bash
# Step 11 on one machine: gene lists, cross-cell-type colocalization of every gene, aggregation.
# On a cluster, use run_coloc_crosscell.sh (one SLURM array task per gene) for the middle part.
set -euo pipefail

N=$(python make_gene_lists.py | awk 'NR == 1 {gsub(",", "", $1); print $1}')
for i in $(seq 1 "$N"); do
    Rscript --vanilla coloc_crosscell.R "$i"
done
python aggregate_crosscell.py
