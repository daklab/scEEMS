#!/bin/bash
# Step 3 on one machine, for one cell type: list its genes, then build the feature table of every gene, those
# fine-mapped by the MEGA analysis (F) and the "other" genes (T). There are 7,000-11,500 genes per cell type;
# on a cluster, run create_gene_datasets.py as an array job, one task per gene.
#   bash run_pipeline.sh Mic_mega_eQTL
set -eo pipefail
cd "$(dirname "$0")"
COHORT=${1:-Mic_mega_eQTL}

python create_gene_lists.py "$COHORT"
for mode in F T; do
    N=$(python -c "
import os, sys
sys.path.insert(0, '../shared')
from config import cohort_name, path
f = os.path.join(path('gene_list_dir', cohort=cohort_name('$COHORT')), 'list_genes.csv' if '$mode' == 'F' else 'list_genes_other.csv')
print(sum(1 for _ in open(f)) - 1 if os.path.exists(f) else 0)")
    for i in $(seq 1 "$N"); do
        python create_gene_datasets.py "$i" "$COHORT" "$mode"
    done
done
