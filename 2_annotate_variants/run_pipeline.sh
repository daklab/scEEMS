#!/bin/bash
# Step 2 on one machine: annotate the unique variants of each chromosome (step 1), then keep the identifier
# columns that steps 3 and 4 use to select annotated variants. Chromosome 21 needs about 60 GB of memory and
# the largest chromosomes several times more; on a cluster, run one job per chromosome.
set -eo pipefail
cd "$(dirname "$0")"

for chr in $(seq 1 22); do
    python annotate_variants.py "$chr"
    python subset_annotated_variants.py "$chr"
done
