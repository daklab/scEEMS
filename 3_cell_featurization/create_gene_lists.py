#!/usr/bin/env python
"""
List the genes step 3 featurizes for one cell type:
  list_genes.csv        the genes fine-mapped by the cell type's MEGA eQTL analysis
  list_genes_other.csv  its "other" genes, fine-mapped only by its DeJager or Kellis analysis
                        (1_process_datasets/create_parquet_files_other.py)
Both are tab-separated with columns gene_id, chr. A gene's row number is its index for
create_gene_datasets.py.

Usage:
    python create_gene_lists.py <cohort>        (e.g. Mic_mega_eQTL or Mic)

Input:   {susie_pips_dir}/PIP_all_parquet/PIP_all.parquet and PIP_all_other_parquet/PIP_all_other.parquet
Output:  {gene_list_dir}/list_genes.csv and list_genes_other.csv
"""
import os
import sys

from dask import dataframe as dd
from dask.diagnostics import ProgressBar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import cohort_name, path

ProgressBar().register()

cohort = cohort_name(sys.argv[1])
pips_dir = path("susie_pips_dir", cohort=cohort)
out_dir = path("gene_list_dir", cohort=cohort)
os.makedirs(out_dir, exist_ok=True)

for out_name, parquet_file in (("list_genes.csv", f'{pips_dir}/PIP_all_parquet/PIP_all.parquet'),
                               ("list_genes_other.csv", f'{pips_dir}/PIP_all_other_parquet/PIP_all_other.parquet')):
    if not os.path.exists(parquet_file):
        print(f"{cohort}: no {parquet_file}; {out_name} not written")
        continue
    genes_df = dd.read_parquet(parquet_file)
    list_genes_df = genes_df[['gene_id', 'chr']].drop_duplicates().compute()
    print(list_genes_df.shape, cohort, out_name)
    list_genes_df.to_csv(f'{out_dir}/{out_name}', index=False, sep='\t')
