#!/usr/bin/env python
"""
Collect the fine-mapping results of a cell type's "other" genes into its MEGA analysis folder: genes that
the cell type's MEGA eQTL analysis did not fine-map but its DeJager or Kellis eQTL analysis did. Step 3
featurizes these genes too (create_gene_datasets.py ... T), so every gene with fine-mapping results in the
cell type is scored, and step 9 includes their credible sets.

  PIP_all_other  PIP_all rows of every gene fine-mapped by DeJager or Kellis and not by MEGA (a gene
                 fine-mapped by both studies contributes the rows of both)
  PIP_top_other  top loci of the genes fine-mapped by exactly one of DeJager and Kellis, and not by MEGA

Usage:
    python create_parquet_files_other.py <cell_type>        (Ast, Exc, Inh, Mic, Oli or OPC)

Input:   {susie_dir}/<cell_type>_{mega,DeJager,Kellis}_eQTL/PIP_all/*_pips.csv, and PIP_top/*_top_loci.csv of
         DeJager and Kellis (get_vars_pips.R)
Output:  {susie_dir}/<cell_type>_mega_eQTL/PIP_all_other_parquet/PIP_all_other.parquet and
         PIP_top_other_parquet/PIP_top_other.parquet, partitioned by chromosome
"""
import os
import sys

from dask import dataframe as dd
from dask.diagnostics import ProgressBar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from config import path

ProgressBar().register()

cell_type = sys.argv[1]
mega_cohort = f'{cell_type}_mega_eQTL'
other_cohorts = [f'{cell_type}_DeJager_eQTL', f'{cell_type}_Kellis_eQTL']
mega_dir = path("susie_pips_dir", cohort=mega_cohort)
other_dirs = [path("susie_pips_dir", cohort=c) for c in other_cohorts]


def get_gene_ids(directory):
    return set([f.split('_')[0] for f in os.listdir(directory) if f.endswith('_pips.csv')])


mega_gene_ids = get_gene_ids(f'{mega_dir}/PIP_all')
cohort_gene_ids = [get_gene_ids(f'{d}/PIP_all') for d in other_dirs]

# ---- PIP_all_other: genes fine-mapped by either study, not by MEGA ----
gene_ids_to_process = (cohort_gene_ids[0] | cohort_gene_ids[1]) - mega_gene_ids
eligible_files = [f'{d}/PIP_all/{gene_id}_pips.csv' for d in other_dirs for gene_id in sorted(gene_ids_to_process)
                  if os.path.exists(f'{d}/PIP_all/{gene_id}_pips.csv')]
print(f"{mega_cohort}: {len(gene_ids_to_process)} other genes, {len(eligible_files)} PIP_all files")

combined_df = dd.read_csv(eligible_files, sep=' ', dtype={'pip': 'float64'})
write_path = f'{mega_dir}/PIP_all_other_parquet'
os.makedirs(write_path, exist_ok=True)
combined_df = combined_df.repartition(partition_size="100MB")
combined_df.to_parquet(f'{write_path}/PIP_all_other.parquet', engine='pyarrow', compression='snappy',
                       partition_on=['chr'], overwrite=True)

# ---- PIP_top_other: genes fine-mapped by exactly one of the two studies, not by MEGA ----
dejager_only = cohort_gene_ids[0] - cohort_gene_ids[1]
kellis_only = cohort_gene_ids[1] - cohort_gene_ids[0]
gene_ids_to_process = (dejager_only | kellis_only) - mega_gene_ids
eligible_files = [f'{d}/PIP_top/{gene_id}_top_loci.csv' for d in other_dirs for gene_id in sorted(gene_ids_to_process)
                  if os.path.exists(f'{d}/PIP_top/{gene_id}_top_loci.csv')]
print(f"{mega_cohort}: {len(gene_ids_to_process)} genes fine-mapped by one study only, {len(eligible_files)} PIP_top files")

combined_df = dd.read_csv(eligible_files, sep=' ', dtype={'pip': 'float64'})
write_path = f'{mega_dir}/PIP_top_other_parquet'
os.makedirs(write_path, exist_ok=True)
dtypes_mapping = {
    'chr': 'string',
    'pos': 'int64',
    'ref': 'string',
    'alt': 'string',
    'betahat': 'double',
    'sebetahat': 'double',
    'maf': 'double',
    'pip': 'float64',
    'cs_coverage_0.95': 'int64',
    'cs_coverage_0.7': 'int64',
    'cs_coverage_0.5': 'int64',
    'gene_id': 'string'
}
combined_df = combined_df.repartition(partition_size="100MB")
combined_df = combined_df.astype(dtypes_mapping)
combined_df.to_parquet(f'{write_path}/PIP_top_other.parquet', engine='pyarrow', compression='snappy',
                       partition_on=['chr'], overwrite=True)
print(f"{mega_cohort}: wrote {mega_dir}/PIP_all_other_parquet and PIP_top_other_parquet")
