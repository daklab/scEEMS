"""
Create training datasets with stratified positive/negative sampling.

Samples positive variants (high PIP) and negative variants (low PIP) from
featurized gene datasets, stratified by variant type (SNP, insertion, deletion).

Usage:
    python create_training_datasets.py <chr_num> <data_split> <cohort> <NPR>

Arguments:
    chr_num: Chromosome number (1-22)
    data_split: "train" or "test"
    cohort: Cell type cohort (e.g., Mic_mega_eQTL)
    NPR: Number of negative samples per positive variant

Requires config.yaml and data_params.yaml.
"""

import os
import sys
import pandas as pd
import numpy as np
from dask import dataframe as dd
from tqdm import tqdm
import yaml
from dask.diagnostics import ProgressBar

ProgressBar().register()

# Get command line arguments
chr_num = sys.argv[1]
data_split = sys.argv[2]
cohort = sys.argv[3]
NPR = int(sys.argv[4])

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

params_data = yaml.safe_load(open('data_params.yaml'))[data_split]
pip_threshold_positive = params_data['positive_class_threshold']
pip_threshold_negative = params_data['negative_class_threshold']

# Define paths
data_dir = config['paths']['data_dir']
susie_dir = os.path.join(data_dir, 'susie_vars_pips')
featurized_dir = os.path.join(data_dir, f'training_data/{cohort}/all_variants')
top_parquet_dir = os.path.join(susie_dir, f'{cohort}/PIP_top_parquet')
all_parquet_dir = os.path.join(susie_dir, f'{cohort}/PIP_all_parquet')
output_dir = os.path.join(data_dir,
    f'training_data/{cohort}/training_data/{data_split}_NPR_{NPR}_PIP_{pip_threshold_positive}_{pip_threshold_negative}')
abc_data_dir = config['paths']['abc_data_dir']

os.makedirs(output_dir, exist_ok=True)

# Load gene TSS information
gene_id_gene_name_df = pd.read_csv(f'{abc_data_dir}/ABC_gene_id_name_mapping.csv', sep=',')

# Load data
top_pip_df_initial = dd.read_parquet(f'{top_parquet_dir}/PIP_top.parquet', engine='pyarrow')

if data_split == 'test':
    top_pip_df = top_pip_df_initial[top_pip_df_initial['pip'] >= pip_threshold_positive]
elif data_split == 'train':
    high_pip_genes = top_pip_df_initial[
        (top_pip_df_initial['pip'] > pip_threshold_positive) &
        (top_pip_df_initial['cs_coverage_0.95'] >= 1)
    ][['gene_id', 'cs_coverage_0.95']].compute()
    high_pip_genes = high_pip_genes.drop_duplicates()
    high_pip_no_coverage = top_pip_df_initial[
        (top_pip_df_initial['pip'] > 0.50) &
        (top_pip_df_initial['cs_coverage_0.95'] == 0)
    ].compute()
    covered_variants = top_pip_df_initial.merge(
        high_pip_genes, on=['gene_id', 'cs_coverage_0.95'], how='inner')
    covered_variants = covered_variants[covered_variants['pip'] > 0.05]
    top_pip_df = dd.concat([covered_variants, high_pip_no_coverage])

# Filter by chromosome
top_pip_df = top_pip_df[
    top_pip_df['chr'] == f'chr{chr_num}'
][['variant_id', 'chr', 'pos', 'ref', 'alt', 'pip', 'gene_id']].reset_index(drop=True)

# Load annotations
annotation_path = os.path.join(susie_dir, 'variant_list/annotated_variants_just_variants')
annotation_file = f'{annotation_path}/annotated_just_variants_chr{chr_num}.parquet'
annotation_df = dd.read_parquet(annotation_file, engine='pyarrow')

# Merge annotations
top_pip_df = top_pip_df.set_index('variant_id')
top_pip_df_annotated = top_pip_df.merge(
    annotation_df, left_index=True, right_index=True, how='inner'
).compute().reset_index()

positive_variants = top_pip_df_annotated[['variant_id', 'chr', 'pos', 'ref', 'alt', 'gene_id', 'pip']]


def get_variant_type(ref, alt):
    """Determine if variant is deletion (-1), SNP (0), or insertion (1)."""
    length_diff = len(alt) - len(ref)
    return -1 if length_diff < 0 else (1 if length_diff > 0 else 0)


def get_length_difference(ref, alt):
    """Calculate exact length difference between alleles."""
    return len(alt) - len(ref)


# Add variant type information
positive_variants['variant_type'] = positive_variants.apply(
    lambda row: get_variant_type(row['ref'], row['alt']), axis=1)
positive_variants['length_diff'] = positive_variants.apply(
    lambda row: get_length_difference(row['ref'], row['alt']), axis=1)
positive_variants['is_snp'] = (positive_variants['variant_type'] == 0).map({True: 1, False: 0})
positive_variants['var_type_snp'] = (positive_variants['variant_type'] == 0).map({True: 1, False: 0})
positive_variants['var_type_ins'] = (positive_variants['variant_type'] == 1).map({True: 1, False: 0})
positive_variants['var_type_del'] = (positive_variants['variant_type'] == -1).map({True: 1, False: 0})


def sample_variants(candidates, n_samples, var_type, gene_id):
    """Sample variants with replacement if necessary."""
    if len(candidates) == 0:
        print(f"Gene {gene_id}: No {var_type} candidates available")
        return pd.DataFrame()
    if n_samples == 0:
        return pd.DataFrame()
    if len(candidates) < n_samples:
        print(f"Gene {gene_id}: Not enough {var_type} candidates ({len(candidates)} available, "
              f"{n_samples} needed), sampling with replacement")
        return candidates.sample(n=n_samples, replace=True, random_state=42)
    else:
        return candidates.sample(n=n_samples, replace=False, random_state=42)


def get_stratified_negative_samples(gene_id, positive_vars, num_samples_per_variant, chr_str, pip_threshold=0.01):
    """Get negative examples stratified by variant type."""
    chr_int = int(chr_str.replace('chr', ''))
    annotation_file = f'{featurized_dir}/{gene_id}/annotated_data_{cohort}_chr{chr_int}.parquet'
    if not os.path.exists(annotation_file):
        print(f"Gene {gene_id}: No annotation file found")
        return pd.DataFrame()
    annotations = dd.read_parquet(annotation_file, engine='pyarrow').compute()
    valid_variant_ids = set(annotations['variant_id'].tolist())
    all_variants = dd.read_parquet(
        f'{all_parquet_dir}/PIP_all.parquet/chr=chr{chr_int}', engine='pyarrow')
    gene_variants_all = all_variants[
        (all_variants['gene_id'] == gene_id) &
        (all_variants['pip'] < pip_threshold)
    ][['variant_id', 'ref', 'alt', 'pos']].compute()
    gene_variants = gene_variants_all[gene_variants_all['variant_id'].isin(valid_variant_ids)].copy()
    if len(gene_variants) == 0:
        print(f"Gene {gene_id}: No negative samples available with annotations")
        return pd.DataFrame()
    gene_variants['variant_type'] = gene_variants.apply(
        lambda row: get_variant_type(row['ref'], row['alt']), axis=1)
    # Count positive variants by type
    pos_snp_count = sum(positive_vars['variant_type'] == 0)
    pos_ins_count = sum(positive_vars['variant_type'] == 1)
    pos_del_count = sum(positive_vars['variant_type'] == -1)
    needed_snps = pos_snp_count * num_samples_per_variant
    needed_ins = pos_ins_count * num_samples_per_variant
    needed_del = pos_del_count * num_samples_per_variant
    # Sample stratified by type
    snp_candidates = gene_variants[gene_variants['variant_type'] == 0]
    ins_candidates = gene_variants[gene_variants['variant_type'] == 1]
    del_candidates = gene_variants[gene_variants['variant_type'] == -1]
    sampled_snps = sample_variants(snp_candidates, needed_snps, "SNP", gene_id)
    sampled_ins = sample_variants(ins_candidates, needed_ins, "insertion", gene_id)
    sampled_del = sample_variants(del_candidates, needed_del, "deletion", gene_id)
    sampled_data = pd.concat([sampled_snps, sampled_ins, sampled_del], ignore_index=True)
    if len(sampled_data) == 0:
        return pd.DataFrame()
    result_data = pd.DataFrame({
        'variant_id': sampled_data['variant_id'],
        'gene_id': gene_id
    })
    result = result_data.merge(annotations, on='variant_id', how='inner')
    result = result.drop_duplicates(subset=['variant_id']).reset_index(drop=True)
    return result


def get_positive_features(variant_row):
    """Get features for a positive variant."""
    gene_id = variant_row['gene_id']
    chr_int = int(variant_row['chr'].replace('chr', ''))
    variant_id = variant_row['variant_id']
    annotation_file = f'{featurized_dir}/{gene_id}/annotated_data_{cohort}_chr{chr_int}.parquet'
    if not os.path.exists(annotation_file):
        return pd.DataFrame()
    annotations = dd.read_parquet(annotation_file, engine='pyarrow').compute()
    result = pd.DataFrame({'variant_id': [variant_id], 'gene_id': [gene_id]})
    result = result.merge(annotations, on='variant_id', how='inner')
    return result.drop_duplicates(subset=['variant_id']).reset_index(drop=True)


# Process negative examples
genes_to_process = positive_variants.groupby('gene_id')
print("Processing negative examples with stratified sampling...")
negative_data_list = []

for gene_id, gene_variants in tqdm(genes_to_process, desc="Processing genes"):
    negative_samples = get_stratified_negative_samples(
        gene_id=gene_id, positive_vars=gene_variants,
        num_samples_per_variant=NPR, chr_str=gene_variants['chr'].iloc[0],
        pip_threshold=pip_threshold_negative)
    if not negative_samples.empty:
        negative_data_list.append(negative_samples)

# Process positive examples
print("Processing positive examples...")
positive_data_list = []

for idx, row in tqdm(positive_variants.iterrows(), total=len(positive_variants)):
    positive_sample = get_positive_features(row)
    if not positive_sample.empty:
        positive_data_list.append(positive_sample)

if not negative_data_list or not positive_data_list:
    print("Error: Either no positive or no negative examples were found!")
    sys.exit(1)

negative_examples = pd.concat(negative_data_list)
variant_type_cols = ['length_diff', 'is_snp', 'var_type_snp', 'var_type_ins', 'var_type_del', 'variant_type']
cols_to_drop = [col for col in variant_type_cols if col in negative_examples.columns]
if cols_to_drop:
    negative_examples.drop(columns=cols_to_drop, inplace=True, errors='ignore')
negative_examples['label'] = 0

positive_examples = pd.concat(positive_data_list)
cols_to_drop = [col for col in variant_type_cols if col in positive_examples.columns]
if cols_to_drop:
    positive_examples.drop(columns=cols_to_drop, inplace=True, errors='ignore')
positive_examples['label'] = 1

# Combine and save
training_data = pd.concat([negative_examples, positive_examples], axis=0)
training_data = dd.from_pandas(training_data, npartitions=100)
training_data = training_data.repartition(partition_size="100MB")

output_file = f'{output_dir}/annotated_data_{cohort}_chr{chr_num}.parquet'
training_data.to_parquet(output_file, engine='pyarrow', compression='snappy', overwrite=True)

print(f"Successfully saved training data to {output_file}")
print(f"Total positive examples: {len(positive_examples)}")
print(f"Total negative examples: {len(negative_examples)}")
