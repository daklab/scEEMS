"""
Score all variants genome-wide using trained CatBoost models.

For each gene, loads the chromosome-matched LOCO model and predicts
the probability that each variant is a causal eQTL.

Usage:
    python model_inference.py <cohort> <gene_index>

Arguments:
    cohort: Cell type cohort (e.g., Mic_mega_eQTL)
    gene_index: 0-based index into the gene list

Requires config.yaml with paths configured.
"""

import os
import sys
import pandas as pd
import numpy as np
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from os import walk
import pickle
import yaml
import joblib
import dask

ProgressBar().register()

np.random.seed(9448)

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

data_dir = config['paths']['data_dir']
scratch_dir = config['paths'].get('scratch_dir', '/tmp')
dask.config.set({'temporary_directory': scratch_dir})

params_data = yaml.safe_load(open('../5_model_training/data_params.yaml'))

cohort = sys.argv[1]
idx = int(sys.argv[2])

NPR_tr = 10

training_data_dir = os.path.join(data_dir, f'training_data/{cohort}')
write_dir = os.path.join(training_data_dir, 'model_results')
variant_dir = os.path.join(training_data_dir, 'all_variants')
prediction_dir = os.path.join(training_data_dir, 'predictions_catboost')

os.makedirs(write_dir, exist_ok=True)
os.makedirs(os.path.join(write_dir, 'predictions_parquet_catboost'), exist_ok=True)
os.makedirs(prediction_dir, exist_ok=True)

# Load column dictionary
columns_dict_file = config['paths']['columns_dict_file']
with open(columns_dict_file, 'rb') as f:
    column_dict = pickle.load(f)


def make_variant_features(df):
    """Extract variant type features from variant_id."""
    df[['chr', 'pos', 'ref', 'alt']] = df['variant_id'].str.split(':', expand=True)
    df['length_diff'] = df['ref'].str.len() - df['alt'].str.len()
    df['is_SNP'] = df['length_diff'].apply(lambda x: 1 if x == 0 else 0)
    df['is_indel'] = df['length_diff'].apply(lambda x: 1 if x != 0 else 0)
    df['is_insertion'] = df['length_diff'].apply(lambda x: 1 if x < 0 else 0)
    df['is_deletion'] = df['length_diff'].apply(lambda x: 1 if x > 0 else 0)
    df.drop(columns=['chr', 'pos', 'ref', 'alt'], inplace=True, errors='ignore')
    return df


# Load gene list
list_genes_df = pd.read_csv(os.path.join(training_data_dir, 'list_genes.csv'), sep='\t')
gene_id = list_genes_df.iloc[idx]['gene_id']
chromosome_out = list_genes_df.iloc[idx]['chr']

# Load auxiliary data
gene_lof_file = config['paths']['gene_lof_file']
gene_lof_df = pd.read_excel(gene_lof_file, "Supplementary Table 1")
gene_lof_df = gene_lof_df[['ensg', 'post_mean']]
gene_lof_df = gene_lof_df.rename(columns={'ensg': 'gene_id', 'post_mean': 'gene_lof'})
gene_lof_df['gene_lof'] = np.log2(gene_lof_df['gene_lof'])

gnomad_dir = config['paths']['gnomad_maf_dir']
maf_files = os.path.join(gnomad_dir, f'gnomad_MAF_{chromosome_out}.tsv')
maf_df = dd.read_csv(maf_files, sep='\t')
maf_df = maf_df[['variant_id', 'gnomad_MAF']].compute()

# Load training data to get feature column alignment
chromosomes = [f'chr{x}' for x in range(1, 23)]
train_chromosomes = [x for x in chromosomes if x != chromosome_out]

train_files = []
for i in train_chromosomes:
    dir_path = os.path.join(training_data_dir,
        f'training_data/train_NPR_{NPR_tr}_PIP_{params_data["train"]["positive_class_threshold"]}'
        f'_{params_data["train"]["negative_class_threshold"]}/annotated_data_{cohort}_{i}.parquet')
    for (dirpath, dirnames, filenames) in walk(dir_path):
        for file in filenames:
            train_files.append(os.path.join(dir_path, file))

train_df = dd.read_parquet(train_files, engine='pyarrow').compute()
train_df = make_variant_features(train_df)
train_df = train_df.merge(gene_lof_df, on='gene_id', how='left')
train_df = train_df.merge(maf_df, on='variant_id', how='left')
train_df['gene_lof'] = train_df['gene_lof'].fillna(train_df['gene_lof'].median())
train_df['gnomad_MAF'] = train_df['gnomad_MAF'].fillna(train_df['gnomad_MAF'].median())

meta_data = ['variant_id', 'pip', 'CHR', 'BP', 'REF', 'ALT', 'SNP', 'label', 'weight']
X_train = train_df.drop(columns=meta_data, errors='ignore')
X_train = X_train.replace([np.inf, -np.inf], 0).fillna(0)
if 'gene_id' in X_train.columns:
    X_train = X_train.drop(columns=['gene_id'])

# Create feature subset
subset_keys = ['distance', 'ABC', 'celltype', 'baseline', 'chrombpnet_positive', 'diff', 'tf_positive']
subset_cols = []
for key in subset_keys:
    if key in column_dict:
        subset_cols.extend(column_dict[key])
subset_cols = [col for col in subset_cols if col in X_train.columns]

variant_features = ['length_diff', 'is_SNP', 'is_indel', 'is_insertion', 'is_deletion', 'gene_lof', 'gnomad_MAF']
for feature in variant_features:
    if feature in X_train.columns and feature not in subset_cols:
        subset_cols.append(feature)

columns_to_abs = []
for key in ['diff', 'tf_positive', 'chrombpnet_positive']:
    if key in column_dict:
        columns_to_abs.extend([col for col in column_dict[key] if col in X_train.columns])

X_train_subset = X_train[subset_cols].copy()
for col in columns_to_abs:
    if col in X_train_subset.columns:
        X_train_subset[col] = X_train_subset[col].abs()

columns_to_drop = ['abs_distance_TSS', 'distance_TSS']
X_train_subset = X_train_subset.drop(columns=columns_to_drop, errors='ignore')

##############################################################################################################
# Load and process variant data

variant_df = pd.read_parquet(os.path.join(variant_dir, gene_id))
variant_df = variant_df.reset_index(drop=True)
variant_df['gene_id'] = gene_id

variant_cols = ['variant_id', 'chr', 'pos', 'ref', 'alt', 'pip', 'gene_id']
variant_info_df = variant_df[variant_cols].copy()

variant_df = make_variant_features(variant_df)
variant_df = pd.merge(variant_df, gene_lof_df, on='gene_id', how='left')
variant_df = pd.merge(variant_df, maf_df, on='variant_id', how='left')
variant_df['gene_lof'] = variant_df['gene_lof'].fillna(variant_df['gene_lof'].median())
variant_df['gnomad_MAF'] = variant_df['gnomad_MAF'].fillna(variant_df['gnomad_MAF'].median())

variant_features_df_subset = variant_df[subset_cols].copy()
for col in columns_to_abs:
    if col in variant_features_df_subset.columns:
        variant_features_df_subset[col] = variant_features_df_subset[col].abs()
variant_features_df_subset = variant_features_df_subset.drop(columns=columns_to_drop, errors='ignore')
variant_features_df_subset = variant_features_df_subset.reindex(columns=X_train_subset.columns)

if X_train_subset.shape[1] != variant_features_df_subset.shape[1]:
    raise ValueError(f'Column mismatch: train={X_train_subset.shape[1]}, '
                     f'variant={variant_features_df_subset.shape[1]}')

##############################################################################################################
# Load model and predict

from catboost import CatBoostClassifier

model_file = os.path.join(write_dir,
    f'model_standard_subset_conservative_weighted_chr_{chromosome_out}_NPR_10.joblib')
model = joblib.load(model_file)

y_pred_proba = model.predict_proba(variant_features_df_subset)[:, 1]

variant_info_df['pred_prob'] = y_pred_proba
variant_info_df = variant_info_df.reset_index(drop=True)

output_file = os.path.join(prediction_dir, f'{gene_id}_predictions.tsv')
variant_info_df.to_csv(output_file, sep='\t', index=False)

print(f'Predictions for {gene_id} written to {output_file}')
