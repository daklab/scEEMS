"""
Compute SHAP values for high-confidence predictions.

Computes TreeSHAP explanations for variants with high prediction probabilities,
split into enhancer (>10kb from TSS) and promoter (<=10kb from TSS) categories.

Usage:
    python shap_analysis.py <cohort> <chromosome>

Arguments:
    cohort: Cell type cohort (e.g., Mic_mega_eQTL)
    chromosome: Chromosome number (1-22)

Requires config.yaml with paths configured.
"""

import os
import sys
import pandas as pd
import numpy as np
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from tqdm import tqdm
from os import walk
import pickle
import yaml
import joblib
import shap
import dask

ProgressBar().register()
tqdm.pandas()

np.random.seed(9448)

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

data_dir = config['paths']['data_dir']
scratch_dir = config['paths'].get('scratch_dir', '/tmp')
dask.config.set({'temporary_directory': scratch_dir})

params_data = yaml.safe_load(open('../5_model_training/data_params.yaml'))

cohort = sys.argv[1]
chr_num = int(sys.argv[2])
chromosome_out = f'chr{chr_num}'

pred_prob_filter = 0.95
NPR_tr = 10

training_data_dir = os.path.join(data_dir, f'training_data/{cohort}')
write_dir = os.path.join(training_data_dir, 'model_results')
variant_dir = os.path.join(training_data_dir, 'all_variants')

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


# Load predictions and filter to high-confidence
predictions_df = dd.read_parquet(
    os.path.join(training_data_dir, 'predictions_parquet_catboost/predictions.parquet'))
predictions_df = predictions_df[predictions_df['chr'] == chromosome_out]
predictions_df = predictions_df[predictions_df['pred_prob'] > pred_prob_filter]
predictions_df = predictions_df.compute()

# Load auxiliary data
gene_lof_file = config['paths']['gene_lof_file']
gene_lof_df = pd.read_excel(gene_lof_file, "Supplementary Table 1")
gene_lof_df = gene_lof_df[['ensg', 'post_mean']]
gene_lof_df = gene_lof_df.rename(columns={'ensg': 'gene_id', 'post_mean': 'gene_lof'})
gene_lof_df['gene_lof'] = np.log2(gene_lof_df['gene_lof'])

gnomad_dir = config['paths']['gnomad_maf_dir']
maf_df = dd.read_csv(os.path.join(gnomad_dir, f'gnomad_MAF_{chromosome_out}.tsv'), sep='\t')
maf_df = maf_df[['variant_id', 'gnomad_MAF']].compute()

# Load training data for feature alignment
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
X_train = train_df.drop(columns=meta_data, errors='ignore').replace([np.inf, -np.inf], 0).fillna(0)
if 'gene_id' in X_train.columns:
    X_train = X_train.drop(columns=['gene_id'])

subset_keys = ['distance', 'ABC', 'celltype', 'baseline', 'chrombpnet_positive', 'diff', 'tf_positive']
subset_cols = []
for key in subset_keys:
    if key in column_dict:
        subset_cols.extend(column_dict[key])
subset_cols = [col for col in subset_cols if col in X_train.columns]

variant_features_list = ['length_diff', 'is_SNP', 'is_indel', 'is_insertion', 'is_deletion', 'gene_lof', 'gnomad_MAF']
for feature in variant_features_list:
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

# Process high-confidence variants
gene_id_list = predictions_df['gene_id'].unique()
high_pip_df = predictions_df[['variant_id', 'gene_id']].copy(deep=True)

variant_features_combined = []
variant_distance_combined = []

for gene_id in tqdm(gene_id_list):
    variant_df = pd.read_parquet(os.path.join(variant_dir, gene_id))
    variant_df = variant_df.reset_index(drop=True)
    variant_df['gene_id'] = gene_id
    variant_df = make_variant_features(variant_df)
    variant_df = pd.merge(variant_df, gene_lof_df, on='gene_id', how='left')
    variant_df = pd.merge(variant_df, maf_df, on='variant_id', how='left')
    variant_df = variant_df.merge(high_pip_df, on=['gene_id', 'variant_id'], how='inner')
    variant_df['gene_lof'] = variant_df['gene_lof'].fillna(variant_df['gene_lof'].median())
    variant_df['gnomad_MAF'] = variant_df['gnomad_MAF'].fillna(variant_df['gnomad_MAF'].median())

    variant_distance_df = variant_df[['variant_id', 'distance_TSS']]
    variant_features_df = variant_df[subset_cols].copy()
    for col in columns_to_abs:
        if col in variant_features_df.columns:
            variant_features_df[col] = variant_features_df[col].abs()
    variant_features_df = variant_features_df.drop(columns=columns_to_drop, errors='ignore')
    variant_features_df = variant_features_df.reindex(columns=X_train_subset.columns)

    variant_features_combined.append(variant_features_df)
    variant_distance_combined.append(variant_distance_df)

variant_features_all = pd.concat(variant_features_combined, axis=0)
variant_distance_all = pd.concat(variant_distance_combined, axis=0).reset_index(drop=True)

# Compute SHAP values
from catboost import CatBoostClassifier

model_file = os.path.join(write_dir,
    f'model_standard_subset_conservative_weighted_chr_{chromosome_out}_NPR_10.joblib')
model = joblib.load(model_file)

explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(variant_features_all)
shap_values_df = pd.DataFrame(shap_values, columns=variant_features_all.columns).reset_index(drop=True)

variant_info_df_shap = pd.concat([variant_distance_all, shap_values_df], axis=1)

# Split by enhancer/promoter
enhancer_threshold = 10000
enhancer_shap = variant_info_df_shap[np.abs(variant_info_df_shap['distance_TSS']) > enhancer_threshold]
promoter_shap = variant_info_df_shap[np.abs(variant_info_df_shap['distance_TSS']) <= enhancer_threshold]

# Save
shap_dir = os.path.join(write_dir, 'shap_values')
os.makedirs(shap_dir, exist_ok=True)
enhancer_shap.to_csv(os.path.join(shap_dir, f'enhancer_shap_values_{chromosome_out}.csv'), index=False)
promoter_shap.to_csv(os.path.join(shap_dir, f'promoter_shap_values_{chromosome_out}.csv'), index=False)

print(f"SHAP analysis complete for {cohort} {chromosome_out}")
