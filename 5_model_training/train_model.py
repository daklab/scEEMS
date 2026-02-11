"""
Train CatBoost classifier for eQTL variant prediction.

Trains the conservative weighted CatBoost model using Leave-One-Chromosome-Out
cross-validation. This is the model variant used for final inference.

Usage:
    python train_model.py <cohort> <chromosome> --gene_lof_file <path> [--yaml_path <path>]

Arguments:
    cohort: Cell type cohort (e.g., Mic_mega_eQTL)
    chromosome: Held-out chromosome number for LOCO CV
    --gene_lof_file: Path to gene loss-of-function scores Excel file
    --yaml_path: Path to data_params.yaml (default: data_params.yaml)

Requires config.yaml with paths configured.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from os import walk
from sklearn import metrics
import yaml
import pickle
import random
import joblib
import dask

ProgressBar().register()

np.random.seed(9448)
random.seed(9448)

# Parse arguments
parser = argparse.ArgumentParser(description="Train eQTL prediction model")
parser.add_argument("cohort", type=str, help="Cohort name (e.g., Mic_mega_eQTL)")
parser.add_argument("chromosome", type=str, help="Held-out chromosome number (e.g., 2)")
parser.add_argument("--gene_lof_file", type=str, required=True,
                    help="Path to gene LoF scores Excel file")
parser.add_argument("--yaml_path", type=str, default="data_params.yaml",
                    help="Path to data_params.yaml")
parser.add_argument(
    "--single_chromosome_demo",
    action="store_true",
    help="Use only the selected chromosome for both train/test (public minimal demo mode).",
)
args = parser.parse_args()

cohort = args.cohort
chromosome = args.chromosome
gene_lof_file = args.gene_lof_file
params_data = yaml.safe_load(open(args.yaml_path))

# Load config
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

data_dir = config['paths']['data_dir']
scratch_dir = config['paths'].get('scratch_dir', '/tmp')
dask.config.set({'temporary_directory': scratch_dir})

chromosome_out = f'chr{chromosome}'
NPR_tr = 10
NPR_te = 10

chromosomes = [f'chr{x}' for x in range(1, 23)]
train_chromosomes = [x for x in chromosomes if x != chromosome_out]
test_chromosomes = [chromosome_out]
if args.single_chromosome_demo:
    train_chromosomes = [chromosome_out]
    test_chromosomes = [chromosome_out]
    print(f"Running in single-chromosome demo mode on {chromosome_out}.")

##############################################################################################################
# Load auxiliary data

gene_lof_df = pd.read_excel(gene_lof_file, "Supplementary Table 1")
gene_lof_df = gene_lof_df[['ensg', 'post_mean']]
gene_lof_df = gene_lof_df.rename(columns={'ensg': 'gene_id', 'post_mean': 'gene_lof'})
gene_lof_df['gene_lof'] = np.log2(gene_lof_df['gene_lof'])

gnomad_dir = config['paths']['gnomad_maf_dir']
maf_files = os.path.join(gnomad_dir, f'gnomad_MAF_chr{chromosome}.tsv')
maf_df = dd.read_csv(maf_files, sep='\t')
maf_df = maf_df[['variant_id', 'gnomad_MAF']].compute()

training_data_dir = os.path.join(data_dir, f'training_data/{cohort}')
write_dir = os.path.join(training_data_dir, 'model_results')
os.makedirs(write_dir, exist_ok=True)
os.makedirs(os.path.join(write_dir, 'predictions_parquet_catboost'), exist_ok=True)

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
    df['is_insertion'] = df['length_diff'].apply(lambda x: 1 if x > 0 else 0)
    df['is_deletion'] = df['length_diff'].apply(lambda x: 1 if x < 0 else 0)
    df.drop(columns=['chr', 'pos', 'ref', 'alt'], inplace=True)
    cols = df.columns.tolist()
    cols.insert(len(cols) - 1, cols.pop(cols.index('label')))
    df = df.loc[:, cols]
    return df


##############################################################################################################
# Load training data

train_files = []
for i in train_chromosomes:
    file_path = os.path.join(training_data_dir,
        f'training_data/train_NPR_{NPR_tr}_PIP_{params_data["train"]["positive_class_threshold"]}'
        f'_{params_data["train"]["negative_class_threshold"]}/annotated_data_{cohort}_{i}.parquet')
    if not os.path.exists(file_path):
        print(f"Warning: File {file_path} does not exist. Skipping chromosome {i}.")
        continue
    train_files.append(file_path)

if not train_files:
    raise ValueError("No training files found for any chromosome.")

train_df = dd.read_parquet(train_files, engine='pyarrow').compute()
train_df = make_variant_features(train_df)
train_df = train_df.merge(gene_lof_df, on='gene_id', how='left')
train_df = train_df.merge(maf_df, on='variant_id', how='left')
train_df['gene_lof'] = train_df['gene_lof'].fillna(train_df['gene_lof'].median())
train_df['gnomad_MAF'] = train_df['gnomad_MAF'].fillna(train_df['gnomad_MAF'].median())

##############################################################################################################
# Load test data

test_files = []
for i in test_chromosomes:
    file_path = os.path.join(training_data_dir,
        f'training_data/test_NPR_{NPR_te}_PIP_{params_data["test"]["positive_class_threshold"]}'
        f'_{params_data["test"]["negative_class_threshold"]}/annotated_data_{cohort}_{i}.parquet')
    if os.path.exists(file_path):
        test_files.append(file_path)

if not test_files:
    raise ValueError("No test files found.")

test_df = dd.read_parquet(test_files, engine='pyarrow').compute()
test_df = make_variant_features(test_df)
test_df = test_df.merge(gene_lof_df, on='gene_id', how='left')
test_df = test_df.merge(maf_df, on='variant_id', how='left')
test_df['gene_lof'] = test_df['gene_lof'].fillna(test_df['gene_lof'].median())
test_df['gnomad_MAF'] = test_df['gnomad_MAF'].fillna(test_df['gnomad_MAF'].median())

##############################################################################################################
# Calculate sample weights

train_class_0 = train_df[train_df['label'] == 0].shape[0]
train_total_pip = train_df[train_df['label'] == 1].pip.sum()
train_pip_percent = train_class_0 / train_total_pip if train_total_pip > 0 else 1
train_df['weight'] = np.where(train_df['label'] == 0, 1, train_df['pip'] * train_pip_percent)

test_class_0 = test_df[test_df['label'] == 0].shape[0]
test_total_pip = test_df[test_df['label'] == 1].pip.sum()
test_pip_percent = test_class_0 / test_total_pip if test_total_pip > 0 else 1
test_df['weight'] = np.where(test_df['label'] == 0, 1, test_df['pip'] * test_pip_percent)

##############################################################################################################
# Prepare features

meta_data = ['variant_id', 'pip', 'CHR', 'BP', 'REF', 'ALT', 'SNP', 'label', 'weight']

X_train = train_df.drop(columns=meta_data, errors='ignore')
Y_train = train_df['label']
weight_train = train_df['weight']
X_train = X_train.replace([np.inf, -np.inf], 0).fillna(0)

X_test = test_df.drop(columns=meta_data, errors='ignore')
Y_test = test_df['label']
X_test = X_test.replace([np.inf, -np.inf], 0).fillna(0)

if 'gene_id' in X_train.columns:
    X_train = X_train.drop(columns=['gene_id'])
    X_test = X_test.drop(columns=['gene_id'])

# Create feature subset
subset_keys = ['distance', 'ABC', 'celltype', 'baseline', 'chrombpnet_positive', 'diff', 'tf_positive']
subset_cols = []
for key in subset_keys:
    if key in column_dict:
        subset_cols.extend(column_dict[key])
subset_cols = [col for col in subset_cols if col in X_train.columns]

variant_features = ['length_diff', 'is_SNP', 'is_indel', 'is_insertion', 'is_deletion', 'gene_lof', 'gnomad_MAF']
subset_cols.extend(variant_features)

columns_to_abs = []
for key in ['diff', 'tf_positive', 'chrombpnet_positive']:
    if key in column_dict:
        columns_to_abs.extend([col for col in column_dict[key] if col in X_train.columns])

X_train_subset = X_train[subset_cols].copy()
X_test_subset = X_test[subset_cols].copy()

for col in columns_to_abs:
    if col in X_train_subset.columns:
        X_train_subset[col] = X_train_subset[col].abs()
        X_test_subset[col] = X_test_subset[col].abs()

X_train_subset = X_train_subset.drop(columns=['abs_distance_TSS', 'distance_TSS'], errors='ignore')
X_test_subset = X_test_subset.drop(columns=['abs_distance_TSS', 'distance_TSS'], errors='ignore')

##############################################################################################################
# Train CatBoost (conservative weighted)

from catboost import CatBoostClassifier

conservative_params = {
    'depth': 5,
    'iterations': 1000,
    'learning_rate': 0.03,
    'l2_leaf_reg': 5.0,
    'min_data_in_leaf': 10,
    'bagging_temperature': 1.0,
    'leaf_estimation_method': 'Newton',
    'leaf_estimation_iterations': 10,
    'verbose': True
}

# Create feature weight dictionary
feature_weights = {}
for col in X_train_subset.columns:
    feature_weights[col] = 1.0
for col in X_train_subset.columns:
    if col in columns_to_abs:
        if any(key in col for key in ['chrombpnet_positive', 'tf_positive', 'diff']):
            feature_weights[col] = 10.0

print(f"Features with weight 10.0: {sum(v == 10.0 for v in feature_weights.values())}")
print(f"Features with weight 1.0: {sum(v == 1.0 for v in feature_weights.values())}")

model = CatBoostClassifier(
    **conservative_params,
    loss_function='Logloss',
    feature_weights=feature_weights,
    name="Conservative-Weighted"
)

print("Training CatBoost model (conservative weighted)...")
model.fit(X_train_subset, Y_train, sample_weight=weight_train)

# Evaluate
preds = model.predict_proba(X_test_subset)[:, 1]
ap = metrics.average_precision_score(Y_test, preds)
auc = metrics.roc_auc_score(Y_test, preds)

print(f"\nTest Set Metrics - AP: {ap:.4f}, AUC: {auc:.4f}")

# Save model
model_file = os.path.join(write_dir,
    f'model_standard_subset_conservative_weighted_chr_{chromosome_out}_NPR_{NPR_tr}.joblib')
joblib.dump(model, model_file)
print(f"Model saved to: {model_file}")

# Save feature importances
importances = model.feature_importances_
features_df = pd.DataFrame({
    'feature': X_train_subset.columns,
    'importance': importances
}).sort_values(by='importance', ascending=False)
features_df.to_csv(os.path.join(write_dir,
    f'features_importance_chr_{chromosome_out}_NPR_{NPR_tr}.csv'), index=False)

# Save predictions
test_df['pred_prob'] = preds
test_df['pred_label'] = model.predict(X_test_subset)
test_df['actual_label'] = Y_test
test_df.to_csv(os.path.join(write_dir,
    f'predictions_parquet_catboost/predictions_chr{chromosome}.tsv'), sep='\t', index=False)

# Save summary
summary_dict = {
    'AP_test': ap,
    'AUC_test': auc,
    'params': conservative_params,
    'test_num_positive': int(Y_test.value_counts().get(1, 0)),
    'test_num_negative': int(Y_test.value_counts().get(0, 0)),
    'train_num_positive': int(Y_train.value_counts().get(1, 0)),
    'train_num_negative': int(Y_train.value_counts().get(0, 0))
}
with open(os.path.join(write_dir, f'summary_dict_chr_{chromosome_out}_NPR_{NPR_tr}.pkl'), 'wb') as f:
    pickle.dump(summary_dict, f)

print("\nTraining complete.")
