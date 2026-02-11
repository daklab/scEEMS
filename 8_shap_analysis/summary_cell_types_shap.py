"""
Cross-cell-type SHAP comparison.

Combines SHAP summaries across all cell types for comparative analysis.

Usage:
    python summary_cell_types_shap.py

Requires config.yaml with paths configured.
"""

import os
import pandas as pd
from tqdm import tqdm
import yaml

# Load configuration
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

data_dir = config['paths']['data_dir']
cohorts = config['cell_types']

aggregate_dir = os.path.join(data_dir, 'training_data/aggregate_results')
aggregate_shap_dir = os.path.join(aggregate_dir, 'shap_values')

shap_scores = [
    'enformer_shap_sum', 'chrombpnet_shap_sum', 'abc_score_shap_sum',
    'brain_CRE_shap_sum', 'distance_shap_sum', 'variant_shap_sum',
    'gene_shap_sum', 'baseline_shap_sum', 'TF_microglia_shap_sum',
    'TF_neuron_shap_sum', 'TF_oligodendrocyte_shap_sum',
    'TF_astrocyte_shap_sum', 'all_TF_shap_sum'
]

combined_df_list = []

for cohort in tqdm(cohorts):
    print(cohort)
    enhancer_summary = pd.read_csv(
        os.path.join(aggregate_shap_dir, f'{cohort}_enhancer_shap_values_summary.tsv'), sep='\t')
    promoter_summary = pd.read_csv(
        os.path.join(aggregate_shap_dir, f'{cohort}_promoter_shap_values_summary.tsv'), sep='\t')

    enhancer_subset = enhancer_summary[shap_scores].copy()
    promoter_subset = promoter_summary[shap_scores].copy()

    enhancer_subset['genomic_region'] = 'enhancer'
    promoter_subset['genomic_region'] = 'promoter'

    combined_df = pd.concat([enhancer_subset, promoter_subset], axis=0)
    combined_df = combined_df.groupby(['genomic_region']).mean().reset_index()
    combined_df['cohort'] = cohort
    combined_df_list.append(combined_df)

combined_df_all = pd.concat(combined_df_list, axis=0)
combined_df_all.to_csv(os.path.join(aggregate_shap_dir, 'all_cohorts_shap_summary.tsv'),
                       sep='\t', index=False)

print(f"Cross-cell-type SHAP summary saved")
