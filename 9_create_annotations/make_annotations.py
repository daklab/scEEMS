import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import dask
import yaml

import pyranges as pr




#abc_name = 'microglia'

# Set up progress bar
tqdm.pandas()
pbar = ProgressBar(dt=1)
pbar.register()

# Configure dask temporary directory
import datetime


def merge_baseline(baseline_df, predictions_df, column_name, cell_type, continuous_threshold=0.80, binary_threshold=0.90):
    baseline_df_initial = baseline_df.copy(deep=True)
    predictions_df_max = (predictions_df
                          .sort_values(['variant_id', column_name], ascending=[True, False])
                          .drop_duplicates('variant_id'))
    predictions_df_max_subset = predictions_df_max[['pos', column_name, 'ref', 'alt']]
    predictions_df_max_subset = predictions_df_max_subset.rename(columns={'pos': 'BP',
                                                                          column_name: f'{cell_type}_{column_name}'})
    predictions_df = baseline_df_initial.merge(predictions_df_max_subset, left_on=['BP', 'A1', 'A2'],
                                               right_on=['BP', 'ref', 'alt'], how='inner')
    # drop ref and alt
    predictions_df = predictions_df.drop(columns=['ref', 'alt'])
    predictions_df_flipped = baseline_df_initial.merge(predictions_df_max_subset, left_on=['BP', 'A1', 'A2'],
                                                       right_on=['BP', 'alt', 'ref'], how='inner')
    predictions_df_flipped = predictions_df_flipped.drop(columns=['ref', 'alt'])
    predictions_df_combined = pd.concat([predictions_df, predictions_df_flipped], axis=0)
    baseline_df_merged = pd.merge(baseline_df_initial, predictions_df_combined, on=['CHR', 'SNP', 'BP', 'A1', 'A2'],
                                  how='left')
    # fill NaN values with 0
    baseline_df_merged[f'{cell_type}_{column_name}'] = baseline_df_merged[f'{cell_type}_{column_name}'].fillna(0)
    # Create binary columns for continuous and binary thresholds
    baseline_df_merged[f'{cell_type}_{column_name}_binary'] = baseline_df_merged[f'{cell_type}_{column_name}'].apply(
        lambda x: 1 if x > binary_threshold else 0)
    # Apply linear scaling from continuous_threshold to 1
    baseline_df_merged[f'{cell_type}_{column_name}_continuous'] = baseline_df_merged[
        f'{cell_type}_{column_name}'].apply(lambda x: 0 if x < continuous_threshold else min(1, (x - continuous_threshold) / (1 - continuous_threshold))
                                            )
    #drop the original column
    baseline_df_merged = baseline_df_merged.drop(columns=[f'{cell_type}_{column_name}'])
    return(baseline_df_merged)


# Load config
with open(os.path.join(os.path.dirname(__file__), '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

dask.config.set({'temporary_directory': config['paths']['scratch_dir']})

chr = sys.argv[1]


baseline_dir = config['paths']['baseline_annot_dir']

baseline_file = f'{baseline_dir}/baseline_chr{chr}.annot.gz'

baseline_df = pd.read_csv(baseline_file, sep='\t',  usecols = ['CHR', 'SNP', 'BP', 'A1', 'A2'])

annotation_df = baseline_df.copy(deep=True)



gene_id_gene_name_dir = config['paths']['abc_data_dir']
gene_id_gene_name_df = pd.read_csv(f'{gene_id_gene_name_dir}/ABC_gene_id_name_mapping.csv', sep=',')



cell_type = sys.argv[2]

pred_prob_thresholds = [0.85, 0.86, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]

#cell_type = 'Mic_mega_eQTL'

out_dir = os.path.join(config['paths']['output_dir'], 'MLxQTL_pareto_nocontrol', cell_type) + '/'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)



for pred_prob_threshold in tqdm(pred_prob_thresholds):
# Load the predictions dataframe with high probability variants
    predictions_dir = os.path.join(config['paths']['data_dir'], cell_type)
    predictions_df = dd.read_parquet(f'{predictions_dir}/predictions_parquet_catboost_revised_nocontrol/predictions.parquet')
    predictions_df = predictions_df[predictions_df['chr'] == f'chr{chr}']
    predictions_df = predictions_df.persist()
    predictions_df_top_pred = predictions_df[predictions_df[f'pred_prob'] > pred_prob_threshold].compute()
    #rename pred_prob to pred_prob_{pred_prob_threshold}
    predictions_df_top_pred = predictions_df_top_pred.rename(columns={'pred_prob': f'{cell_type}_pred_prob_{pred_prob_threshold}',
                                                                      'pos': 'BP',
                                                                        'ref': 'A2',
                                                                        'alt': 'A1'}
                                                             )
    predictions_df_top_pred = predictions_df_top_pred[['BP', 'A1', 'A2', f'{cell_type}_pred_prob_{pred_prob_threshold}']]
    # get top prediction for each BP, A1, A2 combination
    predictions_df_top_pred = predictions_df_top_pred.groupby(['BP', 'A1', 'A2']).max().reset_index()
    predictions_df_top_pred[f'{cell_type}_pred_prob_{pred_prob_threshold}'] = 1
    #cast to integer
    predictions_df_top_pred[f'{cell_type}_pred_prob_{pred_prob_threshold}'] = predictions_df_top_pred[f'{cell_type}_pred_prob_{pred_prob_threshold}'].astype(int)
    pred_df = baseline_df.merge(predictions_df_top_pred, on = ['BP', 'A1', 'A2'], how='left')
    # fill NaN values with 0
    pred_df[f'{cell_type}_pred_prob_{pred_prob_threshold}'] = pred_df[f'{cell_type}_pred_prob_{pred_prob_threshold}'].fillna(0)
    annotation_df = pd.merge(annotation_df, pred_df, on=['CHR','SNP','BP', 'A1', 'A2'], how='left')
    annotation_df.head()






preamble_columns = ['CHR', 'SNP', 'BP', 'A1', 'A2']


# Save to file
annotation_df.to_csv(
    out_dir + 'MLxQTL_chr' + str(chr) + ".annot.gz",
    sep="\t",
    index=False,
    compression='gzip')

# Calculate the sum for the .l2.M file
sum_df = annotation_df.drop(columns=preamble_columns)
num_cols = sum_df.shape[1]
sum_np = sum_df.to_numpy()
sum = np.sum(sum_np, axis=0).reshape((1, num_cols))
np.savetxt(out_dir + 'MLxQTL_chr' + str(chr) + ".l2.M",
           sum,
           fmt='%f',
           delimiter=" ")


