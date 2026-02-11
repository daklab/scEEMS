import os
import sys
import pandas as pd
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
import yaml

# Load config
with open(os.path.join(os.path.dirname(__file__), '..', '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

pbar = ProgressBar(dt=1)
pbar.register()


cell_type = 'Mic_mega_eQTL'

ldsc_results = config['paths']['ldsc_results_dir']

ldscore_revised_dir = os.path.join(ldsc_results, 'MLxQTL_pareto', 'bellenguez_2022', cell_type)


ldscore_original_dir = os.path.join(ldsc_results, 'MLxQTL_pareto_original', 'bellenguez_2022', cell_type)


ldscore_unweighted_dir = os.path.join(ldsc_results, 'MLxQTL_pareto_unweighted', 'bellenguez_2022', cell_type)


ldscore_revised_df = dd.read_csv(f'{ldscore_revised_dir}/*.results', sep='\t')

ldscore_original_df = dd.read_csv(f'{ldscore_original_dir}/*.results', sep='\t')

ldscore_unweighted_df = dd.read_csv(f'{ldscore_unweighted_dir}/*.results', sep='\t')


# keep rows where Category starts with Mic_mega_eQTL_pred_prob
ldscore_revised_df = ldscore_revised_df[ldscore_revised_df['Category'].str.startswith(f'{cell_type}_pred_prob')].compute()

ldscore_original_df = ldscore_original_df[ldscore_original_df['Category'].str.startswith(f'{cell_type}_original_pred_prob')].compute()
ldscore_unweighted_df = ldscore_unweighted_df[ldscore_unweighted_df['Category'].str.startswith(f'{cell_type}_unweighted_pred_prob')].compute()


#scrub Mic_mega_eQTL_pred_prob_ from Category
ldscore_revised_df['pred_prob_pct'] = ldscore_revised_df['Category'].str.replace(f'{cell_type}_pred_prob_', '', regex=False)
ldscore_original_df['pred_prob_pct'] = ldscore_original_df['Category'].str.replace(f'{cell_type}_original_pred_prob_', '', regex=False)

ldscore_unweighted_df['pred_prob_pct'] = ldscore_unweighted_df['Category'].str.replace(f'{cell_type}_unweighted_pred_prob_', '', regex=False)




#scrub _1 from pred_prob_pct and conver to float
ldscore_revised_df['pred_prob_pct'] = ldscore_revised_df['pred_prob_pct'].str.replace('_1', '', regex=False).astype(float)
ldscore_original_df['pred_prob_pct'] = ldscore_original_df['pred_prob_pct'].str.replace('_1', '', regex=False).astype(float)
ldscore_unweighted_df['pred_prob_pct'] = ldscore_unweighted_df['pred_prob_pct'].str.replace('_1', '', regex=False).astype(float)


#sort by pred_prob_pct
ldscore_revised_df = ldscore_revised_df.sort_values(by='pred_prob_pct')
ldscore_original_df = ldscore_original_df.sort_values(by='pred_prob_pct')
ldscore_unweighted_df = ldscore_unweighted_df.sort_values(by='pred_prob_pct')


ldscore_revised_df[['Category', 'pred_prob_pct', 'Prop._SNPs', 'Prop._h2', 'Enrichment']]

ldscore_original_df[['Category', 'pred_prob_pct', 'Prop._SNPs', 'Prop._h2', 'Enrichment']]

ldscore_unweighted_df[['Category', 'pred_prob_pct', 'Prop._SNPs', 'Prop._h2', 'Enrichment']]


ldscore_revised_df['model'] = 'revised model'



ldscore_original_df['model'] = 'weighted'


ldscore_unweighted_df['model'] = 'unweighted'



ldscore_combined = pd.concat([ldscore_revised_df, ldscore_original_df, ldscore_unweighted_df], axis=0)




out_dir = os.path.join(config['paths']['output_dir'], 'aggregate_results')

ldscore_combined.to_csv(f'{out_dir}/ldscore_pareto_data.tsv', sep='\t', index=False)





cell_types = ['Mic_mega_eQTL', 'Ast_mega_eQTL', 'Oli_mega_eQTL', 'Exc_mega_eQTL', 'Inh_mega_eQTL', 'OPC_mega_eQTL']


ldscore_dataframes = []

for cell_type in cell_types:
    ldscore_original_dir = os.path.join(ldsc_results, 'MLxQTL_pareto_original', 'bellenguez_2022', cell_type)
    ldscore_original_df = dd.read_csv(f'{ldscore_original_dir}/*.results', sep='\t')
    ldscore_original_df = ldscore_original_df[ldscore_original_df['Category'].str.startswith(f'{cell_type}_original_pred_prob')].compute()
    ldscore_original_df['pred_prob_pct'] = ldscore_original_df['Category'].str.replace(f'{cell_type}_original_pred_prob_', '', regex=False)
    ldscore_original_df['pred_prob_pct'] = ldscore_original_df['pred_prob_pct'].str.replace('_1', '', regex=False).astype(float)
    ldscore_original_df = ldscore_original_df.sort_values(by='pred_prob_pct')
    ldscore_original_df['model'] = 'weighted'
    ldscore_original_df['cell_type'] = cell_type
    ldscore_dataframes.append(ldscore_original_df)

ldscore_combined_all = pd.concat(ldscore_dataframes, axis=0)

ldscore_combined_all.to_csv(f'{out_dir}/ldscore_pareto_data_all_cell_types.tsv', sep='\t', index=False)