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

# Load config
with open(os.path.join(os.path.dirname(__file__), '..', '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

pbar = ProgressBar(dt=1)
pbar.register()



write_sumstats_dir = config['paths']['magma_dir']

ADGC_path = config['paths']['adgc_sumstats_dir']



#sumstats_bellenguez_df = pd.read_csv(f'{write_sumstats_dir}/bellenguez_MAGMA_sumstats.txt', sep=' ')

#sumstats_kunkle_AFR_df = pd.read_csv(f'{write_sumstats_dir}/kunkle_AFR_MAGMA_sumstats.txt', sep=' ')

#sumstats_ADSP_AFR_df = pd.read_csv(f'{write_sumstats_dir}/ADSP_AFR_MAGMA_sumstats.txt', sep=' ')


#sumstats_ADSP_AMR_df = pd.read_csv(f'{write_sumstats_dir}/ADSP_AMR_MAGMA_sumstats.txt', sep=' ')



#combined_df = pd.concat([sumstats_bellenguez_df[['SNP', 'CHR', 'BP']],
#                         sumstats_kunkle_AFR_df[['SNP', 'CHR', 'BP']],
#                         sumstats_ADSP_AFR_df[['SNP', 'CHR', 'BP']],
#                         sumstats_ADSP_AMR_df[['SNP', 'CHR', 'BP']]],
 #                       ignore_index=True)

#combined_df = combined_df.drop_duplicates().reset_index(drop=True)



plink_dir = config['paths']['plink_ref_dir']


pop='AFR'

AFR_bim = dd.read_csv(f'{plink_dir}/{pop}/plink/ADSP_{pop}_chr*.bim', sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'ALT', 'REF'])


#drop CM column
AFR_bim = AFR_bim.drop(columns=['CM'])

#AFR

ADGC_AFR_path = f'{ADGC_path}/AFA_common_apoe_adj_p-valueOnly.txt'

ADGC_AFR_df = pd.read_csv(ADGC_AFR_path, sep='\t')

ADGC_AFR_df[['CHR' , 'BP', 'REF', 'ALT']] = ADGC_AFR_df['MarkerName'].str.split(':', expand=True)

#rename P-value to P
ADGC_AFR_df = ADGC_AFR_df.rename(columns={'P-value': 'P'})

ADGC_AFR_df['CHR'] = ADGC_AFR_df['CHR'].apply(lambda x: x.replace('chr',''))

ADGC_AFR_df = ADGC_AFR_df.astype({"CHR": int, "BP": int})


ADGC_AFR_df = ADGC_AFR_df[['CHR', 'BP', 'P', 'REF', 'ALT']].sort_values(by=['CHR', 'BP']).reset_index(drop=True)

ADGC_AFR_df.shape

ADGC_AFR_df = ADGC_AFR_df.merge(AFR_bim.compute(), how = 'inner', on = ['CHR','BP','REF','ALT'])

ADGC_AFR_df.shape


ADGC_AFR_df['N'] = 6728


ADGC_AFR_df = ADGC_AFR_df[['SNP', 'CHR', 'BP', 'P', 'N']]


ADGC_AFR_df = ADGC_AFR_df.drop_duplicates(subset=['SNP'], keep='first').reset_index(drop=True)

ADGC_AFR_df.to_csv(f'{write_sumstats_dir}/ADGC_AFR_MAGMA_sumstats.txt', sep = ' ', header=True, index=False)




ADGC_AFR_df[ADGC_AFR_df['P'] < 0.005]




#AMR


pop='AMR'

AMR_bim = dd.read_csv(f'{plink_dir}/{pop}/plink/ADSP_{pop}_chr*.bim', sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'ALT', 'REF'])


#drop CM column
AMR_bim = AMR_bim.drop(columns=['CM'])


ADGC_AMR_path = f'{ADGC_path}/HISP_common_apoe_adj_p-valueOnly.txt'

ADGC_AMR_df = pd.read_csv(ADGC_AMR_path, sep='\t')

ADGC_AMR_df[['CHR' , 'BP', 'REF', 'ALT']] = ADGC_AMR_df['MarkerName'].str.split(':', expand=True)

#rename P-value to P
ADGC_AMR_df = ADGC_AMR_df.rename(columns={'P-value': 'P'})

ADGC_AMR_df['CHR'] = ADGC_AMR_df['CHR'].apply(lambda x: x.replace('chr',''))

ADGC_AMR_df = ADGC_AMR_df.astype({"CHR": int, "BP": int})


ADGC_AMR_df = ADGC_AMR_df [['CHR', 'BP', 'P', 'REF', 'ALT']].sort_values(by=['CHR', 'BP']).reset_index(drop=True)

ADGC_AMR_df.shape

ADGC_AMR_df = ADGC_AMR_df.merge(AMR_bim.compute(), how = 'inner', on = ['CHR','BP', 'REF', 'ALT'])

ADGC_AMR_df.shape

ADGC_AMR_df['N'] = 8899


ADGC_AMR_df = ADGC_AMR_df[['SNP', 'CHR', 'BP', 'P', 'N']]

ADGC_AMR_df = ADGC_AMR_df.drop_duplicates(subset=['SNP'], keep='first').reset_index(drop=True)


ADGC_AMR_df.to_csv(f'{write_sumstats_dir}/ADGC_AMR_MAGMA_sumstats.txt', sep = ' ', header=True, index=False)





#EAS


pop='EAS'

EAS_bim = dd.read_csv(f'{plink_dir}/{pop}/plink/ADSP_{pop}_chr*.bim', sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'ALT', 'REF'])


#drop CM column
EAS_bim = EAS_bim.drop(columns=['CM'])



ADGC_EAS_path = f'{ADGC_path}/EAS_common_apoe_adj_p-valueOnly.txt'

ADGC_EAS_df = pd.read_csv(ADGC_EAS_path, sep='\t')

ADGC_EAS_df[['CHR' , 'BP', 'REF', 'ALT']] = ADGC_EAS_df['MarkerName'].str.split(':', expand=True)

#rename P-value to P
ADGC_EAS_df = ADGC_EAS_df.rename(columns={'P-value': 'P'})

ADGC_EAS_df['CHR'] = ADGC_EAS_df['CHR'].apply(lambda x: x.replace('chr',''))

ADGC_EAS_df = ADGC_EAS_df.astype({"CHR": int, "BP": int})


ADGC_EAS_df = ADGC_EAS_df [['CHR', 'BP', 'P', 'REF', 'ALT']].sort_values(by=['CHR', 'BP']).reset_index(drop=True)

ADGC_EAS_df.shape

ADGC_EAS_df = ADGC_EAS_df.merge(EAS_bim.compute(), how = 'inner', on = ['CHR','BP', 'REF', 'ALT'])

ADGC_EAS_df.shape

ADGC_EAS_df['N'] = 3232


ADGC_EAS_df = ADGC_EAS_df[['SNP', 'CHR', 'BP', 'P', 'N']]


ADGC_EAS_df = ADGC_EAS_df.drop_duplicates(subset=['SNP'], keep='first').reset_index(drop=True)


ADGC_EAS_df.to_csv(f'{write_sumstats_dir}/ADGC_EAS_MAGMA_sumstats.txt', sep = ' ', header=True, index=False)


