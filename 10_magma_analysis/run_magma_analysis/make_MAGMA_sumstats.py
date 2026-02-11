import os
import pandas as pd
import yaml

# Load config
with open(os.path.join(os.path.dirname(__file__), '..', '..', 'config.yaml'), 'r') as f:
    config = yaml.safe_load(f)

path_sumstats = config['paths']['sumstats_file']


#SNP CHR BP P Neff


write_sumstats_dir = config['paths']['magma_dir']


sumstats_df = pd.read_parquet(path_sumstats, engine='pyarrow')


sumstats_df_subset = sumstats_df[['SNP', 'CHR', 'BP', 'P', 'N']].copy(deep=True)


sumstats_df_subset.to_csv(f'{write_sumstats_dir}/bellenguez_MAGMA_sumstats.txt', sep=' ', header=True, index=False)



sumstats_dir = os.path.dirname(config['paths']['sumstats_file'])
path_sumstats_kunkle = os.path.join(sumstats_dir, 'Kunkle_et_al_2019_hg38.tsv.gz')




sumstats_df_kunkle = pd.read_csv(path_sumstats_kunkle, sep='\t', compression='gzip')


sumstats_df_kunkle_subset = sumstats_df_kunkle[['SNP', 'CHR', 'BP', 'P', 'N']].copy(deep=True)


sumstats_df_kunkle_subset.to_csv(f'{write_sumstats_dir}/kunkle_MAGMA_sumstats.txt', sep=' ', header=True, index=False)






path_sumstats_kunkle_AA = os.path.join(sumstats_dir, 'multi_ancestry', 'kunkle_AA_hg38.tsv.gz')




sumstats_df_kunkle_AFR = pd.read_csv(path_sumstats_kunkle_AA, sep='\t', engine='pyarrow', compression='gzip')


sumstats_df_kunkle_AFR_subset = sumstats_df_kunkle_AFR[['SNP', 'CHR', 'BP', 'P', 'N']].copy(deep=True)

#count rows with P < 1e-5
num_significant = sumstats_df_kunkle_AFR_subset[sumstats_df_kunkle_AFR_subset['P'] < 1e-5]

sumstats_df_kunkle_AFR_subset.to_csv(f'{write_sumstats_dir}/kunkle_AFR_MAGMA_sumstats.txt', sep=' ', header=True, index=False)







path_sumstats_ADSP_AFR = os.path.join(sumstats_dir, 'multi_ancestry', 'ADSP_AFR_hg38.tsv.gz')




sumstats_df_ADSP_AFR = pd.read_csv(path_sumstats_ADSP_AFR, sep='\t', engine='pyarrow', compression='gzip')


sumstats_df_ADSP_AFR_subset = sumstats_df_ADSP_AFR[['SNP', 'CHR', 'BP', 'P', 'N']].copy(deep=True)

#count rows with P < 1e-5

sumstats_df_ADSP_AFR_subset.to_csv(f'{write_sumstats_dir}/ADSP_AFR_MAGMA_sumstats.txt', sep=' ', header=True, index=False)








path_sumstats_ADSP_AMR = os.path.join(sumstats_dir, 'multi_ancestry', 'ADSP_AMR_hg38.tsv.gz')




sumstats_df_ADSP_AMR = pd.read_csv(path_sumstats_ADSP_AMR, sep='\t', engine='pyarrow', compression='gzip')


sumstats_df_ADSP_AMR_subset = sumstats_df_ADSP_AMR[['SNP', 'CHR', 'BP', 'P', 'N']].copy(deep=True)

#count rows with P < 1e-5

sumstats_df_ADSP_AMR_subset.to_csv(f'{write_sumstats_dir}/ADSP_AMR_MAGMA_sumstats.txt', sep=' ', header=True, index=False)
