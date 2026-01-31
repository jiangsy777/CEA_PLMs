import pandas as pd
import numpy as np
import scipy.stats as st

df_csv = pd.read_csv('df_DMS_mutants_H1_HA_with_preference_score.csv')
df_all_csv1 = pd.read_csv ('./generate/HA_H1_site_entropy.csv')

site_AA_WT_lst = list(df_all_csv1['site_WT_AA'])
H1_HA_DMS_WT_seq_from_site2 = ''.join(site_AA_WT_lst)
H1_HA_DMS_WT_seq_from_site1 = 'M' + H1_HA_DMS_WT_seq_from_site2
print(H1_HA_DMS_WT_seq_from_site1)
mutation_lst = []

for i in range(len(df_csv)):
    mut_site = int(df_csv.loc[i,'mutant_site'])-1
    ori_aa = H1_HA_DMS_WT_seq_from_site1[mut_site]
    mut_aa = df_csv.loc[i,'replace_amino']
    mutation = ori_aa+str(mut_site+1)+mut_aa
    mutation_lst.append(mutation)
df_csv['mutation'] = mutation_lst
df_csv.to_csv('./generate/HA_DMS_for_esm2_logist_matrix.csv',index=False)