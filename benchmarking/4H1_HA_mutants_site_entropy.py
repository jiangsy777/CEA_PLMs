import pandas as pd
import numpy as np

with open('./viruses-08-00155-s002/Supplemental_File_2_HApreferences.txt','r') as f1:
    lines = f1.readlines()
    # print(type(lines))
site_entropy_lst = []
site_AA_lst = []
site_lst = [a for a in range(2,566)]
for i in range(1,len(lines)):
    line = str(lines[i])
    line_split = line.split(' ')
    site_AA = str(line_split[1])
    site_entropy = line_split[2]
    site_AA_lst.append(site_AA)
    site_entropy_lst.append(site_entropy)
print(len(site_entropy_lst))
df_csv1 = pd.DataFrame()
df_csv1['site_from1'] = site_lst
df_csv1['site_entropy'] = site_entropy_lst
df_csv1['site_WT_AA'] = site_AA_lst
df_csv1.to_csv('./generate/HA_H1_site_entropy.csv',index=False)
