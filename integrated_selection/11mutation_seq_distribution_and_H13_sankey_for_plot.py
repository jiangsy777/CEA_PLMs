import pandas as pd
import numpy as np

def parse_fasta(fasta_file):
    fr = open(fasta_file, 'r')
    contents = fr.readlines()
    idlst, seqlst, seq, seqnum = [], [], '', 0

    for line in contents:
        if line[:1] != '>':
            seq = seq + line[:-1]
        else:
            seqnum += 1
            seqlst.append(seq)
            seq = ''
            seqid = line[1:-1]
            idlst.append(seqid)
    seqlst = seqlst[1:]

    seqlast, line_num_all = '', len(contents)

    for j in range(line_num_all - 1, 0, -1):
        line = contents[j]
        if line[:1] != '>':
            seqlast = line[:-1] + seqlast
        else:
            break

    seqlst.append(seqlast)
    return idlst, seqlst

df = pd.read_csv('../data/df_HA_with_aa_seq_H13579_human_avian_without_wrong_seq_with_mafft_del_empty.csv')
df = df[df['Year']!='Unknown'].reset_index(drop=True)
df = df[(df['Continent']!='unknown')&(df['Continent']!='Unknown')].reset_index(drop=True)

area_lst = list(df['Continent'])
sero_lst = list(df['Serotype_H'])
#year_lst = list(df['Year'])
host_lst = list(df['Host1'])

id_lst1 = []
year_lst = []

for i in range(len(df)):
    id0 = df.loc[i,'strain_name']
    host0 = df.loc[i,'Host1']
    sero0 = df.loc[i,'Serotype_H']
    id_lst1.append(id0+'_'+host0+'_'+sero0)
    year = int(df.loc[i,'Year'])
    if year<2000:
        year_lst.append('before_2000')
    elif (year>=2000)&(year<2010):
        year_lst.append('2000-2009')
    elif (year>=2010)&(year<=2020):
        year_lst.append('2010-2020')
    elif year>2020:
        year_lst.append('after_2020')

dict_id_area = dict(zip(id_lst1,area_lst))
dict_id_sero = dict(zip(id_lst1,sero_lst))
dict_id_year = dict(zip(id_lst1,year_lst))
dict_id_host = dict(zip(id_lst1,host_lst))

df1 = pd.DataFrame()
fasta_file = '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta'
fasta1 = parse_fasta(fasta_file)
id_lst = fasta1[0]
seq_lst = fasta1[1]
df1['seq_id'] = id_lst
df1['seq'] = seq_lst
id_lst2 = []
seq_lst2 = []
area_lst2 = []
year_lst2 = []
host_lst2 = []
sero_lst2 = []
for i in range(len(df1)):
    print(i)
    seq_id = df1.loc[i,'seq_id']
    seq = df1.loc[i,'seq']

    if seq_id in id_lst1:
        id_lst2.append(seq_id)
        seq_lst2.append(seq)
        area_lst2.append(dict_id_area[seq_id])
        year_lst2.append(dict_id_year[seq_id])
        host_lst2.append(dict_id_host[seq_id])
        sero_lst2.append(dict_id_sero[seq_id])
df2 = pd.DataFrame()
df2['id'] = id_lst2
df2['seq'] = seq_lst2
df2['area'] = area_lst2
df2['year'] = year_lst2
df2['host'] = host_lst2
df2['Sero_H'] = sero_lst2
print(len(id_lst))
print(len(df2))
df2.to_csv('./generate/H13579_with_plosone_seq_and_annotations.csv',index=False)

in_lst = []
out_lst = []
for i in range(len(df2)):
    host = df2.loc[i,'host']
    sero = df2.loc[i,'Sero_H']
    area = df2.loc[i,'area']
    year = df2.loc[i,'year']
    in_lst.append(host)
    out_lst.append(sero)
    in_lst.append(sero)
    out_lst.append(area)
    in_lst.append(area)
    out_lst.append(year)
df3 = pd.DataFrame()
df3['in'] = in_lst
df3['out'] = out_lst
df3.to_csv('./generate/for_sankey_H13579_distribution.csv',index=False)

in_lst1 = []
out_lst1 = []

sero_lst = ['H1','H3']
sero_ref_lst = ['H1_ref','H3_ref']
dict_ref_lst = [{'H1':0,'H3':0,'H5':0,'H7':0,'H9':0} for a in range(len(sero_lst))]
for i in range(len(sero_lst)):
    sero = sero_lst[i]
    sero_ref = sero_ref_lst[i]

    df_sero = pd.read_csv('./final_selected/'+sero+'_vBERT_and_ESM2_ft_selected.csv')
    df_sero = df_sero.sort_values(by='site12_mutant_score')
    scores = df_sero['site12_mutant_score'].tolist()
    threshold_high = np.quantile(scores, 2 / 3)  # 66.6% 分位数
    df_sero = df_sero[df_sero['site12_mutant_score']>=threshold_high].reset_index(drop=True)
    df_sero['mut12'] = df_sero['site1_from1_plosone'].astype(str) + df_sero['site1_mut'].str.split('_').str[1]+'&'+df_sero['site2_from1_plosone'].astype(str) + df_sero['site2_mut'].str.split('_').str[1]

    df_sero['aa1'] = df_sero['site1_mut'].str.split('_').str[1]
    df_sero['aa2'] = df_sero['site2_mut'].str.split('_').str[1]
    for j in range(len(df_sero)):
        site1_from0 = int(df_sero.loc[j,'site1_from1_plosone']) - 1
        site2_from0 = int(df_sero.loc[j, 'site2_from1_plosone']) - 1
        aa1 = df_sero.loc[j,'aa1']
        aa2 = df_sero.loc[j, 'aa2']
        mut12 = df_sero.loc[j,'mut12']
        for k in range(len(df2)):
            seq = df2.loc[k,'seq']
            if (seq[site1_from0]==aa1)&(seq[site2_from0]==aa2):
                host = df2.loc[k, 'host']
                sero = df2.loc[k, 'Sero_H']
                area = df2.loc[k, 'area']
                year = df2.loc[k, 'year']

                dict_ref_lst[i][sero]+=1

                in_lst1.append(sero_ref)
                out_lst1.append(mut12)
                in_lst1.append(mut12)
                out_lst1.append(sero)
                in_lst1.append(sero)
                out_lst1.append(host)
                in_lst1.append(host)
                out_lst1.append(area)
                in_lst1.append(area)
                out_lst1.append(year)
df4 = pd.DataFrame()
df4['in'] = in_lst1
df4['out'] = out_lst1
df4.to_csv('./generate/for_sankey_H13_final_selected_mutations_distribution_high_group.csv',index=False)

df5 = pd.DataFrame()
df5['H1_ref'] = list(dict_ref_lst[0].values())
df5['H3_ref'] = list(dict_ref_lst[1].values())
# df5['H5_ref'] = list(dict_ref_lst[2].values())
df5.to_csv('./generate/for_duijizhuzhuangtu_H13_final_selected_mutations_distribution_high_group.csv',index=False)

