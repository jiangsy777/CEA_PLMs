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
df = df[df['Year'] != 'Unknown'].reset_index(drop=True)
df = df[(df['Continent'] != 'unknown') & (df['Continent'] != 'Unknown')].reset_index(drop=True)

area_lst = list(df['Continent'])
sero_lst = list(df['Serotype_H'])
# year_lst = list(df['Year'])
host_lst = list(df['Host1'])

id_lst1 = []
year_lst = []

for i in range(len(df)):
    id0 = df.loc[i, 'strain_name']
    host0 = df.loc[i, 'Host1']
    sero0 = df.loc[i, 'Serotype_H']
    id_lst1.append(id0 + '_' + host0 + '_' + sero0)
    year = int(df.loc[i, 'Year'])
    if year < 2000:
        year_lst.append('before_2000')
    elif (year >= 2000) & (year < 2010):
        year_lst.append('2000-2009')
    elif (year >= 2010) & (year <= 2020):
        year_lst.append('2010-2020')
    elif year > 2020:
        year_lst.append('after_2020')

dict_id_area = dict(zip(id_lst1, area_lst))
dict_id_sero = dict(zip(id_lst1, sero_lst))
dict_id_year = dict(zip(id_lst1, year_lst))
dict_id_host = dict(zip(id_lst1, host_lst))

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
    seq_id = df1.loc[i, 'seq_id']
    seq = df1.loc[i, 'seq']

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
df2.to_csv('./generate/H13579_with_plosone_seq_and_annotations.csv', index=False)

in_lst = []
out_lst = []
for i in range(len(df2)):
    host = df2.loc[i, 'host']
    sero = df2.loc[i, 'Sero_H']
    area = df2.loc[i, 'area']
    year = df2.loc[i, 'year']
    in_lst.append(host)
    out_lst.append(sero)
    in_lst.append(sero)
    out_lst.append(area)
    in_lst.append(area)
    out_lst.append(year)
df3 = pd.DataFrame()
df3['in'] = in_lst
df3['out'] = out_lst
df3.to_csv('./generate/for_sankey_H13579_distribution.csv', index=False)

in_lst1 = []
out_lst1 = []

sero_lst = ['H5', 'H7', 'H9']
sero_ref_lst = ['H5_ref', 'H7_ref', 'H9_ref']
# 创建三维字典：{sero: {year: count}}
dict_ref_lst = [{'H1': {}, 'H3': {}, 'H5': {}, 'H7': {}, 'H9': {}} for _ in range(len(sero_lst))]
for i in range(len(sero_lst)):
    sero = sero_lst[i]
    sero_ref = sero_ref_lst[i]

    df_sero = pd.read_csv('./final_selected/' + sero + '_vBERT_and_ESM2_ft_selected.csv')
    df_sero = df_sero.sort_values(by='site12_mutant_score')
    scores = df_sero['site12_mutant_score'].tolist()
    threshold_high = np.quantile(scores, 2 / 3)  # 66.6% 分位数
    df_sero = df_sero[df_sero['site12_mutant_score'] >= threshold_high].reset_index(drop=True)
    df_sero['mut12'] = df_sero['site1_from1_plosone'].astype(str) + df_sero['site1_mut'].str.split('_').str[1] + '&' + \
                       df_sero['site2_from1_plosone'].astype(str) + df_sero['site2_mut'].str.split('_').str[1]

    df_sero['aa1'] = df_sero['site1_mut'].str.split('_').str[1]
    df_sero['aa2'] = df_sero['site2_mut'].str.split('_').str[1]
    for j in range(len(df_sero)):
        site1_from0 = int(df_sero.loc[j, 'site1_from1_plosone']) - 1
        site2_from0 = int(df_sero.loc[j, 'site2_from1_plosone']) - 1
        aa1 = df_sero.loc[j, 'aa1']
        aa2 = df_sero.loc[j, 'aa2']
        mut12 = df_sero.loc[j, 'mut12']
        for k in range(len(df2)):
            seq = df2.loc[k, 'seq']
            if (seq[site1_from0] == aa1) & (seq[site2_from0] == aa2):
                host = df2.loc[k, 'host']
                sero = df2.loc[k, 'Sero_H']
                area = df2.loc[k, 'area']
                year = df2.loc[k, 'year']

                # 初始化年份计数
                if year not in dict_ref_lst[i][sero]:
                    dict_ref_lst[i][sero][year] = 0
                dict_ref_lst[i][sero][year] += 1

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
df4.to_csv('./generate/for_sankey_H579_final_selected_mutations_distribution_high_group.csv', index=False)

# 计算按年份的风险指数
all_risk_results = []
sero_set = ['H1', 'H3', 'H5', 'H7', 'H9']
year_groups = ['before_2000', '2000-2009', '2010-2020', 'after_2020']

for sero in sero_set:
    for year_group in year_groups:
        # 统计该亚型在该年份的序列数目
        sero_year_count = len(df2[(df2['Sero_H'] == sero) & (df2['year'] == year_group)])

        # 计算该亚型在该年份的突变总次数
        h5_ref_sum = dict_ref_lst[0][sero].get(year_group, 0)
        h7_ref_sum = dict_ref_lst[1][sero].get(year_group, 0)
        h9_ref_sum = dict_ref_lst[2][sero].get(year_group, 0)
        total_mutation_sum = h5_ref_sum + h7_ref_sum + h9_ref_sum

        # 计算风险指数
        if sero_year_count > 0:
            # 使用该年份组的序列数作为规模因子
            risk_index = (total_mutation_sum / sero_year_count) * np.log(sero_year_count)
        else:
            risk_index = 0

        all_risk_results.append({
            'Sero_H': sero,
            'Year_Group': year_group,
            'Risk_Index': risk_index,
            'Sequence_Count': sero_year_count,
            'Mutation_Count': total_mutation_sum
        })

# 创建DataFrame并保存
df_risk_year = pd.DataFrame(all_risk_results)
df_risk_year.to_csv('./generate/risk_index_by_year_H13579.csv', index=False)

# 创建宽表格式用于热图
df_risk_pivot = df_risk_year.pivot_table(index='Sero_H', columns='Year_Group', values='Risk_Index', aggfunc='first')
df_risk_pivot.to_csv('./generate/risk_index_by_year_pivot_H13579.csv')

print("按年份风险指数计算完成，已保存到 generate/risk_index_by_year_H13579.csv")
print(df_risk_year.head(20))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 读取风险指数数据
df_risk_year = pd.read_csv('./generate/risk_index_by_year_H13579.csv')

# 1. 统计所有数值的最大值，进行标准化处理
max_risk = df_risk_year['Risk_Index'].max()
min_risk = df_risk_year['Risk_Index'].min()
df_risk_year['Risk_Index_Normalized'] = df_risk_year['Risk_Index'] / max_risk*100 / 1.1

# 2. 调整年份分组顺序
year_order = ['before_2000', 'before_2020', '2000-2009', '2010-2020', 'after_2020']
df_risk_year['Year_Group'] = pd.Categorical(df_risk_year['Year_Group'], categories=year_order, ordered=True)
df_risk_year = df_risk_year.sort_values(['Sero_H', 'Year_Group'])

# 创建热图数据透视表
df_risk_pivot = df_risk_year.pivot_table(index='Sero_H', columns='Year_Group',
                                         values='Risk_Index_Normalized', aggfunc='first')

# 绘制热图
plt.figure(figsize=(12, 16))  # 竖长型图，高度大于宽度

# 使用低饱和度的colormap
sns.heatmap(df_risk_pivot,
            annot=True,
            fmt='.3f',
            cmap='coolwarm',
            cbar_kws={'label': 'Normalized Risk Index'},
            linewidths=0.5,
            linecolor='gray',
            vmin=min_risk,
            vmax=max_risk)

# 设置字体为Arial
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

# 设置标题和标签
plt.title('Risk Index by Serotype and Year Group (H13579)', fontsize=16, fontweight='bold')
plt.xlabel('Year Group', fontsize=14, fontweight='bold')
plt.ylabel('Serotype', fontsize=14, fontweight='bold')

# 设置边框粗细为0.5pt
ax = plt.gca()
ax.spines['top'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)

plt.tight_layout()
plt.savefig('./generate/risk_index_heatmap_by_year.svg', dpi=600, bbox_inches='tight')
plt.show()

# 输出标准化后的结果
print("标准化后的风险指数数据：")
print(df_risk_year[['Sero_H', 'Year_Group', 'Risk_Index', 'Risk_Index_Normalized']].head(20))