# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ====================== 1. 配置参数（可自定义区域） ======================
# 7类二级结构的配色 (按你的要求定义，可直接修改)
ss_colors = {
    '130 loop': '#4B0082',
    '190 helix': '#8B0000',
    '220 loop': '#FF8C00',
    'Other helix': '#e59b3f',
    'Other sheet': '#d9d35a',
    'Other loop': '#62ac6c',
    'Other junction': '#2a736f'
}
ss_categories = ['130 loop', '190 helix', '220 loop', 'Other helix', 'Other sheet', 'Other loop', 'Other junction']
rbd_categories = ['130 loop', '190 helix', '220 loop']  # 3个柱子对应的分类

# 数据集配置 (填入你实际的文件路径)
datasets = [
    {
        'name': 'H1',
        'fasta': '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta',
        'csv': './generate/HA_H1_human_CA07_ref_site12_AA_pair_top3_site2_for_each_selected_site1.csv',
        'site_ref': '../data/generate/H1_dict_true_aligned_site.csv',
        'ss_ref': '../data/pdb_AlphaFold2/ref_H1_PlosOne_ref.csv',
        'RBS_name_lst': [118, 123, 126, 116, 122, 132, 148, 147, 117, 149, 193, 189, 190, 228, 217, 224, 225, 226,
                         222, 213],
        'RBD_type_lst': ['130loop'] * 10 + ['190helix'] * 3 + ['220loop'] * 7
    },
    {
        'name': 'H3',
        'fasta': '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta',
        'csv': "./generate/HA_H3_human_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1.csv",
        'site_ref': "../data/generate/H3_dict_true_aligned_site.csv",
        'ss_ref': "../data/pdb_AlphaFold2/ref_H3_PlosOne_ref.csv",
        'RBS_name_lst': [118, 123, 126, 116, 122, 132, 148, 147, 117, 149, 193, 189, 190, 228, 217, 224, 225, 226,
                         222, 213],
        'RBD_type_lst': ['130loop'] * 10 + ['190helix'] * 3 + ['220loop'] * 7
    },
    {
        'name': 'H5',
        'fasta': '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta',
        'csv': "./generate/HA_H5_avian_Texas_ref_site12_AA_pair_top3_site2_for_each_selected_site1.csv",
        'site_ref': "../data/generate/H5_dict_true_aligned_site.csv",
        'ss_ref': "../data/pdb_AlphaFold2/ref_H5_Texas_ref.csv",
        'RBS_name_lst': [118, 123, 126, 116, 122, 132, 148, 147, 117, 149, 193, 189, 190, 228, 217, 224, 225, 226,
                         222, 213],
        'RBD_type_lst': ['130loop'] * 10 + ['190helix'] * 3 + ['220loop'] * 7
    },

    {
        'name': 'H7',
        'fasta': '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta',
        'csv': "./generate/HA_H7_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1.csv",
        'site_ref': "../data/generate/H7_dict_true_aligned_site.csv",
        'ss_ref': "../data/pdb_AlphaFold2/ref_H7.csv",
        'RBS_name_lst': [118, 123, 126, 116, 122, 132, 148, 147, 117, 149, 193, 189, 190, 228, 217, 224, 225, 226,
                         222, 213],
        'RBD_type_lst': ['130loop'] * 10 + ['190helix'] * 3 + ['220loop'] * 7
    },
    {
        'name': 'H9',
        'fasta': '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta',
        'csv': "./generate/HA_H9_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1.csv",
        'site_ref': "../data/generate/H9_dict_true_aligned_site.csv",
        'ss_ref': "../data/pdb_AlphaFold2/ref_H7.csv",
        'RBS_name_lst': [118, 123, 126, 116, 122, 132, 148, 147, 117, 149, 193, 189, 190, 228, 217, 224, 225, 226,
                         222, 213],
        'RBD_type_lst': ['130loop'] * 10 + ['190helix'] * 3 + ['220loop'] * 7
    }
]

# 绘图参数
bar_width = 0.25  # 柱子宽度
pie_size_factor = 0.25  # 环形图大小
save_path = './fig/'  # 图片保存路径


# =========================================================================

# ====================== 2. 数据处理函数（与原热图一致） ======================
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
    fr.close()
    return idlst, seqlst


def get_freq_matrix(data_config):
    fasta_file = data_config['fasta']
    csv_path = data_config['csv']
    site_ref_csv = data_config['site_ref']
    ss_ref_csv = data_config['ss_ref']
    RBS_name_lst = data_config['RBS_name_lst']
    RBD_type_lst = data_config['RBD_type_lst']

    parse_fasta(fasta_file)
    site1_to_rbd = {site: rbd_type.replace('loop', ' loop').replace('helix', ' helix')
                    for site, rbd_type in zip(RBS_name_lst, RBD_type_lst)}

    df_aligned_true_site_ref = pd.read_csv(site_ref_csv)
    dict_aligned_true_site_ref = dict(
        zip(df_aligned_true_site_ref['aligned_site'], df_aligned_true_site_ref['true_site']))
    df_true_site_ss_ref = pd.read_csv(ss_ref_csv)
    dict_true_site_ss_ref = dict(
        zip(df_true_site_ss_ref['Site_True_From_1'], df_true_site_ss_ref['Secondary_Structure']))

    df_target0 = pd.read_csv(csv_path)
    df_target = pd.DataFrame()
    df_target['site1'] = list(df_target0['site1_from1_plosone'])[1:]
    df_target['site2'] = list(df_target0['site2_from1_plosone'])[1:]
    df_target['site_pair'] = df_target['site1'].astype(str) + '_' + df_target['site2'].astype(str)
    df_target_unique = df_target.drop_duplicates(subset=['site_pair'], keep='first').reset_index(drop=True)

    count_matrix = np.zeros((len(rbd_categories), len(ss_categories)), dtype=int)
    for _, row in df_target_unique.iterrows():
        site1 = row['site1']
        site2 = row['site2']
        rbd_label = site1_to_rbd[int(site1)]
        rbd_idx = rbd_categories.index(rbd_label)
        true_site = dict_aligned_true_site_ref.get(int(site2))
        ss_type = dict_true_site_ss_ref.get(true_site)

        if ss_type == 'Helix':
            ss_label = 'Other helix'
        elif ss_type == 'Sheet':
            ss_label = 'Other sheet'
        elif ss_type == 'Loop':
            ss_label = 'Other loop'
        elif ss_type == 'Junction':
            ss_label = 'Other junction'
        else:
            ss_label = ss_type

        if ss_label in ss_categories:
            ss_idx = ss_categories.index(ss_label)
            count_matrix[rbd_idx, ss_idx] += 1

    row_sums = count_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    freq_matrix = count_matrix / row_sums
    return freq_matrix


# ====================== 3. 绘图函数 ======================
def plot_combined_chart(data_config, freq_matrix):
    name = data_config['name']
    # 提取柱状图数据：每列对应一个二级结构，每行对应一个RBD分类
    bar_data = freq_matrix.T  # shape: (7, 3)

    # 创建画布
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(1, 1, left=0.1, right=0.9, top=0.85, bottom=0.2)
    ax_main = fig.add_subplot(gs[0])

    # X轴位置
    x_positions = np.arange(len(ss_categories))

    # 绘制3组柱状图
    for i, rbd in enumerate(rbd_categories):
        x = x_positions + (i - 1) * bar_width
        ax_main.bar(x, bar_data[:, i], width=bar_width, color=ss_colors[rbd],
                    edgecolor='black', linewidth=0.5)

    # 配置坐标轴
    ax_main.set_xlim(-0.5, len(ss_categories) - 0.5)
    ax_main.set_ylim(0, np.max(freq_matrix) * 1.2)
    ax_main.set_ylabel('Relative Proportion of Selected Site SS', fontsize=18, labelpad=15)
    ax_main.set_title(f'{name} Coevolution', fontsize=18, y=1.05)
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    # ax_main.spines['bottom'].set_linewidth(1.5)
    # ax_main.spines['left'].set_linewidth(1.5)
    # 设置坐标轴及刻度线粗细为0.5pt
    ax_main.spines['left'].set_linewidth(0.5)
    ax_main.spines['bottom'].set_linewidth(0.5)
    ax_main.tick_params(axis='both', width=0.5)
    ax_main.set_xticks(x_positions)
    ax_main.set_xticklabels(ss_categories, rotation=0, ha='center', fontsize=16)
    ax_main.tick_params(axis='y', labelsize=12)

    # 绘制底部色块图例
    rect_y_start = -np.max(freq_matrix) * 0.15
    rect_height = np.max(freq_matrix) * 0.05
    for i, ss in enumerate(ss_categories):
        rect_x_start = x_positions[i] - bar_width * 1.5
        rect_width = bar_width * 3
        rect = plt.Rectangle((rect_x_start, rect_y_start), rect_width, rect_height,
                             facecolor=ss_colors[ss], edgecolor='black', linewidth=0.15, clip_on=False)
        ax_main.add_patch(rect)

    # ========== 添加7类颜色的图例 ==========
    # 创建图例的子图区域
    ax_legend = fig.add_axes([0.65, 0.6, 0.2, 0.3])  # [左, 下, 宽, 高]
    ax_legend.axis('off')
    # 图例标题
    ax_legend.text(0.51, 0.97, 'Coevolution Site SS:', ha='center', va='top', fontsize=16)
    # 添加图例项
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=ss_colors[cat], edgecolor='black', linewidth=0.15)
                       for cat in ss_categories]
    ax_legend.legend(legend_elements, ss_categories, loc='center', ncol=1, frameon=True,
                     framealpha=1, edgecolor='black', fontsize=14, bbox_to_anchor=(0.5, 0.45))

    # ========== 添加环形图左侧的标注 ==========
    # 环形图位置参数
    pie_x_norm = [0.2, 0.35, 0.5]
    pie_y_norm = 0.75
    # 添加总标题
    fig.text(0.2, 0.87, 'Selected Site SS:', ha='right', va='center', fontsize=16)
    # 绘制3个环形图 (对应3行数据)
    for i, rbd in enumerate(rbd_categories):
        row_data = freq_matrix[i, :]
        ax_pie = fig.add_axes([pie_x_norm[i] - pie_size_factor / 2, pie_y_norm - pie_size_factor / 2,
                               pie_size_factor, pie_size_factor], zorder=10)
        ax_pie.pie(row_data, colors=[ss_colors[ss] for ss in ss_categories], radius=1, startangle=90,
                   wedgeprops=dict(width=0.4, edgecolor='white', linewidth=0.5))
        ax_pie.set_aspect('equal')
        ax_pie.text(0, 0, rbd, ha='center', va='center', fontsize=16)

    # 保存图片
    import os
    os.makedirs(save_path, exist_ok=True)
    plt.savefig(f'{save_path}/{name}_selected_site1_and_Top3_coevo_site2_each_serotype_selected_two_sites_mutation_not_empty_site.svg', dpi=600, bbox_inches='tight')
    # plt.savefig(
    #     f'{save_path}/{name}_selected_site1_and_Top3_coevo_site2_each_serotype_selected_two_sites_mutation_not_empty_site.png',
    #     dpi=600, bbox_inches='tight')

    plt.show()


# ====================== 4. 循环生成H1/H3/H5图表 ======================
for data in datasets:
    freq_matrix = get_freq_matrix(data)
    plot_combined_chart(data, freq_matrix)