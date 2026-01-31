# -*- coding: utf-8 -*-
import math
import pickle
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# 设置全局字体为 Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']
# 关闭自动布局调整
plt.rcParams['figure.autolayout'] = False


################################
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


################################
def high_occur_pair_visualize(ax, path, sort_num, site1, site2, seqnum, RBS_name,
                              bg_color, bar_color, ss_color):
    colseq1 = colseq_lst[site1]
    colseq2 = colseq_lst[site2]
    thresh_freq = 0.01
    colseq_12 = [colseq1[i] + colseq2[i] for i in range(len(colseq1))]

    dict_Ccodonpairs = Counter(colseq_12)
    dict_Ccodonpairs = dict(sorted(dict_Ccodonpairs.items(), key=lambda x: x[1], reverse=True))

    pair_ref = list(dict_Ccodonpairs.keys())[0] if dict_Ccodonpairs else ""
    # title = f'Site {RBS_name} and {site2 + 1}; AA pair ref: \'{pair_ref}\''
    title = f'Site {RBS_name} and {site2 + 1}'

    dict_high_occur_pair_freq = {k: v for k, v in dict_Ccodonpairs.items() if v / seqnum > thresh_freq}
    names = list(dict_high_occur_pair_freq.keys())
    values = list(dict_high_occur_pair_freq.values())

    ax.set_facecolor(bg_color)
    ax.bar(names, values, color=bar_color)

    # 删除上、右边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # 设置坐标轴及刻度线粗细为0.5pt
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='both', width=0.5)

    # 字号适配，无加粗
    ax.set_title(title, fontsize=18, pad=5)
    ax.tick_params(axis='x', labelsize=16, rotation=90)
    ax.tick_params(axis='y', labelsize=17)

    # ========== 顶部横条：长度一致+宽度一致+不挡图 ==========
    rect = mpatches.Rectangle(
        (0, 0.95),  # 下边缘在子图95%高度处
        1.0,  # 长度=子图宽度100%，所有子图一致
        0.4,  # 宽度=子图高度的15%
        color=ss_color,
        transform=ax.transAxes  # 相对坐标保证全局尺寸统一
    )
    ax.add_patch(rect)
    # 扩展顶部空间，适配加宽的横条，不遮挡标题
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.4)


################################
def coevo_matrix_calculation(colseq1, colseq2, seqnum):
    thresh_freq = 0.01
    thresh_freq_high = 0.05
    colseq_12 = [colseq1[i] + colseq2[i] for i in range(len(colseq1))]
    colseq_1 = colseq1
    colseq_2 = colseq2

    dict_Ccodonpairs = Counter(colseq_12)
    pair_number0 = len(dict_Ccodonpairs)

    if pair_number0 == 1:
        coevo_score = 1
    else:
        max_count = max(dict_Ccodonpairs.values())
        pair_ref = [k for k, v in dict_Ccodonpairs.items() if v == max_count][0]

        two_sites_mutation_freq = 0
        for pair00 in dict_Ccodonpairs.keys():
            freq00 = dict_Ccodonpairs[pair00] / seqnum
            if (pair00[0] != pair_ref[0]) & (pair00[1] != pair_ref[1]) & (freq00 > thresh_freq):
                two_sites_mutation_freq += freq00

        if two_sites_mutation_freq == 0:
            coevo_score = 0
        else:
            site_pair_entropy = -sum([(v / seqnum) * math.log2(v / seqnum) for v in dict_Ccodonpairs.values()])
            dict_Csite1 = Counter(colseq_1)
            dict_Csite2 = Counter(colseq_2)

            co_mut_score = 1
            high_occur_pair_num = 0
            high_occur_pair_num1 = 0
            for pair000, count000 in dict_Ccodonpairs.items():
                freq000 = count000 / seqnum
                high_occur_pair_num1 += 1
                if freq000 > thresh_freq:
                    high_occur_pair_num += 1
                if freq000 > thresh_freq_high:
                    count_col1 = dict_Csite1[pair000[0]]
                    count_col2 = dict_Csite2[pair000[1]]
                    co_mut_score *= ((count000 ** 2) / (count_col1 * count_col2))
            coevo_score = (two_sites_mutation_freq ** 2 * co_mut_score ** 5) / (
                    math.sqrt(2 ** site_pair_entropy) * high_occur_pair_num * math.log2(high_occur_pair_num1))

    sitepair = (site1, site2)
    site_score = (sitepair, coevo_score)
    return site_score, dict_Ccodonpairs


################################
if __name__ == "__main__":
    file = '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta'
    path = ''

    # 解析fasta
    seq_id_lst, seq_seq_lst = parse_fasta(file)
    seqnum, colnum0 = len(seq_id_lst), len(seq_seq_lst[1])
    colseq_lst = []
    for j in range(colnum0):
        colseq = ''
        for k in range(seqnum):
            colseq += seq_seq_lst[k][j]
        colseq_lst.append(colseq)

    # 加载位点映射和二级结构
    df_aligned_true_site_ref = pd.read_csv('../data/generate/H1_dict_true_aligned_site.csv')
    dict_aligned_true_site_ref = dict(
        zip(df_aligned_true_site_ref['aligned_site'], df_aligned_true_site_ref['true_site']))
    df_true_site_ss_ref = pd.read_csv('../data/pdb_AlphaFold2/ref_H1_PlosOne_ref.csv')
    dict_true_site_ss_ref = dict(
        zip(df_true_site_ss_ref['Site_True_From_1'], df_true_site_ss_ref['Secondary_Structure']))

    # 定义RBS列表和对应RBD类型
    RBS_name_lst = [118, 123, 126, 116, 122, 132, 148, 147, 117, 149, 193, 189, 190, 228, 217, 224, 225, 226, 222,
                    213]
    RBD_type_lst = ['130loop'] * 10 + ['190helix'] * 3 + ['220loop'] * 7

    # ========== 扩展配色方案：7类二级结构 ==========
    color_config = {
        '130 loop': '#4B0082',
        '190 helix': '#8B0000',
        '220 loop': '#FF8C00',
        'Other helix': '#e59b3f',
        'Other sheet': '#d9d35a',
        'Other loop': '#62ac6c',
        'Other junction': '#2a736f'
    }

    # 拆分分组索引
    idx_130loop = list(range(10))
    idx_190helix = list(range(10, 13))
    idx_220loop = list(range(13, 20))

    # ========== 关键修改：使用GridSpec精确控制布局 ==========
    fig = plt.figure(figsize=(35, 35), constrained_layout=False)

    # 创建GridSpec，精确控制布局
    gs = fig.add_gridspec(
        nrows=5, ncols=5,
        left=0.15, right=0.95, top=0.90, bottom=0.10,
        wspace=0.35, hspace=1.5,
        width_ratios=[1] * 5, height_ratios=[1] * 5
    )

    # 创建子图
    axes = []
    for i in range(5):
        row_axes = []
        for j in range(5):
            ax = fig.add_subplot(gs[i, j])
            ax.axis('off')
            row_axes.append(ax)
        axes.append(row_axes)
    axes = np.array(axes)

    # 绘制子图
    for col, i in enumerate(idx_130loop[:5]):
        ax = axes[0, col]
        ax.axis('on')
        rbs = RBS_name_lst[i]
        rbd_type = RBD_type_lst[i]
        df_site1 = pd.read_csv('./generate/site_coevo_with_site' + str(rbs) + '.csv')
        site_20 = list(df_site1['site2_from1'])[1]

        # 映射真实位点并获取二级结构
        true_site = dict_aligned_true_site_ref.get(site_20, None)
        ss_type = dict_true_site_ss_ref.get(true_site, 'Junction')

        # 二级结构名称转换
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

        ss_color = color_config[ss_label]
        cfg_bg = '#E6E6FA' if rbd_type == '130loop' else '#FFE4E1' if rbd_type == '190helix' else '#FFEFD5'
        cfg_bar = '#4B0082' if rbd_type == '130loop' else '#8B0000' if rbd_type == '190helix' else '#FF8C00'

        high_occur_pair_visualize(ax, path, 1, rbs - 1, site_20 - 1, seqnum, rbs,
                                  cfg_bg, cfg_bar, ss_color)

    for col, i in enumerate(idx_130loop[5:]):
        ax = axes[1, col]
        ax.axis('on')
        rbs = RBS_name_lst[i]
        rbd_type = RBD_type_lst[i]
        df_site1 = pd.read_csv('./generate/site_coevo_with_site' + str(rbs) + '.csv')
        site_20 = list(df_site1['site2_from1'])[1]

        true_site = dict_aligned_true_site_ref.get(site_20, None)
        ss_type = dict_true_site_ss_ref.get(true_site, 'Junction')
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
        ss_color = color_config[ss_label]
        cfg_bg = '#E6E6FA'
        cfg_bar = '#4B0082'

        high_occur_pair_visualize(ax, path, 1, rbs - 1, site_20 - 1, seqnum, rbs,
                                  cfg_bg, cfg_bar, ss_color)

    for col, i in enumerate(idx_190helix):
        ax = axes[2, col]
        ax.axis('on')
        rbs = RBS_name_lst[i]
        rbd_type = RBD_type_lst[i]
        df_site1 = pd.read_csv('./generate/site_coevo_with_site' + str(rbs) + '.csv')
        site_20 = list(df_site1['site2_from1'])[1]

        true_site = dict_aligned_true_site_ref.get(site_20, None)
        ss_type = dict_true_site_ss_ref.get(true_site, 'Junction')
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
        ss_color = color_config[ss_label]
        cfg_bg = '#FFE4E1'
        cfg_bar = '#8B0000'

        high_occur_pair_visualize(ax, path, 1, rbs - 1, site_20 - 1, seqnum, rbs,
                                  cfg_bg, cfg_bar, ss_color)

    for col, i in enumerate(idx_220loop[:4]):
        ax = axes[3, col]
        ax.axis('on')
        rbs = RBS_name_lst[i]
        rbd_type = RBD_type_lst[i]
        df_site1 = pd.read_csv('./generate/site_coevo_with_site' + str(rbs) + '.csv')
        site_20 = list(df_site1['site2_from1'])[1]

        true_site = dict_aligned_true_site_ref.get(site_20, None)
        ss_type = dict_true_site_ss_ref.get(true_site, 'Junction')
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
        ss_color = color_config[ss_label]
        cfg_bg = '#FFEFD5'
        cfg_bar = '#FF8C00'

        high_occur_pair_visualize(ax, path, 1, rbs - 1, site_20 - 1, seqnum, rbs,
                                  cfg_bg, cfg_bar, ss_color)

    for col, i in enumerate(idx_220loop[4:]):
        ax = axes[4, col]
        ax.axis('on')
        rbs = RBS_name_lst[i]
        rbd_type = RBD_type_lst[i]
        df_site1 = pd.read_csv('./generate/site_coevo_with_site' + str(rbs) + '.csv')
        site_20 = list(df_site1['site2_from1'])[1]

        true_site = dict_aligned_true_site_ref.get(site_20, None)
        ss_type = dict_true_site_ss_ref.get(true_site, 'Junction')
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
        ss_color = color_config[ss_label]
        cfg_bg = '#FFEFD5'
        cfg_bar = '#FF8C00'

        high_occur_pair_visualize(ax, path, 1, rbs - 1, site_20 - 1, seqnum, rbs,
                                  cfg_bg, cfg_bar, ss_color)

    # ========== 左侧标签 ==========
    # fig.text(0.03, 0.95, 'Selected Site SS:', fontsize=18, ha='left', va='center')
    # fig.text(0.03, 0.76, '130 loop', fontsize=16, color='#4B0082', ha='left', va='center')
    # fig.text(0.03, 0.5, '190 helix', fontsize=16, color='#8B0000', ha='left', va='center')
    # fig.text(0.03, 0.23, '220 loop', fontsize=16, color='#FF8C00', ha='left', va='center')

    # ========== 生成7类图例 ==========
    legend_labels = ['130 loop', '190 helix', '220 loop',
                     'Other helix', 'Other sheet', 'Other loop', 'Other junction']
    legend_colors = [color_config[label] for label in legend_labels]
    patches = [mpatches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]

    # 添加图例，固定位置
    legend = fig.legend(
        handles=patches,
        loc='upper right',
        bbox_to_anchor=(0.97, 0.45),
        fontsize=18,
        title='Coevolution Site SS:',
        title_fontsize=18
    )

    # ========== 关键修改：保存前不调整布局 ==========
    # 不使用tight_layout，避免自动调整
    fig.set_constrained_layout(False)

    # 保存图片（不使用bbox_inches='tight'）
    folder_path = os.path.join(path, 'generate')
    os.makedirs(folder_path, exist_ok=True)

    # 保存时使用固定DPI和原始布局
    # save_path = os.path.join(folder_path, 'selected_site1_and_top1_coevo_site_AApair_distribution_and_site12_ss_marked_with_H1_ref.png')
    save_path = os.path.join(folder_path, 'selected_site1_and_top1_coevo_site_AApair_distribution_and_site12_ss_marked_with_H1_ref.svg')

    # 先显示查看效果
    plt.show()

    # 再保存（与显示时完全一致）
    fig.savefig(save_path, dpi=600, bbox_inches=None, pad_inches=0)

    print(f"图片已保存至: {save_path}")
    print("注意：保存的图片与显示效果完全一致，无布局变形")