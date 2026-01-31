# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.patches import Rectangle

# ===================== 全局配置（学术图表规范） =====================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.pad_inches'] = 0.1
plt.rcParams['figure.autolayout'] = False  # 禁用自动布局，避免错位

# 定义二级结构分类
SS_CATEGORIES = ['130 loop', '190 helix', '220 loop',
                 'Other helix', 'Other sheet', 'Other loop', 'Other junction']
# Site1分组配色
SITE1_GROUP_DICT = {
    '130 loop': {'range': (115, 149), 'color': '#7e3a7e'},
    '190 helix': {'range': (188, 195), 'color': '#c02942'},
    '220 loop': {'range': (213, 229), 'color': '#D2691E'}
}

# 热图参数
FILL_VALUE = -0.1
V_MIN = -0.1
HEATMAP_CMAP = 'Reds'
# 三级标签：high在上，low在下
LEVELS = ['high', 'middle', 'low']


def parse_fasta(fasta_file):
    fr = open(fasta_file, 'r')
    contents = fr.readlines()
    idlst, seqlst, seq, seqnum = [], [], '', 0
    for line in contents:
        if line[:1] != '>':
            seq += line[:-1]
        else:
            seqnum += 1
            seqlst.append(seq)
            seq = ''
            idlst.append(line[1:-1])
    seqlst = seqlst[1:]
    seqlast = ''
    for j in range(len(contents)-1, 0, -1):
        line = contents[j]
        if line[:1] != '>':
            seqlast = line[:-1] + seqlast
        else:
            break
    seqlst.append(seqlast)
    fr.close()
    return idlst, seqlst


def process_subtype_heatmap(subtype_config, csv_path):
    subtype = subtype_config['title'].split()[0]
    site_ref_file = subtype_config['site_ref']
    ss_ref_file = subtype_config['ss_ref']
    result_path = subtype_config['result_path']
    output_folder = os.path.join(result_path, f'generate')
    os.makedirs(output_folder, exist_ok=True)

    # 构建映射字典
    df_aligned_true = pd.read_csv(site_ref_file)
    dict_aligned2true = dict(zip(df_aligned_true['aligned_site'].astype(int),
                                 df_aligned_true['true_site'].astype(int)))
    df_true_ss = pd.read_csv(ss_ref_file)
    true_site_col = 'Site_True_From_1' if 'Site_True_From_1' in df_true_ss.columns else df_true_ss.columns[0]
    ss_col = 'Secondary_Structure' if 'Secondary_Structure' in df_true_ss.columns else df_true_ss.columns[1]
    dict_true2ss = dict(zip(df_true_ss[true_site_col].astype(int), df_true_ss[ss_col]))

    # 读取数据
    df = pd.read_csv(csv_path)
    df = df[(df['delta_human_prob'] != 'single') & (df['delta_human_prob'].notna())].reset_index(drop=True)
    df['site12_mutant_score'] = df['site12_mutant_score'].astype(float)
    scores = list(df['site12_mutant_score'])
    # thresh = np.median(scores)
    thresh = np.quantile(scores, 0.7, method='linear')
    df = df[df['site12_mutant_score'] > thresh].reset_index(drop=True)
    df['delta_avian_prob'] = df['delta_avian_prob'].astype(float)
    df = df[df['delta_avian_prob'] > 0].reset_index(drop=True)

    df['site1_from1_plosone'] = df['site1_from1_plosone'].astype(int)
    df['site2_from1_plosone'] = df['site2_from1_plosone'].astype(int)
    # df['delta_avian_prob'] = df['delta_avian_prob'].astype(float)

    # 分组聚合+按score升序排序
    df_grouped = df.groupby(['site1_from1_plosone', 'site2_from1_plosone'])['delta_avian_prob'].mean().reset_index()
    df_grouped = df_grouped.sort_values(by='delta_avian_prob', ascending=True).reset_index(drop=True)

    # 二级结构标签
    def get_ss_category(site2_aligned):
        try:
            true_site = dict_aligned2true[site2_aligned]
            ss = dict_true2ss[true_site]
            if ss == 'Helix':
                return 'Other helix'
            elif ss == 'Sheet':
                return 'Other sheet'
            elif ss == 'Loop':
                return 'Other loop'
            elif ss == 'Junction':
                return 'Other junction'
            else:
                return ss if ss in SS_CATEGORIES else 'Other loop'
        except KeyError:
            return 'Other loop'
    df_grouped['site2_ss'] = df_grouped['site2_from1_plosone'].apply(get_ss_category)

    # ========== 精确分位数阈值（33.3% 和 66.6%） ==========
    scores = df_grouped['delta_avian_prob'].values
    threshold_mid = np.quantile(scores, 1/3)    # 33.3% 分位数
    threshold_high = np.quantile(scores, 2/3)   # 66.6% 分位数

    def assign_level(score):
        if score >= threshold_high:
            return 0  # high
        elif score >= threshold_mid:
            return 1  # middle
        else:
            return 2  # low
    df_grouped['level'] = df_grouped['delta_avian_prob'].apply(assign_level)

    # 构建主热图矩阵
    df_grouped['x_label'] = df_grouped.apply(lambda row: f"{row['site1_from1_plosone']} & {row['site2_from1_plosone']}", axis=1)

    x_labels = df_grouped['x_label'].unique().tolist()
    N_X = len(x_labels)  # X轴格子数量
    N_Y_MAIN = len(SS_CATEGORIES)  # 主热图Y轴格子数量
    N_Y_LEVEL = len(LEVELS)  # 三级热图Y轴格子数量

    heatmap_matrix = np.full((N_Y_MAIN, N_X), FILL_VALUE, dtype=np.float_)
    for i, ss_cat in enumerate(SS_CATEGORIES):
        for j, x_lab in enumerate(x_labels):
            mask = (df_grouped['x_label'] == x_lab) & (df_grouped['site2_ss'] == ss_cat)
            if mask.any():
                heatmap_matrix[i, j] = df_grouped.loc[mask, 'delta_avian_prob'].values[0]

    # 构建三级热图矩阵
    level_matrix = np.full((N_Y_LEVEL, N_X), FILL_VALUE, dtype=np.float_)
    for j, x_lab in enumerate(x_labels):
        mask = df_grouped['x_label'] == x_lab
        if mask.any():
            level = df_grouped.loc[mask, 'level'].values[0]
            score = df_grouped.loc[mask, 'delta_avian_prob'].values[0]
            level_matrix[level, j] = score

    # ========== 核心修改1：调整子图位置右移，解决左边挡字、右边留白 ==========
    fig_width = N_X * 0.8
    fig_height = 6.5
    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=False)

    # 子图整体右移：left从0.1 → 0.15，width从0.8 → 0.75（保证总长度不变，右移后左边不挡字，右边留白减少）
    ax_main = fig.add_axes([0.15, 0.25, 0.75, 0.5])    # 主热图右移
    ax_level = fig.add_axes([0.15, 0.8, 0.623, 0.12]) # 顶部小热图同步右移，保持和主热图left/width一致

    # ========== 核心修改2：强制X轴范围一致，保证格子对齐 ==========
    x_min, x_max = -0.5, N_X - 0.5
    ax_main.set_xlim(x_min, x_max)
    ax_level.set_xlim(x_min, x_max)

    # 绘制三级小热图
    sns.heatmap(
        level_matrix,
        ax=ax_level,
        cmap=HEATMAP_CMAP,
        annot=False,
        cbar=False,
        xticklabels=False,
        yticklabels=LEVELS,
        linewidths=0.2,
        linecolor='gray',
        vmin=V_MIN,
        vmax=scores.max(),
        alpha=0.9,
        square=False
    )
    ax_level.tick_params(axis='y', labelsize=5, rotation=0)
    ax_level.set_ylabel('Level', fontsize=5, labelpad=2)
    ax_level.tick_params(color='gray')

    # 绘制主热图
    sns.heatmap(
        heatmap_matrix,
        ax=ax_main,
        cmap=HEATMAP_CMAP,
        annot=False,
        cbar=True,
        cbar_kws={'shrink': 0.5, 'pad': 0.02, 'location': 'right'},
        xticklabels=x_labels,
        yticklabels=SS_CATEGORIES,
        linewidths=0.2,
        linecolor='gray',
        vmin=V_MIN,
        vmax=scores.max(),
        alpha=0.9,
        square=False
    )
    ax_main.tick_params(color='gray')
    cbar = ax_main.collections[0].colorbar
    cbar.ax.tick_params(labelsize=4,color='gray')

    # 绘制顶部色带
    current_xlim = ax_main.get_xlim()
    top_band_colors = []
    for x_lab in x_labels:
        site1 = int(x_lab.split(' & ')[0])
        band_color = 'white'
        for group, info in SITE1_GROUP_DICT.items():
            if info['range'][0] <= site1 <= info['range'][1]:
                band_color = info['color']
                break
        top_band_colors.append(band_color)
    band_height = 0.2
    band_y_pos = 0
    bar_width = (current_xlim[1] - current_xlim[0]) / N_X
    for j in range(N_X):
        rect = Rectangle(
            (current_xlim[0] + j * bar_width, band_y_pos),
            bar_width, band_height,
            color=top_band_colors[j], linewidth=0.2, edgecolor='black', alpha=1,
            transform=ax_main.transData
        )
        ax_main.add_patch(rect)

    # 图例位置同步右移：left从0.75 → 0.8，和子图右移匹配
    legend_elements = [Rectangle((0, 0), 1, 1, facecolor=info['color'], label=group, edgecolor='black', linewidth=0.3)
                       for group, info in SITE1_GROUP_DICT.items()]
    ax_legend = fig.add_axes([0.76, 0.66, 0.12, 0.15])
    ax_legend.axis('off')
    ax_legend.legend(handles=legend_elements, loc='center', title='Selected Site SS',
                     title_fontsize=4, fontsize=4, frameon=True, framealpha=0.9)

    # 样式调整
    ax_main.set_xlabel('Coevolutionary Site1 & Site2', fontsize=5, labelpad=3)
    ax_main.set_ylabel('Coevolution Site SS', fontsize=5, labelpad=0.5)
    ax_main.set_yticklabels(SS_CATEGORIES, rotation=0, fontsize=5)
    ax_main.set_xticklabels(x_labels, rotation=90, ha='center', fontsize=5)
    fig.set_constrained_layout(False)
    fig.set_tight_layout(False)

    # ========== 核心修改3：保证保存和显示一致 ==========
    # 1. 强制渲染画布，让显示和保存的布局一致
    fig.canvas.draw()
    # 2. 保存时关闭bbox_inches='tight'，使用固定pad_inches，避免裁剪
    #保存会乱，直接存plt.show()的PDF版本
    # save_path = os.path.join(output_folder, f'{subtype}_heatmap_grid_aligned_final.png')
    # plt.savefig(
    #     save_path,
    #     bbox_inches=None,  # 禁用自动裁剪，和显示的画布完全一致
    #     pad_inches=0.1,
    #     dpi=600,
    #     facecolor='white',
    #     edgecolor='none'
    # )
    # print(f"✅ 图片已保存至: {save_path}")
    plt.show()


if __name__ == "__main__":
    INPUT_CSV = "./generate/HA_H3_human_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv"
    datasets = [
        {
            "title": "H3 Coevolution Score",
            "fasta": '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta',
            "site_ref": '../data/generate/H3_dict_true_aligned_site.csv',
            "ss_ref": '../data/pdb_AlphaFold2/ref_H3_PlosOne_ref.csv',
            "result_path": '.'
        }
    ]
    for ds in datasets:
        process_subtype_heatmap(ds, INPUT_CSV)