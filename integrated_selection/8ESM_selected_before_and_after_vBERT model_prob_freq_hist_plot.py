# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, mannwhitneyu
from matplotlib.patches import Rectangle

# -------------------------- 1. 配置参数 --------------------------
serotype_lst = ['H1', 'H3', 'H5','H7','H9']
csv_path_lst = [
    './generate/HA_H1_human_CA07_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H3_human_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H5_avian_Texas_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H7_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H9_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv'
]
# thresh_lst = [-20, -15, -15]
column_lst = ['delta_avian_prob', 'delta_avian_prob', 'delta_human_prob', 'delta_human_prob', 'delta_human_prob']

# -------------------------- 2. 配色配置 --------------------------
color_config = {
    "H1": {"before": "#917bbd", "after": "#7B6C9F", "bg": "#FFFFFF"},
    "H3": {"before": "#F7B2C7", "after": "#F7A6AC", "bg": "#FFFFFF"},
    "H5": {"before": "#EEC78A", "after": "#FAA26F", "bg": "#FFFFFF"},
    "H7": {"before": "#B0D6A9", "after": "#65BDBA", "bg": "#FFFFFF"},
    "H9": {"before": "#7BC0CD", "after": "#3C9BC9", "bg": "#FFFFFF"}
}

# -------------------------- 3. 全局绘图设置 --------------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

# 创建画布：3行2列，共享坐标轴
fig, axes = plt.subplots(5, 2, figsize=(12, 10), sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.15, wspace=0.25)

# -------------------------- 4. 数据处理+绘图+检验 --------------------------
for row_idx, sero in enumerate(serotype_lst):
    csv_path = csv_path_lst[row_idx]
    # thresh = thresh_lst[row_idx]
    col = column_lst[row_idx]

    df = pd.read_csv(csv_path)
    df_before = df[(df['delta_human_prob'] != 'single') & (df['delta_human_prob'].notna())].reset_index(drop=True)
    df_before['site12_mutant_score'] = df_before['site12_mutant_score'].astype(float)
    scores = list(df_before['site12_mutant_score'])
    # thresh = np.median(scores)
    # thresh = np.quantile(scores,0.65,method='linear')
    thresh = np.quantile(scores, 0.7, method='linear')
    print(thresh)
    # 数据清洗
    vBERT_prob_before = df_before[col].astype(np.float64, errors='raise').values
    vBERT_prob_before = vBERT_prob_before[~np.isnan(vBERT_prob_before)]
    vBERT_prob_before = vBERT_prob_before[np.isfinite(vBERT_prob_before)]

    df_after = df_before[df_before['site12_mutant_score'] > thresh].reset_index(drop=True)
    vBERT_prob_after = df_after[col].astype(np.float64, errors='raise').values
    vBERT_prob_after = vBERT_prob_after[~np.isnan(vBERT_prob_after)]
    vBERT_prob_after = vBERT_prob_after[np.isfinite(vBERT_prob_after)]
    print(len(df_after))
    # 数据量校验
    if len(vBERT_prob_before) < 2 or len(vBERT_prob_after) < 2:
        print(f"警告：{sero} 组有效数据量不足，跳过该组绘图")
        continue

    # 曼-惠特尼U检验
    stat, p_val = mannwhitneyu(vBERT_prob_before, vBERT_prob_after, alternative='two-sided')
    p_text = f'P = {p_val:.4f}' if p_val >= 0.0001 else 'P < 0.0001'

    # 计算大于0的占比
    before_ratio = (vBERT_prob_before > 0).sum() / len(vBERT_prob_before) * 100
    after_ratio = (vBERT_prob_after > 0).sum() / len(vBERT_prob_after) * 100
    print((vBERT_prob_after > 0).sum())
    # 绘制 Before 子图
    ax_before = axes[row_idx, 0]
    ax_before.set_facecolor(color_config[sero]["bg"])
    n1, bins1, _ = ax_before.hist(vBERT_prob_before, bins="auto", density=True,
                                  facecolor=color_config[sero]["before"], alpha=0.6,
                                  edgecolor=color_config[sero]["before"], linewidth=0.8)
    kde1 = gaussian_kde(vBERT_prob_before)
    x_range1 = np.linspace(bins1[0], bins1[-1], 200)
    ax_before.plot(x_range1, kde1(x_range1), color=color_config[sero]["before"], linewidth=2.5)
    # x=0灰色竖线
    ax_before.axvline(x=0, color='#888888', linestyle='--', linewidth=1.5, alpha=0.8)
    # 占比标注
    ax_before.text(0.03, 0.65, f'> 0: {before_ratio:.1f}%', fontsize=12, color='black',
                   transform=ax_before.transAxes, verticalalignment='center')
    # 标签+边框
    rect_before = Rectangle((0.03, 0.83), 0.02, 0.04, color=color_config[sero]["before"], transform=ax_before.transAxes)
    ax_before.add_patch(rect_before)
    ax_before.text(0.06, 0.83, f'{sero} (Before)', fontsize=12, color='black', transform=ax_before.transAxes)
    for spine in ax_before.spines.values():
        # spine.set_color(color_config[sero]["before"])
        spine.set_color('black')
        spine.set_linewidth(0.5)
    ax_before.tick_params(colors="#333333")

    # 绘制 After 子图
    ax_after = axes[row_idx, 1]
    ax_after.set_facecolor(color_config[sero]["bg"])
    n2, bins2, _ = ax_after.hist(vBERT_prob_after, bins="auto", density=True,
                                 facecolor=color_config[sero]["after"], alpha=0.6,
                                 edgecolor=color_config[sero]["after"], linewidth=0.8)
    kde2 = gaussian_kde(vBERT_prob_after)
    x_range2 = np.linspace(bins2[0], bins2[-1], 200)
    ax_after.plot(x_range2, kde2(x_range2), color=color_config[sero]["after"], linewidth=2.5)
    # x=0灰色竖线
    ax_after.axvline(x=0, color='#888888', linestyle='--', linewidth=1.5, alpha=0.8)
    # 占比标注
    ax_after.text(0.03, 0.65, f'> 0: {after_ratio:.1f}%', fontsize=12, color='black',
                  transform=ax_after.transAxes, verticalalignment='center')
    # 标签+边框
    rect_after = Rectangle((0.03, 0.83), 0.02, 0.04, color=color_config[sero]["after"], transform=ax_after.transAxes)
    ax_after.add_patch(rect_after)
    ax_after.text(0.06, 0.83, f'{sero} (After)', fontsize=12, color='black', transform=ax_after.transAxes)

    # -------------------------- 核心修改：P值条件标注 --------------------------
    # if p_val < 0.05:
    #     ax_after.text(1.02, 0.5, p_text, fontsize=12, color='black',
    #                   transform=ax_after.transAxes, verticalalignment='center')

    for spine in ax_after.spines.values():
        spine.set_color('black')
        spine.set_linewidth(0.5)
    ax_after.tick_params(colors="#333333")

# -------------------------- 5. 坐标轴标签 --------------------------
for col_idx in range(2):
    axes[-1, col_idx].set_xlabel("vBERT Probability", fontsize=14)
# for row_idx in range(3):
#     axes[row_idx, 0].set_ylabel("Frequency Density", fontsize=12, fontweight='bold')
#     axes[row_idx, 0].yaxis.set_label_coords(-0.1, 0.5)
# 只在整个画布左侧垂直中间位置添加一次Y轴标签
fig.text(0.09, 0.5, "Frequency Density",
         fontsize=14,
         rotation=90,  # 文字垂直显示
         ha='center', # 水平居中
         va='center')  # 垂直居中

# -------------------------- 6. 保存 --------------------------
# plt.savefig("vBERT_prob_hist_freq_before_after_ESM2_ft_score_H135.png", dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig("./generate/vBERT_prob_hist_freq_before_after_ESM2_ft_score_H135.svg", dpi=600, bbox_inches="tight", facecolor="white")
plt.show()