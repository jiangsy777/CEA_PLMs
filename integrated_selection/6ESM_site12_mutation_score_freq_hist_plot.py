# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.patches import Rectangle  # 导入矩形绘制工具

# -------------------------- 1. 配置真实数据路径与参数 --------------------------
serotype_lst = ['H1', 'H3', 'H5','H7','H9']
csv_path_lst = [
    './generate/HA_H1_human_CA07_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H3_human_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H5_avian_Texas_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H7_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv',
    './generate/HA_H9_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv'
]
# 移除固定阈值列表，改为循环内计算中位数

# -------------------------- 2. 马卡龙色系配置 --------------------------
color_config = {
    "H1": {"bg": "#FFFFFF", "main": "#917bbd"},
    "H3": {"bg": "#FFFFFF", "main": "#F7B2C7"},
    "H5": {"bg": "#FFFFFF", "main": "#EEC78A"},
    "H7": {"bg": "#FFFFFF", "main": "#B0D6A9"},
    "H9": {"bg": "#FFFFFF", "main": "#7BC0CD"}
}

# -------------------------- 3. 全局绘图设置（学术标准） --------------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

# 创建画布：3个子图垂直排列，共享x轴+紧凑布局
fig, axes = plt.subplots(5, 1, figsize=(12, 10), sharex=True, sharey=True)
fig.subplots_adjust(hspace=0.15)  # 减小间距实现紧凑排列

# -------------------------- 4. 读取真实数据并绘制图表 --------------------------
for idx, sero in enumerate(serotype_lst):
    ax = axes[idx]
    csv_path = csv_path_lst[idx]

    # 读取CSV并提取有效mutant score
    df = pd.read_csv(csv_path)
    df = df[(df['delta_human_prob'] != 'single') & (df['delta_human_prob'].notna())].reset_index(drop=True)
    site12_mutant_score_lst = df['site12_mutant_score'].tolist()
    scores = np.array(site12_mutant_score_lst, dtype=np.float32)  # 转为数值数组
    print(len(scores))
    # 核心修改：计算当前亚型的中位数作为阈值
    # thresh = np.median(scores)
    thresh = np.quantile(scores, 0.7, method='linear')
    # 打印各亚型中位数阈值，方便查看
    print(f"【{sero} 亚型】 中位数阈值 = {thresh:.4f}")

    # 设置子图背景色
    ax.set_facecolor(color_config[sero]["bg"])

    # 绘制直方图：密度归一化，柱子边缘与主色一致
    n, bins, patches = ax.hist(scores, bins="auto", density=True,
                               facecolor=color_config[sero]["main"], alpha=0.6,
                               edgecolor=color_config[sero]["main"], linewidth=0.8)

    # 绘制核密度拟合曲线
    kde = gaussian_kde(scores)
    x_range = np.linspace(bins[0], bins[-1], 200)
    ax.plot(x_range, kde(x_range), color=color_config[sero]["main"], linewidth=2.5)

    # 绘制stretch阈值虚线：灰色虚线，仅在当前子图范围内
    ax.axvline(x=thresh, color="#888888", linestyle="--", linewidth=2, alpha=0.8,
               label=f"Threshold = {thresh:.2f}")

    # -------------------------- 核心保留：标签+小矩形色块 --------------------------
    # 1. 绘制小矩形色块（尺寸可通过width/height调整）
    rect = Rectangle((0.03, 0.85),  # 矩形左下角坐标（相对轴域）
                     width=0.02, height=0.04,  # 矩形宽高
                     color=color_config[sero]["main"],  # 矩形颜色=主色
                     transform=ax.transAxes)  # 相对轴域坐标
    ax.add_patch(rect)

    # 2. 添加黑色亚型文字（位置在矩形右侧）
    ax.text(0.06, 0.85, sero, fontsize=12, color="black",
            transform=ax.transAxes, verticalalignment="center")

    # 调整子图边框颜色与宽度
    for spine in ax.spines.values():
        # spine.set_color(color_config[sero]["main"])
        spine.set_color('black')
        spine.set_linewidth(0.5)

    # 调整刻度
    ax.tick_params(axis="both", labelsize=12, colors="#333333")
    ax.legend(loc="upper right", fontsize=12, frameon=True, facecolor="white", framealpha=0.8)

# -------------------------- 5. 共享坐标轴标签设置 --------------------------
# 仅在最下方子图添加x轴标签
axes[-1].set_xlabel("Site12 Mutant Score", fontsize=12)
# 仅在中间子图添加y轴标签并居中
axes[1].set_ylabel("Frequency Density", fontsize=12)
axes[1].yaxis.set_label_coords(-0.05, -0.7)

# -------------------------- 6. 保存高清学术图表 --------------------------
plt.savefig("./generate/ESM2_ft_score_hist_freq_H135_all.svg", dpi=600, bbox_inches="tight", facecolor="white")
# plt.savefig("ESM2_ft_score_hist_freq_H135_all.png", dpi=600, bbox_inches="tight", facecolor="white")
plt.show()