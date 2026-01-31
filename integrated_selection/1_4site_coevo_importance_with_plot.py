# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import mannwhitneyu

# 路径与参数
csv_path = "./generate/selected_site_coevo.csv"
# fig_save_path = "./coevo_site_importance_log2_normalized_1e-30_with_violin_selected_unselected.png"
fig_save_path = "./generate/coevo_site_importance_log2_normalized_1e-30_with_violin_selected_unselected.pdf"
# os.makedirs(os.path.dirname(fig_save_path), exist_ok=True)
#
regions = {
    "130 loop": {"range": (115, 149), "color": "#7E3A7E", "fill": "#D8BFD8", "error_fill": "#EBD8EB"},
    "190 Helix": {"range": (188, 195), "color": "#8B0000", "fill": "#F08080", "error_fill": "#FFB6C1"},
    "220 loop": {"range": (213, 229), "color": "#D2691E", "fill": "#FFE4B5", "error_fill": "#FFF0D6"}
}
other_color = "#008080"
other_error_fill = "#D1F0F0"

selected_sites = {
    "130 loop": [118,123,126,116,122,132,148,147,117,149],
    "190 Helix": [189, 190, 193],
    "220 loop": [228,217,224,225,226,222,213]
}
unselected_sites = {
    "130 loop": [124,133,146,134,128,135,144,115,127,137,143,129,131,145],
    "190 Helix": [],
    "220 loop": [214,223,221,227,215,220,229]
}

# 全局设置
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
font_size = 14
window_size = 5

# ========== 无截断的0-1标准化函数 ==========
def min_max_01(data):
    if len(data) == 0 or np.max(data) == np.min(data):
        return np.zeros_like(data) if len(data) > 0 else data
    return (data - np.min(data)) / (np.max(data) - np.min(data))

# 计算滚动标准差函数
def rolling_std(data, window=window_size):
    return pd.Series(data).rolling(window=window, min_periods=1, center=False).std()

# ---------------------- 数据处理：严格按你的顺序 ----------------------
df = pd.read_csv(csv_path)
df = df.sort_values(by='site1_from1').rename(columns={
    "site1_from1": "site",
    "coevo_mean_all_sites": "coevo_mean"
})
df = df[(df["site"] >= 70) & (df["site"] <= 290)].reset_index(drop=True)

# Step1: 所有数据先取log2
df["coevo_log2"] = np.log2(df["coevo_mean"] + 1e-30)

# Step2: 对log2结果做全局0-1标准化 → 数据天然在0-1
df["coevo_01"] = min_max_01(df["coevo_log2"].values)

# Step3: 标准化后的数据做5位点平滑 → 用于折线图
df_sorted = df.sort_values("site").reset_index(drop=True)
df_sorted["coevo_01_smooth"] = df_sorted["coevo_01"].rolling(
    window=window_size, min_periods=1, center=False
).mean()

# Step4: 计算滚动标准差（用于误差区间）
df_sorted["coevo_01_std"] = rolling_std(df_sorted["coevo_01"].values, window=window_size)

# Step5: 标记区域（用于确定不同区段的颜色）
def mark_region(row):
    for reg_name, reg_info in regions.items():
        if reg_info["range"][0] <= row["site"] <= reg_info["range"][1]:
            return reg_name
    return "other"
df_sorted["region"] = df_sorted.apply(mark_region, axis=1)

# 标记选中状态
df["region"] = df.apply(mark_region, axis=1)
df["is_selected"] = df.apply(
    lambda x: 1 if x["site"] in selected_sites.get(x["region"], []) else 0, axis=1
)

# ---------------------- 绘图：保证无数据污染 ----------------------
fig, ax_main = plt.subplots(figsize=(16, 9))
ax_main.set_position([0.1, 0.15, 0.8, 0.55])
ax_main.set_xlim(68, 292)
ax_main.set_ylim(0.25, 1.5)

# Step 1: 首先绘制整条连贯的折线（用统一的颜色）
ax_main.plot(
    df_sorted["site"], df_sorted["coevo_01_smooth"],
    color="#008080", linewidth=2.5, alpha=0.9, zorder=2, label="Smoothed Coevolution Score"
)

# Step 2: 绘制整条连贯的误差区间（±1标准差）
ax_main.fill_between(
    df_sorted["site"],
    df_sorted["coevo_01_smooth"] - df_sorted["coevo_01_std"],
    df_sorted["coevo_01_smooth"] + df_sorted["coevo_01_std"],
    color=other_error_fill, alpha=0.4, zorder=1, label="±1 SD"
)

# Step 3: 在不同区域上叠加对应颜色的误差区间
for reg_name, reg_info in regions.items():
    reg_df = df_sorted[df_sorted["region"] == reg_name]
    if len(reg_df) > 0:
        # 在该区域上叠加对应颜色的误差区间
        ax_main.fill_between(
            reg_df["site"],
            reg_df["coevo_01_smooth"] - reg_df["coevo_01_std"],
            reg_df["coevo_01_smooth"] + reg_df["coevo_01_std"],
            color=reg_info["error_fill"], alpha=0.5, zorder=1
        )

# Step 4: 在不同区域上叠加对应颜色的折线段
for reg_name, reg_info in regions.items():
    reg_df = df_sorted[df_sorted["region"] == reg_name]
    if len(reg_df) > 0:
        # 在该区域上叠加对应颜色的折线段
        ax_main.plot(
            reg_df["site"], reg_df["coevo_01_smooth"],
            color=reg_info["color"], linewidth=3, alpha=0.9, zorder=3
        )

# 小提琴图配置
violin_config = {
    "130 loop": {"ax_pos": [0.18, 0.75, 0.12, 0.2]},
    "190 Helix": {"ax_pos": [0.38, 0.75, 0.12, 0.2]},
    "220 loop": {"ax_pos": [0.58, 0.75, 0.12, 0.2]}
}

# P值标注函数
def add_p_value(ax, data1, data2, reg_color):
    if len(data1) == 0 or len(data2) == 0:
        return
    stat, p_val = mannwhitneyu(data1, data2, alternative="two-sided")
    p_text = "p < 0.001" if p_val < 0.001 else ("p < 0.01" if p_val < 0.01 else ("p < 0.05" if p_val < 0.05 else "p > 0.05"))
    ax.text(0.5, 1.3, p_text, ha="center", va="bottom", fontsize=font_size+2, color=reg_color, fontweight="bold")

# 绘制小提琴图：用【标准化后未平滑】的数据（天然0-1）
for reg_name in ["130 loop", "190 Helix", "220 loop"]:
    reg_info = regions[reg_name]
    # 取原始df中该区域的标准化未平滑数据
    reg_df = df[df["region"] == reg_name]
    data_selected = reg_df[reg_df["site"].isin(selected_sites[reg_name])]["coevo_01"].values
    data_unselected = reg_df[reg_df["site"].isin(unselected_sites[reg_name])]["coevo_01"].values

    ax_violin = fig.add_axes(violin_config[reg_name]["ax_pos"])
    ax_violin.set_xlim(-0.5, 1.5)  # 防止红叉偏移，不影响数值范围

    # 190 Helix 处理：红叉在x=0，小提琴在x=1，恢复深红色轮廓
    if len(data_unselected) == 0 and len(data_selected) > 0:
        sns.violinplot(
            data=[[np.nan], data_selected],
            ax=ax_violin,
            palette=["white", reg_info["fill"]],
            linewidth=2,
            edgecolor=[reg_info["color"], reg_info["color"]],  # 深红色轮廓
            inner=None
        )
        ax_violin.scatter(
            np.full_like(data_selected, 1), data_selected,
            s=35, alpha=0.7, color=reg_info["color"], zorder=3
        )
        ax_violin.text(0, np.mean(data_selected), "X", ha="center", va="center",
                       fontsize=22, color="red", fontweight="bold")
        ax_violin.set_xticks([0, 1])
        ax_violin.set_xticklabels(["Unselected", "Selected"], fontsize=font_size+2,rotation=45)

    else:
        sns.violinplot(
            data=[data_unselected, data_selected],
            ax=ax_violin,
            palette=["white", reg_info["fill"]],
            linewidth=2,
            edgecolor=[reg_info["color"], reg_info["color"]],
            inner=None
        )
        for i, d in enumerate([data_unselected, data_selected]):
            ax_violin.scatter(np.full_like(d, i), d, s=35, alpha=0.7, color=reg_info["color"], zorder=3)
        add_p_value(ax_violin, data_unselected, data_selected, reg_info["color"])
        ax_violin.set_xticks([0, 1])
        ax_violin.set_xticklabels(["Unselected", "Selected"], fontsize=font_size+2,rotation=45)

    # 小提琴图样式
    ax_violin.tick_params(axis="y", labelsize=font_size+2)
    ax_violin.set_ylabel("Normalized Score (log2)", fontsize=font_size+2, labelpad=3)
    ax_violin.spines["top"].set_visible(False)
    ax_violin.spines["right"].set_visible(False)

# 主图美化
ax_main.set_xlabel("Site", fontsize=font_size+4)
ax_main.set_ylabel("Normalized Coevolution Score (log2)", fontsize=font_size+4)
ax_main.tick_params(axis="both", labelsize=font_size+2)
ax_main.spines["top"].set_visible(False)
ax_main.spines["right"].set_visible(False)

# 创建区域颜色块的图例（用于指示不同功能区域）
from matplotlib.patches import Patch
legend_patches = [
    Patch(facecolor="#7E3A7E", alpha=0.7, label="130 loop (115-149)"),
    Patch(facecolor="#8B0000", alpha=0.7, label="190 Helix (188-195)"),
    Patch(facecolor="#D2691E", alpha=0.7, label="220 loop (213-229)"),
    Patch(facecolor="#008080", alpha=0.7, label="Other regions"),
    Patch(facecolor=other_error_fill, alpha=0.4, label="Error band (±1 SD)")
]
ax_main.legend(handles=legend_patches, loc="upper right", fontsize=font_size+3, frameon=True, fancybox=True, shadow=True)

plt.tight_layout()
plt.savefig(fig_save_path, dpi=600, bbox_inches="tight")
# plt.savefig(fig_save_path, dpi=600)

plt.show()
print(f"保存路径: {fig_save_path}")