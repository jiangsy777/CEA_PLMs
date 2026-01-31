import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# 1. 数据准备（不变）
data_dict = {
    'category': ['IAV_HA'] * 5 + ['Spike_RBD'] * 5,
    'Index': ['silhouette', 'CH', 'DBI', 'ARI', 'NMI'] * 2,
    'vBERT': [0.576, 3348.260, 0.580, 0.834, 0.905, 0.870, 14400.057, 0.215, 0.968, 0.955],
    'ESM-2 ft': [0.565, 2901.477, 0.608, 0.859, 0.908, 0.822, 5480.111, 0.384, 0.931, 0.920],
    'ESM-2': [0.562, 3064.684, 0.683, 0.906, 0.935, 0.581, 1165.237, 0.677, 0.771, 0.813],
}
df = pd.DataFrame(data_dict)
model_cols = ['vBERT', 'ESM-2 ft', 'ESM-2']
df_norm = df.copy()
for idx in df_norm.index:
    row_vals = df_norm.loc[idx, model_cols].values
    max_val = row_vals.max()
    df_norm.loc[idx, model_cols] = row_vals / max_val
dbi_mask = df_norm['Index'] == 'DBI'
df_norm.loc[dbi_mask, model_cols] = 1 - df_norm.loc[dbi_mask, model_cols] + 0.1
df_norm.loc[dbi_mask, 'Index'] = '1-DBI+0.1'

# 2. 统计数量（不变）
count_vbert_gt_esm2ft = (df_norm['vBERT'] > df_norm['ESM-2 ft']).sum()
count_esm2ft_gt_esm2 = (df_norm['ESM-2 ft'] > df_norm['ESM-2']).sum()

# 3. 绘图设置【核心修改区】
plt.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(25, 14), dpi=300)  # 画布适度放大适配大字号
gs = gridspec.GridSpec(2, 1, height_ratios=[4, 3], hspace=0.0)
ax_A = fig.add_subplot(gs[0])
ax_B = fig.add_subplot(gs[1], sharex=ax_A)

# 柱子更紧凑：total_bar_width调大到0.6，组内间距缩小，组间更密
colors = {'vBERT': '#d62728', 'ESM-2 ft': '#ffda24', 'ESM-2': '#1f77b4'}
model_order = ['vBERT', 'ESM-2 ft', 'ESM-2']
x = np.arange(len(df_norm))
total_bar_width = 0.6  # 关键：0.4→0.6，柱子更紧凑
single_bar_width = total_bar_width / len(model_order)
bar_offsets = {model: x - single_bar_width + i * single_bar_width for i, model in enumerate(model_order)}

# 4. A图：柱状图+大字号
for model in model_order:
    ax_A.bar(bar_offsets[model], df_norm[model], label=model, color=colors[model], width=single_bar_width,
             edgecolor='black', linewidth=0.3, zorder=3)
ax_A.set_ylim(0, 1.5)
ax_A.set_yticks(np.arange(0, 1.2, 0.2))
ax_A.set_ylabel('Normalized Value', fontsize=30)  # 放大y轴标签
ax_A.tick_params(axis='y', labelsize=26)  # 放大y轴刻度
ax_A.spines['top'].set_visible(False)
ax_A.spines['right'].set_visible(False)
ax_A.spines['left'].set_linewidth(0.8)
ax_A.spines['bottom'].set_linewidth(0.8)
ax_A.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
# 右上角注释放大
ax_A.text(0.98, 0.95, f'vBERT > ESM-2 ft: {count_vbert_gt_esm2ft} indicators', transform=ax_A.transAxes,
          ha='right', va='top', fontsize=32)
ax_A.text(0.98, 0.88, f'ESM-2 ft > ESM-2: {count_esm2ft_gt_esm2} indicators', transform=ax_A.transAxes,
          ha='right', va='top', fontsize=32)
# 图例放大
ax_A.legend(title='Model', loc='upper left', fontsize=28, title_fontsize=28,
            frameon=True, edgecolor='black', facecolor='white')
ax_A.set_xticks(x)
ax_A.tick_params(axis='x', labelbottom=False, length=0)
ax_A.set_xlim(-0.5, len(df_norm)-0.5)  # 适配紧凑柱子，居中不溢出

# 5. B图：指标+最优模型标注（大字号）
ax_B.set_xticks(x)
ax_B.set_xticklabels([])
ax_B.set_ylim(-3, 1)
ax_B.set_yticks([])
ax_B.spines[:].set_visible(False)
ax_B.tick_params(axis='x', length=0)

y_offset_counter = {}
for i, (_, row) in enumerate(df_norm.iterrows()):
    x_pos = i
    idx_name = row['Index']
    # 指标标签放大到14
    ax_B.text(x_pos, 0.4, idx_name, ha='center', va='bottom', rotation=0, fontsize=28)
    max_val = row[model_order].max()
    best_model = row[model_order][row[model_order] == max_val].index[0]
    best_color = colors[best_model]
    y_offset_counter[x_pos] = y_offset_counter.get(x_pos, 0)+1
    line_end_y = -0.8 if y_offset_counter[x_pos] % 2 == 1 else -2.2
    text_y = line_end_y - 0.15
    ax_B.plot([x_pos, x_pos], [0, line_end_y], color='black', linewidth=1, zorder=2)
    ax_B.scatter(x_pos, line_end_y, color=best_color, s=300, edgecolor='black', linewidth=0.3, zorder=4)
    # 最优模型标注放大到12
    ax_B.text(x_pos, text_y, best_model, ha='center', va='top', rotation=90, fontsize=28)

# 大类标签放大到14
iav_end = 4
ax_B.axvline(x=iav_end + 0.5, color='black', linewidth=0.8, linestyle='--', zorder=10)
# 最左端竖虚线：x 坐标设为 -0.5（与 x 轴左边界对齐）
ax_B.axvline(x=-0.5, color='black', linewidth=0.8, linestyle='--', zorder=10)
# 最右端竖虚线：x 坐标设为 len(df_norm)-0.5（与 x 轴右边界对齐）
ax_B.axvline(x=len(df_norm)-0.5, color='black', linewidth=0.8, linestyle='--', zorder=10)
ax_B.text(iav_end/2, -2.8, 'IAV_HA', ha='center', va='center', fontsize=30)
ax_B.text(iav_end + 0.5+(len(df_norm)-iav_end - 1)/2, -2.8, 'Spike_RBD',
          ha='center', va='center', fontsize=30)

# 6. 保存：bbox_inches='tight'确保完整显示大字号
plt.tight_layout()
plt.savefig('./generate/benchmarking_vBERT_esm2_ft_esm2_embedding_GIVAL_HA_S_RBDF_testset.png', dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.savefig('./generate/benchmarking_vBERT_esm2_ft_esm2_embedding_GIVAL_HA_S_RBDF_testset.svg', dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.show()