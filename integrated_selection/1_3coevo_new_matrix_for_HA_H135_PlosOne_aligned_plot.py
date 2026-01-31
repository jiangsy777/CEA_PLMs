import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ----------------------------------------------------
# 1. 加载npy矩阵+精准匹配位点范围
# ----------------------------------------------------
npy_file_path = "./generate/coevo_matrix_H13579_PlosOne_70_290.npy"  # 替换为你的文件路径
coevo_matrix = np.load(npy_file_path)

# 真实位点范围（从1开始计数）
start_loc = 70
end_loc = 290
matrix_size = end_loc - start_loc + 1  # 计算位点总数

# 强制矩阵维度与位点数量一致
if coevo_matrix.shape[0] != matrix_size:
    print(f"警告：矩阵维度{coevo_matrix.shape}与位点数量{matrix_size}不匹配，已裁剪至{matrix_size}x{matrix_size}")
    coevo_matrix = coevo_matrix[:matrix_size, :matrix_size]

# 生成真实位点名称列表（70-290）
gene_names = list(range(start_loc, end_loc + 1))

# ----------------------------------------------------
# 2. 定义位点分组+配色
# ----------------------------------------------------
groups = {
    "130 loop": {"range": (115, 149), "color": "#7e3a7e"},
    "190 helix": {"range": (188, 195), "color": "#c02942"},
    "220 loop": {"range": (213, 229), "color": "#D2691E"},
    "Other": {"color": "#238b45"}
}

def get_site_color(site):
    for name, info in groups.items():
        if "range" in info and info["range"][0] <= site <= info["range"][1]:
            return info["color"]
    return groups["Other"]["color"]

color_list = [get_site_color(site) for site in gene_names]

# 十六进制转RGB（解决imshow报错）
def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return np.array([int(hex_color[i:i+2], 16)/255 for i in (0, 2, 4)])

rgb_color_list = [hex_to_rgb(c) for c in color_list]
color_array_top = np.array(rgb_color_list).reshape(1, matrix_size, 3)
color_array_right = np.array(rgb_color_list).reshape(matrix_size, 1, 3)

# ----------------------------------------------------
# 3. 全局字体设置
# ----------------------------------------------------
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
font_size = 12

# ----------------------------------------------------
# 4. 数据处理
# ----------------------------------------------------
coevo_matrix_log = np.sqrt(np.sqrt(np.sqrt(np.sqrt(np.sqrt(np.sqrt(coevo_matrix))))))
min_val = coevo_matrix_log.min()
max_val = coevo_matrix_log.max()

if max_val > min_val:
    data_normalized = (coevo_matrix_log - min_val) / (max_val - min_val)
else:
    data_normalized = coevo_matrix_log

data_matrix = pd.DataFrame(data_normalized, index=gene_names, columns=gene_names)

# ----------------------------------------------------
# 5. 布局设计：顶部色条精准对齐热图，右侧色条贴合
# ----------------------------------------------------
fig = plt.figure(figsize=(14, 12))
# 主热图轴
ax_main = fig.add_subplot(111)
# 顶部色条：嵌入主热图上方，与热图宽度完全一致
ax_top = ax_main.inset_axes([0, 1.005, 1, 0.03])  # [x,y,width,height] 精准对齐
# 右侧色条：嵌入主热图右侧，完全贴合
ax_right = ax_main.inset_axes([1.005, 0, 0.03, 1])

# ----------------------------------------------------
# 6. 绘制色条
# ----------------------------------------------------
ax_top.imshow(color_array_top, aspect="auto")
ax_top.axis("off")

ax_right.imshow(color_array_right, aspect="auto")
ax_right.axis("off")

# ----------------------------------------------------
# 7. 绘制主热图：调整蓝红色条缩短+右移，避免重叠
# ----------------------------------------------------
im = ax_main.imshow(data_matrix.values, cmap='RdBu_r', vmin=0, vmax=1, aspect='equal')

# 核心调整：shrink从0.8→0.6 缩短色条，pad=0.12增大间距
cbar = fig.colorbar(im, ax=ax_main, shrink=0.5, pad=0.04, location='right')
cbar.set_label('Normalized Score', rotation=270, labelpad=20, fontsize=font_size)
cbar.ax.tick_params(labelsize=font_size)

# ----------------------------------------------------
# 8. 刻度设置：精准匹配位点
# ----------------------------------------------------
tick_interval = 20
xticks_positions = []
xticks_labels = []
yticks_positions = []
yticks_labels = []

for i, site in enumerate(gene_names):
    if site % tick_interval == 0 or site == gene_names[0] or site == gene_names[-1]:
        xticks_positions.append(i)
        xticks_labels.append(str(site))
        yticks_positions.append(i)
        yticks_labels.append(str(site))

ax_main.set_xticks(xticks_positions)
ax_main.set_xticklabels(xticks_labels, fontsize=font_size, rotation=90)
ax_main.set_yticks(yticks_positions)
ax_main.set_yticklabels(yticks_labels, fontsize=font_size)

ax_main.set_xlabel('Site (70-290)', fontsize=font_size+1, labelpad=10)
ax_main.set_ylabel('Site (70-290)', fontsize=font_size+1, labelpad=10)
ax_main.set_xlim(-0.5, matrix_size-0.5)
ax_main.set_ylim(matrix_size-0.5, -0.5)

# ----------------------------------------------------
# 9. 美化图例：上移调整 bbox_to_anchor 从(1.32, 0.95)→(1.32, 1.02)
# ----------------------------------------------------
legend_elements = [
    plt.Rectangle((0, 0), 1, 1, facecolor=info["color"], label=name, edgecolor='gray', linewidth=0.5)
    for name, info in groups.items()
]
legend = ax_main.legend(
    handles=legend_elements,
    bbox_to_anchor=(1.25, 1.02),  # 上移到色条上方
    loc="upper right",
    title="Site Groups",
    title_fontsize=font_size+1,
    fontsize=font_size,
    frameon=True,        # 显示边框
    fancybox=True,       # 圆角边框
    shadow=False,        # 无阴影
    borderpad=0.8,       # 内边距
    framealpha=0.9       # 边框透明度
)

# ----------------------------------------------------
# 10. 保存图片
# ----------------------------------------------------
# os.makedirs('./result/fig', exist_ok=True)
# plt.savefig('./coevo_H13579_PlosOne_sqrt16_70_290.png', dpi=600, bbox_inches='tight')
plt.savefig('./generate/coevo_H13579_PlosOne_sqrt16_70_290.svg', dpi=600)
plt.show()

# print(f"热图已保存至: ./result/fig/coevo_heatmap_final_perfect.png")
# print(f"位点范围: {start_loc}-{end_loc}，色条长度完全匹配")