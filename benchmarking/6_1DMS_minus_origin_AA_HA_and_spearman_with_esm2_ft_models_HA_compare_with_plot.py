import pandas as pd
import numpy as np
import scipy.stats as st

with open('./viruses-08-00155-s002/Supplemental_File_2_HApreferences.txt', 'r') as f1:
    lines = f1.readlines()
    # print(type(lines))

line_aa = str(lines[0])
line_aa = line_aa.split(' ')

line_aa = line_aa[4:-1]
aa_lst = [i[-1] for i in line_aa]
aa_lst.append('Y')
print(len(aa_lst))
WT_seq = 'MKAKLLVLLYAFVATDADTICIGYHANNSTDTVDTILEKNVAVTHSVNLLEDSHNGKLCKLKGIAPLQLGKCNITGWLLGNPECDSLLPARSWSYIVETPNSENGACYPGDLIDYEELREQLSSVSSLERFEIFPKESSWPNHTFNGVTVSCSHRGKSSFYRNLLWLTKKGDSYPKLTNSYVNNKGKEVLVLWGVHHPSSSDEQQSLYSNGNAYVSVASSNYNRRFTPEIAARPKVRDQHGRMNYYWTLLEPGDTIIFEATGNLIAPWYAFALSRGFESGIITSNASMHECNTKCQTPQGAINSNLPFQNIHPVTIGECPKYVRSTKLRMVTGLRNIPSIQYRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNNLEKRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDLNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCDNECMESVRNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI'
preference_new_lst = []
for i in range(1, len(lines)):
    ref_aa = WT_seq[i]
    print(ref_aa)
    line = str(lines[i])
    line_split = line.split(' ')
    preference_score_lst0 = []
    score_line_lst = line_split[3:]
    if len(score_line_lst) != 20:
        print(len(score_line_lst))
    for score0 in score_line_lst:
        score0 = float(score0)
        preference_score_lst0.append(score0)
    dict_aa_score = dict(zip(aa_lst,preference_score_lst0))
    ref_score = float(dict_aa_score[ref_aa])
    # print(len(dict_aa_score))
    preference_score_lst01 = [float(i)-ref_score for i in preference_score_lst0]
    preference_new_lst+=preference_score_lst01
df00 = pd.read_csv('./DMS_esm2_logist_matrix/HA_DMS_for_esm2_logist_matrix_labeled_esm2_t33_650M_UR50D.csv')
df00['preference_new'] = preference_new_lst
# df.to_csv('./DMS_esm2_logist_matrix/HA_DMS_for_esm2_logist_matrix_labeled_esm2_t33_650M_UR50D_with_new_preference.csv',index=False)

esm2_lst00 = list(df00['esm2_t33_650M_UR50D'])
corr10,p10 = st.spearmanr(preference_new_lst,esm2_lst00)

cp_lst = ['esm2_650M_origin']
corr_lst = [corr10]
p_value_lst = [p10]
print('esm2_650M_origin',corr10,p10)
# checkpoints = [i for i in range(200,15280+400,400)]
checkpoints = ['_sorted']+['_ft'+str(i) for i in range(200,1000+200,200)]+['_ft'+'1180']

for cp in checkpoints:
    cp1 = str(cp)
    df = pd.read_csv('../ESM-2_ft/HA_logist_new_cdhit999_freeze32_lr5e-6/HA_mutation_effects_esm2'+cp1+'.csv')
    esm2_lst = list(df['delta_logit'])

    corr1,p1 = st.spearmanr(preference_new_lst,esm2_lst)
    cp_lst.append(cp)
    corr_lst.append(corr1)
    p_value_lst.append(p1)

    print(cp1,corr1,p1)
df_11 = pd.DataFrame()
df_11['checkpoint'] = cp_lst
df_11['spearman_corr'] = corr_lst
df_11['spearman_p'] = p_value_lst
df_11.to_csv('../ESM-2_ft/HA_logist_new_cdhit999_freeze32_lr5e-6/spearman_test_results_esm2_650M_and_ft_models.csv',index=False)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

# -------------------------- 1. 数据加载与预处理 --------------------------
df = pd.read_csv('../ESM-2_ft/HA_logist_new_cdhit999_freeze32_lr5e-6/spearman_test_results_esm2_650M_and_ft_models.csv')

def extract_step(checkpoint):
    if checkpoint == '_sorted':
        return 0
    elif 'ft' in checkpoint:
        return int(checkpoint.replace('_ft', ''))
    else:
        return 0

df['step'] = df['checkpoint'].apply(extract_step)
baseline_corr = df.loc[0, 'spearman_corr']
df_models = df.iloc[1:].reset_index(drop=True)
y_max = max(baseline_corr, df_models['spearman_corr'].max()) + 0.02

# -------------------------- 2. 绘图样式配置 --------------------------
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

fig, ax = plt.subplots(figsize=(10, 6))
# 关键修复：允许显示坐标轴外的元素
ax.set_clip_on(False)

# -------------------------- 3. 自定义渐变色 --------------------------
# colors = ['#8ecae6', '#219ebc', '#023047', '#ffb703', '#fb8500']
# n_bins = len(df_models)
# cmap = LinearSegmentedColormap.from_list('custom_gradient', colors, N=n_bins)
# bar_colors = cmap(np.linspace(0, 1, n_bins))

from matplotlib.colors import LinearSegmentedColormap
# 马卡龙低饱和蓝绿渐变（浅绿→浅蓝）
# start_color = '#c8e6c9'  # 马卡龙浅绿
# end_color = '#bbdefb'    # 马卡龙浅蓝
start_color = '#a8d5a9'  # 马卡龙浅绿
end_color = '#90caf9'    # 马卡龙浅蓝
n_bins = len(df_models)
cmap = LinearSegmentedColormap.from_list('macaron_blue_green', [start_color, end_color], N=n_bins)
bar_colors = cmap(np.linspace(0, 1, n_bins))

# -------------------------- 4. 绘制基准虚线 --------------------------
ax.axhline(y=baseline_corr, color='#808080', linestyle='--', linewidth=2,
           label=f'Baseline (esm2_650M_origin): {baseline_corr:.3f}')

# -------------------------- 5. 绘制渐变柱状图 --------------------------
bars = ax.bar(
    x=df_models['step'],
    height=df_models['spearman_corr'],
    color=bar_colors,
    width=100,
    edgecolor='white',
    linewidth=1
)

# -------------------------- 6. 柱子顶端添加数值标签 --------------------------
for bar, val in zip(bars, df_models['spearman_corr']):
    height = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        height + 0.003,
        f'{val:.3f}',
        ha='center', va='bottom',
        fontsize=12
    )
# -------------------------- 7. 绘制X轴下方灰色横线+显著性点 --------------------------
x_min, x_max = df_models['step'].min() - 150, df_models['step'].max() + 150
ax.set_xlim(x_min, x_max)

# 调整y轴范围，为下方的显著性标记留出空间
y_bottom = 0.48  # 降低底部边界
ax.set_ylim(bottom=y_bottom, top=y_max)

# 横线纵坐标：设置在y轴底部稍微上方一点的位置
sig_line_y = y_bottom + 0.005
ax.plot([x_min, x_max], [sig_line_y, sig_line_y], color='#808080', linewidth=1.2, zorder=3)

# 绘制显著性点
for idx, row in df_models.iterrows():
    x_pos = row['step']
    p_val = row['spearman_p']
    color = '#d64045' if p_val < 0.05 else '#3d5a80'
    # 将点放在横线稍微上方一点
    ax.scatter(x_pos, sig_line_y, color=color, s=50, zorder=5,
               edgecolors='white', linewidth=0.5)

# 添加标注说明
# ax.text(x_min + (x_max - x_min) * 0.5, sig_line_y + 0.003,
#         'Significance: ● = p<0.05  ● = p≥0.05',
#         ha='center', va='bottom', fontsize=9, color='#666666')

# -------------------------- 8. 图表美化与标签 --------------------------
ax.set_xlabel('Training Step', fontsize=14)
ax.set_ylabel('Spearman Correlation Coefficient', fontsize=14)
# ax.set_title('Spearman Correlation of ESM2 Models with Preference Scores',
#              fontsize=16, pad=25)  # 增加顶部边距

# 修改图例位置，避免与新增内容重叠
custom_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#d64045',
           markersize=8, markeredgecolor='white', markeredgewidth=0.5, label='p < 0.05'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#3d5a80',
           markersize=8, markeredgecolor='white', markeredgewidth=0.5, label='p ≥ 0.05')
]
legend = ax.legend(handles=ax.get_legend_handles_labels()[0] + custom_handles,
                   loc='upper right', fontsize=12, frameon=True,
                   shadow=True, borderpad=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 调整底部边距，确保下方内容完全显示
plt.subplots_adjust(bottom=0.2)  # 增加底部边距

# 确保x轴位置保持在0
# ax.axhline(y=0, color='black', linewidth=1, zorder=1)

# -------------------------- 9. 保存图片 --------------------------
plt.tight_layout()
# plt.savefig('./spearman_correlation_plot_final_fixed_HA_freeze32_lr5e-6.png', dpi=600, bbox_inches='tight')
plt.savefig('./generate/spearman_correlation_plot_final_fixed_HA_freeze32_lr5e-6.svg', dpi=600)
plt.show()