import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_1 = pd.read_csv('./generate/preference_value_each_site_H1_HA_DMS.csv')
df_2 = pd.read_csv('./generate/site_coevo_mean_H1_HA_DMS.csv')
df_3 = pd.read_csv('./DCA/site_dca_mean_H1_HA_DMS.csv')
df_4 = pd.read_csv('./DMS_esm2_logist_matrix/esm2_logist_each_site_HA_DMS.csv')

site_coevo_lst = list(df_2['site_coevo_mean'])
site_dca_lst = list(df_3['site_DCA_mean'])
site_esm2_lst = list(df_4['score_each_site_sum'])

site_coevo_lst = site_coevo_lst[1:]#remove site1 score
site_DCA_lst = site_dca_lst[1:]#remove site1 score

df_1['site_coevo_mean'] = site_coevo_lst
df_1['site_DCA_mean'] = site_DCA_lst
df_1['site_esm2'] = site_esm2_lst

site_coevo_lst_normalized = [(i-min(site_coevo_lst))/(max(site_coevo_lst)-min(site_coevo_lst)) for i in site_coevo_lst]
site_DCA_lst_normalized = [(i-min(site_DCA_lst))/(max(site_DCA_lst)-min(site_DCA_lst)) for i in site_DCA_lst]
site_esm2_lst_normalized = [(i-min(site_esm2_lst))/(max(site_esm2_lst)-min(site_esm2_lst)) for i in site_esm2_lst]

df_1['site_coevo_mean_normalized'] = site_coevo_lst_normalized
df_1['site_DCA_mean_normalized'] = site_DCA_lst_normalized
df_1['site_esm2_normalized'] = site_esm2_lst_normalized

# df_1 = df_1[df_1['site_coevo_mean']!=0]#remove sites with low mutation rate
df_1 = df_1.sort_values(by='preference_each_site_sum',ascending=False)
df_1.to_csv('site_preference_coevo_DCA_esm2_H1_HA_DMS_without_sites_with_low_mutation_rate.csv',index=False)

df_csv = pd.read_csv('site_preference_coevo_DCA_esm2_H1_HA_DMS_without_sites_with_low_mutation_rate.csv')
site_coevo_lst1 = list(df_csv['site_coevo_mean_normalized'])
site_DCA_lst1 = list(df_csv['site_DCA_mean_normalized'])
site_esm2_lst1 = list(df_csv['site_esm2_normalized'])
preference_each_site_sum_lst1 = list(df_csv['preference_each_site_sum'])

site_coevo_high_lst = []
site_coevo_low_lst = []
site_DCA_high_lst = []
site_DCA_low_lst = []
site_esm2_high_lst = []
site_esm2_low_lst = []
preference_high_lst = []
preference_low_lst = []

num_lst = []

for i in range(10,250,10):
    num_lst.append(i)
    site_coevo_high_lst.append(np.mean(site_coevo_lst1[:i]))
    site_coevo_low_lst.append(np.mean(site_coevo_lst1[-i:]))

    site_DCA_high_lst.append(np.mean(site_DCA_lst1[:i]))
    site_DCA_low_lst.append(np.mean(site_DCA_lst1[-i:]))

    site_esm2_high_lst.append(np.mean(site_esm2_lst1[:i]))
    site_esm2_low_lst.append(np.mean(site_esm2_lst1[-i:]))

    preference_high_lst.append(np.mean(preference_each_site_sum_lst1[:i]))
    preference_low_lst.append(np.mean(preference_each_site_sum_lst1[-i:]))
    # print(len(site_coevo_lst1[:i]))
    # print(len(site_coevo_lst1[-i:]))
df_final = pd.DataFrame()
df_final['comparing_site_num'] = num_lst
df_final['site_coevo_high'] = site_coevo_high_lst
df_final['site_coevo_low'] = site_coevo_low_lst
df_final['site_DCA_high'] = site_DCA_high_lst
df_final['site_DCA_low'] = site_DCA_low_lst
df_final['site_esm2_high'] = site_esm2_high_lst
df_final['site_esm2_low'] = site_esm2_low_lst
df_final['preference_high'] = preference_high_lst
df_final['preference_low'] = preference_low_lst
df_final.to_csv('./generate/DMS_H1_HA_benchmarking_coevo_vs_DCA_esm2_result.csv',index=False)

# from scipy import stats
# s1,p1 = stats.kruskal(preference_high_lst,preference_low_lst)
# s2,p2 = stats.kruskal(site_coevo_high_lst,site_coevo_low_lst)
# s3,p3 = stats.kruskal(site_DCA_high_lst,site_DCA_low_lst)
# s4,p4 = stats.kruskal(site_esm2_high_lst,site_esm2_low_lst)
#
# print(p1,p2,p3,p4)


def statistical_test(before,after,name):
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt

    # 可视化
    plt.boxplot([before, after], labels=['Before', 'After'])
    plt.show()

    # 正态性检验
    differences = after - before
    shapiro_stat, shapiro_p = stats.shapiro(differences)
    # print(f"Shapiro-Wilk检验p值: {shapiro_p:.4f}")

    if shapiro_p > 0.05:
        # 配对t检验
        test_method = 'ttest_rel'
        t_stat, p_value = stats.ttest_rel(after, before,alternative='less')
        # print(f"t统计量: {t_stat:.4f}, p值: {p_value:.4f}")
        # if p_value < 0.05:
        #     print("差异显著")
        # else:
        #     print("差异不显著")

        # 效应量
        cohen_d = np.mean(differences) / np.std(differences, ddof=1)
        #print(f"Cohen's d: {cohen_d:.4f}")
        print(name + '_calculation_end!!!!!!!!!!!!!!!!!!!!!!!')
        return (test_method,p_value)
    else:
        # Wilcoxon检验
        test_method = 'wilcoxon'
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(after, before,alternative='less')
        # print(f"Wilcoxon检验p值: {wilcoxon_p:.4f}")
        print(name+'_calculation_end!!!!!!!!!!!!!!!!!!!!!!!')
        return (test_method,wilcoxon_p)





import numpy as np
from scipy import stats
###################################DMS################################
# 示例数据（同一组对象的前后测值）
before = np.array(preference_high_lst)
after = np.array(preference_low_lst)
p1 = statistical_test(before,after,'preference')
print(p1)

###################################coevo################################
before = np.array(site_coevo_high_lst)
after = np.array(site_coevo_low_lst)
p1 = statistical_test(before,after,'coevo')
print(p1)

###################################DCA################################
before = np.array(site_DCA_high_lst)
after = np.array(site_DCA_low_lst)
p1 = statistical_test(before,after,'dca')
print(p1)

###################################esm2################################
before = np.array(site_esm2_high_lst)
after = np.array(site_esm2_low_lst)
p1 = statistical_test(before,after,'esm2')
print(p1)


# -------------------------- 组内独立0-1标准化 --------------------------
def min_max_normalize(high_lst, low_lst):
    combined = np.array(high_lst + low_lst)
    min_val = combined.min()
    max_val = combined.max()
    if max_val == min_val:
        return np.zeros_like(high_lst), np.zeros_like(low_lst)
    high_norm = (high_lst - min_val) / (max_val - min_val)
    low_norm = (low_lst - min_val) / (max_val - min_val)
    return high_norm, low_norm


# 每组特征单独标准化
preference_high_norm, preference_low_norm = min_max_normalize(preference_high_lst, preference_low_lst)
coevo_high_norm, coevo_low_norm = min_max_normalize(site_coevo_high_lst, site_coevo_low_lst)
dca_high_norm, dca_low_norm = min_max_normalize(site_DCA_high_lst, site_DCA_low_lst)
esm2_high_norm, esm2_low_norm = min_max_normalize(site_esm2_high_lst, site_esm2_low_lst)

# 整理标准化后数据
high_data = [preference_high_norm, coevo_high_norm, dca_high_norm, esm2_high_norm]
low_data = [preference_low_norm, coevo_low_norm, dca_low_norm, esm2_low_norm]
features = ['DMS', 'Coevo', 'DCA', 'ESM2']


# -------------------------- 原统计检验函数（完全未修改） --------------------------
def statistical_test(before, after, name):
    differences = [i - j for i in after for j in before]
    shapiro_stat, shapiro_p = stats.shapiro(differences)
    if shapiro_p > 0.05:
        test_method = 'ttest_rel'
        t_stat, p_value = stats.ttest_rel(after, before, alternative='less')
        cohen_d = np.mean(differences) / np.std(differences, ddof=1)
        print(name + '_calculation_end!!!!!!!!!!!!!!!!!!!!!!!')
        return (test_method, p_value)
    else:
        test_method = 'wilcoxon'
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(after, before, alternative='less')
        print(name + '_calculation_end!!!!!!!!!!!!!!!!!!!!!!!')
        return (test_method, wilcoxon_p)


# -------------------------- 原统计检验执行代码（存储p值） --------------------------
p_values_list = []
# Preference
p1 = statistical_test(preference_high_lst, preference_low_lst, 'preference')
p_values_list.append(p1[1])
print(p1)
# Coevo
p1 = statistical_test(site_coevo_high_lst, site_coevo_low_lst, 'coevo')
p_values_list.append(p1[1])
print(p1)
# DCA
p1 = statistical_test(site_DCA_high_lst, site_DCA_low_lst, 'dca')
p_values_list.append(p1[1])
print(p1)
# ESM2
p1 = statistical_test(site_esm2_high_lst, site_esm2_low_lst, 'esm2')
p_values_list.append(p1[1])
print(p1)

# -------------------------- 定制化小提琴图（核心：配色修改） --------------------------
# 全局字体设置：Arial 不加粗
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.weight'] = 'normal'  # 强制不加粗

# 四组配色方案（严格按要求：深橙/浅橙、深红/浅红、深紫/浅紫、深青/浅青）
# 格式：[深颜色(边框), 浅颜色(填充)]
color_schemes = [
    ['#d95319', '#feb236'],  # Preference: 深橙+浅橙
    ['#c02942', '#f1a9a0'],  # Coevo: 深红+浅红
    ['#7e3a7e', '#cbc0d3'],  # DCA: 深紫+浅紫
    ['#238b45', '#bae4bc']  # ESM2: 深青+浅青
]

fig, ax = plt.subplots(figsize=(14, 9))
x_pos = np.arange(len(features))
width = 0.25  # 小提琴宽度
inner_width = 0.5  # 组内两个小提琴的间距系数

# 绘制小提琴图 + 同色散点
for i in range(len(features)):
    # 增大组内两个小提琴的间距
    pos_high = x_pos[i] - inner_width * width
    pos_low = x_pos[i] + inner_width * width

    # 提取当前组的深浅颜色
    dark_color = color_schemes[i][0]
    light_color = color_schemes[i][1]

    # 绘制小提琴图
    violin_parts = ax.violinplot([high_data[i], low_data[i]],
                                 positions=[pos_high, pos_low],
                                 widths=width,
                                 showmedians=True,
                                 showextrema=False)
    # 设置小提琴样式
    # High组小提琴：浅颜色填充 + 深颜色边框
    violin_parts['bodies'][0].set_facecolor(light_color)
    violin_parts['bodies'][0].set_edgecolor(dark_color)
    violin_parts['bodies'][0].set_alpha(0.7)
    # Low组小提琴：无填充 + 深颜色边框
    violin_parts['bodies'][1].set_facecolor('none')
    violin_parts['bodies'][1].set_edgecolor(dark_color)
    violin_parts['bodies'][1].set_alpha(0.7)
    # 中位数线：黑色加粗
    violin_parts['cmedians'].set_color('#000000')
    violin_parts['cmedians'].set_linewidth(2)

    # 绘制散点（在小提琴上方）
    # High组散点：深颜色填充 + 深颜色轮廓
    ax.scatter(np.random.normal(pos_high, 0.015, len(high_data[i])),
               high_data[i] + 0.03,
               color=dark_color,
               s=30,
               edgecolor=dark_color,
               linewidth=0.5,
               alpha=0.8,
               label='High (Filled)' if i == 0 else "")
    # Low组散点：无填充 + 深颜色轮廓
    ax.scatter(np.random.normal(pos_low, 0.015, len(low_data[i])),
               low_data[i] + 0.03,
               facecolor='none',
               edgecolor=dark_color,
               s=30,
               linewidth=1.2,
               alpha=0.8,
               label='Low (Unfilled)' if i == 0 else "")

    # 标注p值（Arial 不加粗）
    max_val = max(np.max(high_data[i]), np.max(low_data[i]))
    ax.text(x_pos[i], max_val + 0.05, f'p = {p_values_list[i]:.4f}',
            ha='center', va='bottom', fontsize=18, fontweight='normal')

# 绘制组间灰色虚线：精准位于横坐标刻度中点
for i in range(1, len(features)):
    mid_pos = x_pos[i - 1] + 0.5
    ax.axvline(x=mid_pos, color='#999999', linestyle='--', linewidth=1.5)

# 坐标轴设置
ax.set_xticks(x_pos)
ax.set_xticklabels(features, fontsize=20, fontweight='normal')
ax.set_ylabel('Normalized Value (0-1)', fontsize=22, fontweight='normal')
ax.tick_params(axis='y', labelsize=18, width=1.5)
ax.set_ylim(-0.15, 1.3)

# 移除顶部和右侧边框
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# -------------------------- 图例修改：样式+位置 --------------------------
# 获取图例对象并修改样式
# legend = ax.legend(fontsize=18, loc='upper right', bbox_to_anchor=(1, 0.99), frameon=False)
# # 修改图例标记样式：深灰色边框，High填充浅灰，Low无填充
# for handle, label_text in zip(legend.legendHandles, legend.get_texts()):
#     handle.set_edgecolor('#333333')  # 深灰色边框
#     handle.set_linewidth(1.5)
#     if 'High' in label_text.get_text():
#         handle.set_facecolor('#dddddd')  # 浅灰色填充
#     else:
#         handle.set_facecolor('none')  # 无填充
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 自定义图例方块：大小、颜色、边框
high_patch = mpatches.Rectangle((0, 0), 1, 1,
                                facecolor='#bbbbbb',  # 加深的浅灰色填充
                                edgecolor='#333333',  # 深灰色边框
                                linewidth=1.5)        # 边框粗细
low_patch = mpatches.Rectangle((0, 0), 1, 1,
                               facecolor='none',     # 无填充
                               edgecolor='#333333',
                               linewidth=1.5)

# 添加图例到图表，位置可通过bbox_to_anchor调整高低
ax.legend([high_patch, low_patch],
          ['High (Filled)', 'Low (Unfilled)'],
          fontsize=18,
          loc='upper right',
          bbox_to_anchor=(1, 1.0),  # y=1.0 贴紧顶部，可微调
          frameon=False,
          handlelength=1.3,         # 控制方块大小，值越大方块越长
          handleheight=0.7)         # 控制方块高度，和length一致就是正方形
plt.tight_layout()
# plt.savefig('high_low_violin_scatter_final_color_legend_HA.png', dpi=600, bbox_inches='tight')
# plt.savefig('high_low_violin_scatter_final_color_legend_HA.svg',dpi=600)
plt.savefig('./generate/high_low_violin_scatter_final_color_legend_HA.pdf',dpi=600)

plt.show()