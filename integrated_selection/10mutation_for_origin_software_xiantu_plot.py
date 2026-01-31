import pandas as pd

# 替换为你指定的文件路径
input_file_path = "./final_selected/H1_vBERT_and_ESM2_ft_selected.csv"
output_file_path = "./generate/H1_vBERT_and_ESM2_ft_selected_for_xiantu.csv"

# input_file_path = "./final_selected/H3_vBERT_and_ESM2_ft_selected.csv"
# output_file_path = "./generate/H3_vBERT_and_ESM2_ft_selected_for_xiantu.csv"

# input_file_path = "./final_selected/H5_vBERT_and_ESM2_ft_selected.csv"
# output_file_path = "./generate/H5_vBERT_and_ESM2_ft_selected_for_xiantu.csv"

# input_file_path = "./final_selected/H7_vBERT_and_ESM2_ft_selected.csv"
# output_file_path = "./generate/H7_vBERT_and_ESM2_ft_selected_for_xiantu.csv"

# input_file_path = "./final_selected/H9_vBERT_and_ESM2_ft_selected.csv"
# output_file_path = "./generate/H9_vBERT_and_ESM2_ft_selected_for_xiantu.csv"

# 1. 读取原始CSV文件
# 若读取时出现编码错误，可将 encoding='utf-8' 改为 encoding='gbk'
df = pd.read_csv(input_file_path, encoding='utf-8')

# 2. 数据转换：拼接突变位点格式（I116M、K155N）
# 拆分site1_mut并拼接位点数字
df['mut1'] = df['site1_mut'].str.split('_').str[0] + df['site1_from1_plosone'].astype(str) + df['site1_mut'].str.split('_').str[1]
# 拆分site2_mut并拼接位点数字
df['mut2'] = df['site2_mut'].str.split('_').str[0] + df['site2_from1_plosone'].astype(str) + df['site2_mut'].str.split('_').str[1]

# 3. 添加第三列，所有值固定为1
df['weight'] = 1

# 4. 筛选出需要的三列并保存为新CSV
result_df = df[['mut1', 'mut2', 'weight']]
result_df.to_csv(output_file_path, index=False, encoding='utf-8')

print(f"转换完成！新文件已保存至：{output_file_path}")
print("\n转换后的数据前5行预览：")
print(result_df.head())