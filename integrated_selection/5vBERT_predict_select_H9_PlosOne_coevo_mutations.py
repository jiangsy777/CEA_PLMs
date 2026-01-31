import pandas as pd
import warnings
warnings.filterwarnings('ignore')  # 屏蔽空值均值警告

# ===================== 配置参数 =====================
CSV1_PATH = "../vBERT_classifier_establishment/Predicted_prob_H13579_PlosOne_aligned_cdhit98_test_seq_split_random30.csv"
CSV2_PATH = "./generate/HA_H9_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score.csv"
SAVE_PATH = "./generate/HA_H9_avian_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score_with_vBERT_predicted_prob_mean.csv"

# ===================== 读取数据 =====================
df1 = pd.read_csv(CSV1_PATH)
df2 = pd.read_csv(CSV2_PATH)
df2 = df2.drop(index=0).reset_index(drop=True)

# 关键：位点编号从 1 转 0 起始（向量化操作提升效率）
df2["site1_from1_plosone"] = df2["site1_from1_plosone"].astype(int) - 1
df2["site2_from1_plosone"] = df2["site2_from1_plosone"].astype(int) - 1

# ===================== 提取突变前后的氨基酸 =====================
def get_aa(mut_str):
    """从mut_str（如 A_F）提取突变前后氨基酸：前为野生型，后为突变型"""
    parts = mut_str.split("_")
    return parts[0], parts[-1]

# 向量化提取，替代zip(*apply)避免tuple解包警告
df2[["site1_aa_wt", "site1_aa_mut"]] = df2["site1_mut"].apply(get_aa).tolist()
df2[["site2_aa_wt", "site2_aa_mut"]] = df2["site2_mut"].apply(get_aa).tolist()

# ===================== 统计突变前后的概率均值 + 计算差值 =====================
# 初始化结果列表，全部填充为 "single"
total_rows = len(df2)
avian_wt_mean_list = ["single"] * total_rows
human_wt_mean_list = ["single"] * total_rows
avian_mut_mean_list = ["single"] * total_rows
human_mut_mean_list = ["single"] * total_rows

# 缓存序列列，避免循环内重复索引
seq_series = df1["seq"]
avian_score = df1["risk_score_avian"]
human_score = df1["risk_score_human"]

# 从第3行（索引2）开始，每隔2行计算一次（步长3：索引2,5,8...）
for idx in range(2, total_rows, 3):
    row = df2.iloc[idx]
    pos1, pos2 = row["site1_from1_plosone"], row["site2_from1_plosone"]
    aa1_wt, aa2_wt = row["site1_aa_wt"], row["site2_aa_wt"]
    aa1_mut, aa2_mut = row["site1_aa_mut"], row["site2_aa_mut"]

    # --- 突变前均值计算 ---
    mask_wt = (seq_series.str[pos1] == aa1_wt) & (seq_series.str[pos2] == aa2_wt)
    avian_wt_mean = avian_score[mask_wt].mean() if not mask_wt.sum() == 0 else float("nan")
    human_wt_mean = human_score[mask_wt].mean() if not mask_wt.sum() == 0 else float("nan")

    # --- 突变后均值计算 ---
    mask_mut = (seq_series.str[pos1] == aa1_mut) & (seq_series.str[pos2] == aa2_mut)
    avian_mut_mean = avian_score[mask_mut].mean() if not mask_mut.sum() == 0 else float("nan")
    human_mut_mean = human_score[mask_mut].mean() if not mask_mut.sum() == 0 else float("nan")

    # 仅更新当前索引行的结果
    avian_wt_mean_list[idx] = avian_wt_mean
    human_wt_mean_list[idx] = human_wt_mean
    avian_mut_mean_list[idx] = avian_mut_mean
    human_mut_mean_list[idx] = human_mut_mean

# ===================== 赋值并计算差值 =====================
# 突变前后均值列
df2["Prob_Avian_Before_Mutation"] = avian_wt_mean_list
df2["Prob_Human_Before_Mutation"] = human_wt_mean_list
df2["Prob_Avian_After_Mutation"] = avian_mut_mean_list
df2["Prob_Human_After_Mutation"] = human_mut_mean_list

# 位点编号还原为 1 起始（向量化操作）
df2["site1_from1_plosone"] = df2["site1_from1_plosone"].astype(int) + 1
df2["site2_from1_plosone"] = df2["site2_from1_plosone"].astype(int) + 1

# 差值列：后 - 前，仅对数值行计算，非数值行填充single
df2["delta_avian_prob"] = df2.apply(
    lambda x: x["Prob_Avian_After_Mutation"] - x["Prob_Avian_Before_Mutation"]
    if isinstance(x["Prob_Avian_After_Mutation"], (int, float)) else "single", axis=1
)
df2["delta_human_prob"] = df2.apply(
    lambda x: x["Prob_Human_After_Mutation"] - x["Prob_Human_Before_Mutation"]
    if isinstance(x["Prob_Human_After_Mutation"], (int, float)) else "single", axis=1
)

# ===================== 保存结果 =====================
df2.to_csv(SAVE_PATH, index=False, encoding="utf-8-sig")
print(f"结果已保存至 {SAVE_PATH}")
print(f"共处理 {len(df2)} 条位点突变数据")
# 输出计算行的索引，方便验证
calc_indices = list(range(2, total_rows, 3))
print(f"计算的行索引为: {calc_indices[:10]}...")  # 仅打印前10个