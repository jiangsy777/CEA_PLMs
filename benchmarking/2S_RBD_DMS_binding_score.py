import pandas as pd
import numpy as np
import scipy.stats as st

df_csv = pd.read_csv('./dms_SARS2binding/bind_expr_BA1_Nature2022.csv')
score_lst0 = list(df_csv['bind_avg'])

#score_lst = []
score_each_site_sum_lst = []
score_each_site_var_lst = []
score_each_site_CV_lst = []
score_each_site_ptp_lst = []
score_each_site_kurt_lst = []

for i in range(0,len(score_lst0),20):
    score_lst00 = score_lst0[i:(i+20)]
    score_each_site_sum_lst.append(np.sum(score_lst00))  #sum
    score_each_site_var_lst.append(np.var(score_lst00))  # 方差
    score_each_site_CV_lst.append(np.std(score_lst00) / np.mean(score_lst00))  # 变异系数
    score_each_site_ptp_lst.append(np.ptp(score_lst00))  # 极差
    score_each_site_kurt_lst.append(st.kurtosis(score_lst00))  # 峰度

score_each_site_sum_lst1 = [
    (i - min(score_each_site_sum_lst)) / (max(score_each_site_sum_lst) - min(score_each_site_sum_lst))
    for i in score_each_site_sum_lst]
score_each_site_var_lst1 = [
    (i - min(score_each_site_var_lst)) / (max(score_each_site_var_lst) - min(score_each_site_var_lst))
    for i in score_each_site_var_lst]
score_each_site_CV_lst1 = [
    (i - min(score_each_site_CV_lst)) / (max(score_each_site_CV_lst) - min(score_each_site_CV_lst)) for i
    in score_each_site_CV_lst]
score_each_site_ptp_lst1 = [
    (i - min(score_each_site_ptp_lst)) / (max(score_each_site_ptp_lst) - min(score_each_site_ptp_lst))
    for i in score_each_site_ptp_lst]
score_each_site_kurt_lst1 = [
    (i - min(score_each_site_kurt_lst)) / (max(score_each_site_kurt_lst) - min(score_each_site_kurt_lst))
    for i in score_each_site_kurt_lst]

score_each_site_sum = []
for i in range(len(score_each_site_var_lst1)):
    score_each_site_sum.append(score_each_site_var_lst1[i]+score_each_site_CV_lst1[i]+score_each_site_ptp_lst1[i]+score_each_site_kurt_lst1[i])
print(len(score_each_site_sum))
df_ = pd.DataFrame()
df_['site_from1'] = [k for k in range(331,532)]
df_['score_each_site_var'] = score_each_site_var_lst1
df_['score_each_site_CV'] = score_each_site_CV_lst1
df_['score_each_site_ptp'] = score_each_site_ptp_lst1
df_['score_each_site_kurt'] = score_each_site_kurt_lst1
df_['score_each_site_sum'] = score_each_site_sum
df_['score_each_site_sum_pos_nag'] = score_each_site_sum_lst1
df_.to_csv('./generate/binding_avg_value_each_site_S_RBD_DMS.csv',index=False)
