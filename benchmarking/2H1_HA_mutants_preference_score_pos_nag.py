import pandas as pd
import numpy as np
import scipy.stats as st

with open('./viruses-08-00155-s002/Supplemental_File_2_HApreferences.txt', 'r') as f1:
    lines = f1.readlines()
    # print(type(lines))
preference_score_lst = []
preference_each_site_sum_lst = []
preference_each_site_var_lst = []
preference_each_site_CV_lst = []
preference_each_site_ptp_lst = []
preference_each_site_kurt_lst = []

line_aa = str(lines[0])
line_aa = line_aa.split(' ')

line_aa = line_aa[4:-1]
aa_lst = [i[-1] for i in line_aa]
aa_lst.append('Y')
print(len(aa_lst))
WT_seq = 'MKAKLLVLLYAFVATDADTICIGYHANNSTDTVDTILEKNVAVTHSVNLLEDSHNGKLCKLKGIAPLQLGKCNITGWLLGNPECDSLLPARSWSYIVETPNSENGACYPGDLIDYEELREQLSSVSSLERFEIFPKESSWPNHTFNGVTVSCSHRGKSSFYRNLLWLTKKGDSYPKLTNSYVNNKGKEVLVLWGVHHPSSSDEQQSLYSNGNAYVSVASSNYNRRFTPEIAARPKVRDQHGRMNYYWTLLEPGDTIIFEATGNLIAPWYAFALSRGFESGIITSNASMHECNTKCQTPQGAINSNLPFQNIHPVTIGECPKYVRSTKLRMVTGLRNIPSIQYRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNNLEKRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDLNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCDNECMESVRNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI'

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
        preference_score_lst.append(score0)
        preference_score_lst0.append(score0)
    dict_aa_score = dict(zip(aa_lst,preference_score_lst0))
    ref_score = float(dict_aa_score[ref_aa])
    # print(len(dict_aa_score))
    preference_score_lst01 = [float(i)-ref_score for i in preference_score_lst0]

    preference_each_site_sum_lst.append(np.sum(preference_score_lst01))  # sum
    preference_each_site_var_lst.append(np.var(preference_score_lst01))  # 方差
    preference_each_site_CV_lst.append(np.std(preference_score_lst01) / np.mean(preference_score_lst01))  # 变异系数
    preference_each_site_ptp_lst.append(np.ptp(preference_score_lst01))  # 极差
    preference_each_site_kurt_lst.append(st.kurtosis(preference_score_lst01))  # 峰度


preference_each_site_sum_lst1 = [
    (i - min(preference_each_site_sum_lst)) / (max(preference_each_site_sum_lst) - min(preference_each_site_sum_lst))
    for i in preference_each_site_sum_lst]
preference_each_site_var_lst1 = [
    (i - min(preference_each_site_var_lst)) / (max(preference_each_site_var_lst) - min(preference_each_site_var_lst))
    for i in preference_each_site_var_lst]
preference_each_site_CV_lst1 = [
    (i - min(preference_each_site_CV_lst)) / (max(preference_each_site_CV_lst) - min(preference_each_site_CV_lst)) for i
    in preference_each_site_CV_lst]
preference_each_site_ptp_lst1 = [
    (i - min(preference_each_site_ptp_lst)) / (max(preference_each_site_ptp_lst) - min(preference_each_site_ptp_lst))
    for i in preference_each_site_ptp_lst]
preference_each_site_kurt_lst1 = [
    (i - min(preference_each_site_kurt_lst)) / (max(preference_each_site_kurt_lst) - min(preference_each_site_kurt_lst))
    for i in preference_each_site_kurt_lst]



preference_each_site_sum = []
for i in range(len(preference_each_site_var_lst1)):
    preference_each_site_sum.append(preference_each_site_var_lst1[i]+preference_each_site_CV_lst1[i]+preference_each_site_ptp_lst1[i]+preference_each_site_kurt_lst1[i])
print(len(preference_each_site_sum))
df_ = pd.DataFrame()
df_['site_from1'] = [k for k in range(2,566)]
df_['preference_each_site_var'] = preference_each_site_var_lst1
df_['preference_each_site_CV'] = preference_each_site_CV_lst1
df_['preference_each_site_ptp'] = preference_each_site_ptp_lst1
df_['preference_each_site_kurt'] = preference_each_site_kurt_lst1
df_['preference_each_site_sum'] = preference_each_site_sum
df_['preference_each_site_sum_pos_nag'] = preference_each_site_sum_lst1
df_.to_csv('./generate/preference_value_each_site_H1_HA_DMS.csv',index=False)
# if np.sum(preference_score_lst0) != 1:
#     print(np.sum(preference_score_lst0))
#     # print(preference_score_lst0)
#     print('not1!!!')

# mean_score = np.mean(preference_score_lst)
# median_score = np.median(preference_score_lst)
# label_mean_lst = []
# label_median_lst = []
# for score_1 in preference_score_lst:
#     if score_1 < mean_score:
#         label_mean_lst.append(0)
#     else:
#         label_mean_lst.append(1)
# for score_2 in preference_score_lst:
#     if score_2 < median_score:
#         label_median_lst.append(0)
#     else:
#         label_median_lst.append(1)
# df_preference_score = pd.DataFrame()
# df_preference_score['preference_score'] = preference_score_lst
# df_preference_score['label_mean'] = label_mean_lst
# df_preference_score['label_median'] = label_median_lst
# df_preference_score.to_csv('preference_score_H1_HA_mutants.csv',index=False)
