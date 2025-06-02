# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import pearsonr, mannwhitneyu, ks_2samp
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({
'axes.titlesize': 15,     # 제목 글꼴 크기
'axes.labelsize': 14,     # x, y축 라벨 글꼴 크기
'xtick.labelsize': 12,    # x축 틱 라벨 글꼴 크기
'ytick.labelsize': 12,    # y축 틱 라벨 글꼴 크기
'legend.fontsize': 12,    # 범례 글꼴 크기
'figure.titlesize': 15    # figure 제목 글꼴 크기
})

data_type = "MYC8"
if data_type == "MYC8":
    DIR = "/home/kangh/lab-server/Tool/Splice-decoder_git/example/MYC8h/result/"
## From SpliceDecoder
data = pd.read_csv(DIR+"Whole_DS_score_Whole_w_meanPs.txt",
                   sep="\t")
data2 = pd.read_csv(DIR+"Whole_DS_score_Whole.txt",
                   sep="\t")
data = data[data["Domain_change_rate"]>0]  # Collect only functional events
data2 = data2[data2["Domain_change_rate"]>0]  # Collect only functional events
# %%
merged = pd.merge(data,
                  data2,
                  left_on="LongID",
                  right_on="LongID")
# %%
plt.figure(figsize=(4,4))
plt.scatter(x="Effect_Score_x",
            y="Effect_Score_y",
            data=merged)
# %%
