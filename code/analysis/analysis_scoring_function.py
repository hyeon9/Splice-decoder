#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import sys
import argparse
import os
import HGfig
from scipy.stats import pearsonr, mannwhitneyu, spearmanr
# %%
data = pd.read_csv("/home/kangh/lab-server/Tool/Splice-decoder_git/example/MYC8h/result/Whole_DS_score_Whole.txt",
                   sep="\t")
data["Delta_PSI"] = np.abs(data["Delta_PSI"])
print(data.columns)
# %%
## Check unknown changes
sub = (data[["Domain_change_rate","Delta_Amino_acid"]])
print(sub[(sub["Domain_change_rate"]==0) & (sub["Delta_Amino_acid"]!=1.0)])
# %%
sns.set_style("darkgrid")
for i in ["Domain_change_rate","Transcript_usage","Delta_PSI","Max_PSI"]:
    plt.figure(figsize=(4,4))
    plt.scatter(x="Effect_Score",
                y=i,
                data=data)
    r, p = pearsonr(data[i],data["Effect_Score"])
    plt.ylabel(i)
    plt.xlabel("Effect_Score")
    plt.text(0.4,0.2,f"r= {round(r,3)}",
             fontsize=15)
# %%
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

# %%
# 독립 변수와 종속 변수 설정
X = data[["Domain_change_rate","Transcript_usage","Delta_PSI","Max_PSI"]]
X = sm.add_constant(X)  # 절편 추가
print(X)
y = data['Effect_Score']

# 회귀 모델 적합
model = sm.OLS(y, X).fit()
print(model.summary())  # 회귀 요약 결과

# VIF (다중공선성 확인)
vif_df = pd.DataFrame()
vif_df['variable'] = X.columns
vif_df['VIF'] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
print(vif_df)
# %%
