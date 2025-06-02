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
from scipy.stats import pearsonr, mannwhitneyu

data = pd.read_csv("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/250530_Whole_DS_score_Whole.txt",
                   sep="\t")
cano = pd.read_csv("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/gencode_info/canonical_gencode_length.txt",
                   sep="\t")
non_cano = pd.read_csv("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/gencode_info/non-cano_gencode_length.txt",
                       sep="\t")

# %%
## Make cano/non-cano DS_score_matrix
cano_df = data[data["Reference_transcript"].isin(cano["0"])]
non_cano_df = data[data["Reference_transcript"].isin(non_cano["0"])]
cano_df = cano_df[["LongID","Effect_Score","Domain_change_rate","Transcript_usage"]]
cano_df["type"] = "Cano"
non_cano_df = non_cano_df[["LongID","Effect_Score","Domain_change_rate","Transcript_usage"]]
non_cano_df["type"] = "Non-cano"

# %%
## Merge them
df = pd.concat([cano_df, non_cano_df])
uncorr = []
corr = []
k = 0
j = 0
## Check the category of representative pair of Effect_Score or delta E
for i in df["LongID"].unique():
    sub = df[df["LongID"]==i]
    if len(sub["type"].unique()) > 1:
        j += 1  # number of DS with canonical and non-canonical pairs
        if (sub[sub["Domain_change_rate"]==sub["Domain_change_rate"].max()]["type"].unique()[0]) == "Non-cano" and\
            (sub[sub["Effect_Score"]==sub["Effect_Score"].max()]["type"].unique()[0]) == "Cano": 
                ## Corr: Delta E pointed non-canonical as representative while Effect_score pointed cacnonical as representative
            corr.append(i)
        elif (sub[sub["Domain_change_rate"]==sub["Domain_change_rate"].max()]["type"].unique()[0]) == "Non-cano" and\
            (sub[sub["Effect_Score"]==sub["Effect_Score"].max()]["type"].unique()[0]) == "Non-cano":
                ## Uncorr: Delta E and Effect_score pointed non-canonical as representative
            uncorr.append(i)
        elif (sub[sub["Domain_change_rate"]==sub["Domain_change_rate"].max()]["type"].unique()[0]) == "Cano":
                ## Delta E pointed canonical as representative
            k += 1

print(len(corr),len(uncorr),k,j)

# %%
corr_df = df[df["LongID"].isin(uncorr+corr)]
print(corr_df)

# %%
figure,axs = plt.subplots(1,2,figsize=(5,3))
plt.subplots_adjust(wspace=1)
sns.boxplot(x="type",y="Domain_change_rate",
            data=corr_df,
            showfliers=False,
            width=.6,
            ax=axs[0])
sns.boxplot(x="type",y="Effect_Score",
            data=corr_df,
            showfliers=False,
            width=.6,
            ax=axs[1])
sns.despine()
# plt.savefig('/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/figure/effect_score_vs_deltaL_in_can_non_cano.pdf',
#             bbox_inches="tight")

print(mannwhitneyu(corr_df[corr_df["type"]=="Cano"]["Domain_change_rate"],
                   corr_df[corr_df["type"]!="Cano"]["Domain_change_rate"],
                   alternative="less"))
print(mannwhitneyu(corr_df[corr_df["type"]=="Cano"]["Effect_Score"],
                   corr_df[corr_df["type"]!="Cano"]["Effect_Score"],
                   alternative="greater"))

# %%
# plt.figure(figsize=(3,3))
# sns.kdeplot(df[df["type"]=="Cano"]["Effect_Score"],label="Cano",)
# sns.kdeplot(df[df["type"]!="Cano"]["Effect_Score"],label="Non-cano",)
# # plt.xlim(0,0.25)
# # sns.stripplot(x="type",y="Effect_Score",
# #               data=df)
# plt.figure(figsize=(3,3))
# plt.figure(figsize=(3,3))
# sns.kdeplot(df[df["type"]=="Cano"]["Domain_change_rate"],label="Cano",)
# sns.kdeplot(df[df["type"]!="Cano"]["Domain_change_rate"],label="Non-cano",)

# # sns.stripplot(x="type",y="Domain_change_rate",
# #               data=df, alpha=.7)
# print(mannwhitneyu(df[df["type"]=="Cano"]["Effect_Score"],
#                    df[df["type"]!="Cano"]["Effect_Score"],
#                    alternative="greater"))
# print(mannwhitneyu(df[df["type"]=="Cano"]["Domain_change_rate"],
#                    df[df["type"]!="Cano"]["Domain_change_rate"],
#                    alternative="greater"))
# %%
mean_df = []
mean_df.append([df[df["type"]=="Cano"]["Domain_change_rate"].mean(), df[df["type"]!="Cano"]["Domain_change_rate"].mean()])
mean_df.append([df[df["type"]=="Cano"]["Effect_Score"].mean(), df[df["type"]!="Cano"]["Effect_Score"].mean()])
mean_df = pd.DataFrame(mean_df).T
mean_df.columns = ["Domain_change_rate","Effect_Score"]
mean_df.index = ["Cano","Non-cano"]
print(df)

# %%
gene1 = data[(data["LongID"].str.contains("MUTYH")) &\
             (data["Simulated_event"].str.contains("RI"))]
gene2 = data[(data["LongID"].str.contains("BCLAF1")) &\
             (data["Simulated_event"].str.contains("ES"))]
print(gene2)
# %%
fig, ax = plt.subplots(figsize = (4,3))
plt.title('BCLAF1',
          fontstyle="italic")
plt.xticks(rotation=90)
ax2 = ax.twinx() 
# using the twinx() for creating another
# axes object for secondary y-Axis
gene2 = gene2.sort_values(by="Effect_Score",
                          ascending=False)
sns.scatterplot(x="Reference_transcript",
                y="Effect_Score",
                data=gene2,
                ax=ax,
                color="#2B75FFFF",
                s=100)
sns.scatterplot(x="Reference_transcript",
                y="Domain_change_rate",
                data=gene2,
                ax=ax2,
                color="#CECECEFF",
                s=60)

ax.set_ylabel("Effect Score")
ax2.set_ylabel("Delta L")
ax.set_xlabel(None)
sns.despine(right=False)
# plt.savefig(f"/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/figure/BCLAF1_effect.pdf",
#             bbox_inches="tight")
# %%
