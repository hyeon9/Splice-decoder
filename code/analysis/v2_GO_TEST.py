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
    DIR = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/"
else:
    DIR = "/Users/kangh/Documents/Analysis_note/isoform_translation/MYC_24h/result/"
OUT = DIR+"figure/"

## From SpliceDecoder
data = pd.read_csv(DIR+"250530_Whole_DS_score_Whole.txt",
                   sep="\t")
data = data[data["Domain_change_rate"]>0]  # Collect only functional events
# %%
# %%
## From rMATS
data2 = pd.read_csv("/".join(DIR.split("/")[:-2])+"/rmat.csv",
                    sep=",")
data2 = data2[(data2["FDR"] < 0.05) &
              (np.abs(data2["dPSI_2_minus_1"]) > 0.1)]

data2 = data2[["geneSymbol"]].astype(str)
data2.columns = ["gene"]
print(f'SpliceDecoder Genes: {len(data["gene"].unique())}')
print(f'rMATS Genes: {len(data2["gene"].unique())}')

# %%
## Effect scores from all DSEs
data3 = pd.read_csv(DIR+"250530_Whole_DS_score_Whole.txt",
                    sep="\t")
gene_score_whole = {}
for gene in data3["gene"].unique():
    gene_score_whole[gene] = []
    temp = data3[data3["gene"]==gene]["Effect_Score"].max() # Select Maximum effect score for each gene
    gene_score_whole[gene].append(temp)

gene_score_whole = pd.DataFrame.from_dict(gene_score_whole, orient="index")

# %%
## Effect scores from NMD and DOA
gene_score = {}
for gene in data["gene"].unique():
    gene_score[gene] = []
    temp = data[(data["gene"]==gene)]["Effect_Score"].max()
    gene_score[gene].append(temp)

gene_score = pd.DataFrame.from_dict(gene_score, orient="index")
 
# %%
def GO(target_list, cutoff, term, title):
    key_set = {"mf":"GO_Molecular_Function_2021",
               "bp":"GO_Biological_Process_2018", 
               "react":"Reactome_2022",
               "Msig":"MSigDB_Hallmark_2020",
               "cc":"GO_Cellular_Component_2021"}

    import gseapy as gp
    enr = gp.enrichr(gene_list=target_list,
                     gene_sets=key_set[term],  #GO_Biological_Process_2023, MSigDB_Hallmark_2020, GO_Molecular_Function_2023, GO_Cellular_Component_2023
                     organism="human", # don't forget to set organism to the one you desired! e.g. Yeast
                     outdir=DIR+"GO_{}".format(title), # don't write to disk
                     cutoff=cutoff)
    
    input = enr.results[enr.results["Adjusted P-value"]<cutoff][["Term", "Adjusted P-value","Combined Score","Genes"]]
    input["-log(q)"] = -np.log10(input["Adjusted P-value"])
    input["Term"] = input["Term"].str.split("(").str[0]
    plt.figure(figsize=(5,1.5))
    plt.scatter(y="-log(q)", x="Term", 
                data=input.sort_values(by="-log(q)", ascending=False).iloc[:15,:],
                c="Combined Score",
                s=80)
    plt.ylabel("-Log10(FDR)")
    # plt.yticks(fontsize=8)
    plt.xticks(rotation=90)
    sns.despine()
    plt.savefig(DIR+f"figure/{term}_{title}2023.pdf",
                bbox_inches="tight")

    return input, key_set[term]


## Set Specific GOterm
terms = "react"
result, title = GO(data["gene"], 0.05, terms, "{}_decode".format(terms))
result2, title = GO(data2["gene"], 0.05, terms, "{}_rmat".format(terms))

# %%
## Get an information for uniquely captured or shared terms
SD = list(set(result["Term"].tolist()) - set(result2["Term"].tolist()))
RM = list(set(result2["Term"].tolist()) - set(result["Term"].tolist()))
shared = list(set(result["Term"].tolist()) & set(result2["Term"].tolist()))
print(f"SpliceDecoder only:({len(SD)}){SD}")
print(f"rMATS only:({len(RM)}){RM}")
print(f"Shared:({len(shared)}){shared}")

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

def plot_venn(sets, labels=None, title=None, 
              label_fontsize=14,
              number_fontsize=20):
    """
    Plots a Venn diagram for two or three sets.
    
    Parameters:
        sets (list of set): A list containing two or three sets.
        labels (list of str): Labels for the sets. Default is None.
        title (str): Title of the plot. Default is None.
    """
    if len(sets) == 2:
        venn = venn2(sets, set_labels=labels,
                     alpha=.9)
    elif len(sets) == 3:
        venn = venn3(sets, set_labels=labels,
                     alpha=.9)
    else:
        raise ValueError("This function supports only two or three sets.")
    
    # Adjust font size for set labels
    if labels:
        for label in venn.set_labels:
            if label is not None:  # Check for empty labels
                label.set_fontsize(label_fontsize)
    
    # Adjust font size for numbers inside the Venn diagram
    for subset in venn.subset_labels:
        if subset is not None:  # Check for empty subsets
            subset.set_fontsize(number_fontsize)

    if title:
        plt.title(title)


# Plot a Venn diagram for two sets
plot_venn([set(result["Term"].tolist()), set(result2["Term"].tolist())], labels=["SpliceDecoder\nonly", "rMATS\nonly"], title=f"{title}")
plt.savefig(DIR+f"figure/{data_type}_{title}_venn_diagram.pdf",
            bbox_inches="tight")
# %%
for i in SD:
    print(i)
# %%
genes = []
for i in result["Genes"].tolist():
    for gene in (i.split(";")):
        genes.append(gene)
genes = list(set(genes))
print(genes)

# %%
## Comp GO result between original rMATS vs filtered result
def Comp_GO(a,b,title):
    result["method"] = a
    result2["method"] = b
    shared_term = list(set(result["Term"].tolist()) & set(result2["Term"].tolist()))
    rmat_only_term = list(set(result2["Term"].tolist()) - set(result["Term"].tolist()))
    SD_only_term = list(set(result["Term"].tolist()) - set(result2["Term"].tolist()))
    inputa = result[["method","Term","-log(q)","Combined Score"]].set_index("Term")
    inputa = inputa.reindex(shared_term)
    inputa = inputa.sort_values("-log(q)",
                                ascending=False)
    ###
    ## To minimize figure size
    top15 = inputa.index.tolist()[:15]
    ###
    
    inputb = result2[["method","Term","-log(q)","Combined Score"]].set_index("Term")
    inputa = inputa.reindex(top15)
    inputb = inputb.reindex(top15)
    
    ### Sort by adjp and p_diff
    merged_df2 = pd.merge(inputa, inputb,
                          left_index=True,
                          right_index=True)
    merged_df2["p_diff"] = np.abs(merged_df2["-log(q)_x"] / merged_df2["-log(q)_y"])
    merged_df2 = merged_df2.sort_values(by=["-log(q)_x","p_diff"],
                                        ascending=[False,False])
    merged_df2 = merged_df2.iloc[:10,:]

    ## Apply sorted order
    merged_df = pd.concat([inputa.reindex(merged_df2.index), inputb.reindex(merged_df2.index)])
    merged_df["Term"] = merged_df.index.str.split("R-HSA").str[0]

    ### TEMP
    # figure,ax = plt.subplots(1,2,
    #                          figsize=(2.5,4),
    #                          sharey=True,
    #                          width_ratios=[2, 1])
    # plt.figure(figsize=(1.8,4))
    plt.figure(figsize=(6,1.8))
    custom_pal = {a:"#FF7E7E",
                  b:"#9F9F9F"}
    sns.scatterplot(y="-log(q)", x="Term", 
                    data=merged_df, hue="method",
                    palette=custom_pal, s=100,
                    edgecolor=None,)
                    # ax=ax[0])
    # sns.barplot(x="p_diff",
    #             y=merged_df2.index,
    #             data=merged_df2,
    #             width=.5,
    #             ax=ax[1],
    #             color="#46C4FF")
    # sns.despine(ax=ax[0])
    # sns.despine(ax=ax[1],
    #             left=False)
    sns.despine()

    inc_terms = []  # Terms with increased significance
    dec_terms = []  # Terms with decreased significance
    for term in shared_term:
        sc = merged_df[(merged_df.index==term) &
                        (merged_df["method"]==a)]["-log(q)"].values
        mat = merged_df[(merged_df.index==term) &
                        (merged_df["method"]==b)]["-log(q)"].values
        if (sc - mat) > 0.8:
            inc_terms.append(term)
        elif (mat - sc) > 0.2:
            dec_terms.append(term)

    plt.xticks(fontsize=10,
               rotation=90,
               ha='right',
               rotation_mode='anchor')
    plt.legend(frameon=False,
               bbox_to_anchor=(1.02,1.02))
    plt.xlabel(title)
    plt.ylabel("$-Log_{10}$(q)")
    plt.savefig(DIR+f"figure/{data_type}_{terms}_{a}_vs_{b}.pdf",
                bbox_inches="tight")
    
    print("Shared: {}, MATS_specific: {}, Decoder_specific: {}".format(len(shared_term),
                                                                       len(list(set(result2["Term"].tolist()) - set(result["Term"].tolist()))),
                                                                       len(list(set(result["Term"].tolist()) - set(result2["Term"].tolist())))
                                                                       )
                                                                       )
    return shared_term, inc_terms, rmat_only_term
    
shared_term, inc_terms, dec_terms = Comp_GO("SpliceDecoder","rMATS", title)

# %%
## Make effect score dataframe per each gene
## gene_score_whole contains all genes' maximum effect score, gene_score contains NMD and DOA genes' maximum effect score
nf_gene = set(gene_score_whole.index.tolist()) - set(gene_score.index.tolist())
gene_score_whole = gene_score_whole.reindex(nf_gene)    # non-functional changes genes only
merged_gene = pd.concat([gene_score, gene_score_whole])

## Make Functional change figure in each term
decreased = []  # rMATS only # of functional genes
de_score = []   # rMATS only effect score
## Check effect score distribution of decreased terms (rMATS only term)
for term in dec_terms:
    RMATS = set(result2[result2["Term"]==term]["Genes"].str.split(";").iloc[0]) # Extract gene list from rMATS result
    score_temp = len(RMATS & set(gene_score.index.tolist())) / len(RMATS)   # Calculate a proportion of functionally changed genes
    decreased.append(score_temp)
    for i in merged_gene.reindex(RMATS).dropna()[0]:    # Collect all effect score of each term
        de_score.append(i)

increased=[]    # Increased # of functional genes
in_score = []   # Increased effect score
for term in inc_terms:
    RMATS = set(result2[result2["Term"]==term]["Genes"].str.split(";").iloc[0]) # Extract gene list from rMATS result to confirm the differences between two categories (rMATS and increasd)
    score_temp = len(RMATS & set(gene_score.index.tolist())) / len(RMATS)
    increased.append(score_temp)
    for i in merged_gene.reindex(RMATS).dropna()[0]:
        in_score.append(i)

## Make a figure for proportion of functionally changed genes
plt.figure(figsize=(5,3))
sns.kdeplot(decreased, label="rMATS only", fill=True, color="#9F9F9F")
sns.kdeplot(increased, label="Increased", fill=True, color="#FF7E7E")
plt.xlabel("Proportion of funtional change events")
plt.legend(frameon=False,
           bbox_to_anchor=(1.02,1.02))
sns.despine()
plt.savefig(OUT+f"{data_type}_{terms}_increased_term_kde.pdf",
            bbox_inches="tight")
print(ks_2samp(decreased, increased))

## Make a figure for Maximum effect score dist
plt.figure(figsize=(1.5,3))
input=(pd.DataFrame([de_score,in_score]).T)
sns.boxplot(input, showfliers=False, 
            palette=["#9F9F9F","#FF7E7E"],
            width=.6)
plt.legend(frameon=False)
plt.xticks([0,1],["rMATS only","Increased"])
plt.ylabel("Effect Score")
plt.xticks(rotation=90)
sns.despine()
plt.savefig(OUT+f"{data_type}_{terms}_increased_term_box.pdf",
            bbox_inches="tight")
print(mannwhitneyu(de_score, in_score))

# %%
## Draw splicing consequence of increased term and rMATS only term as a control
## Select an interesting term
def Direction(df):
    """ Assign directions for 2 possible classes

    Args:
        df (_type_): dataframe

    Returns:
        _type_: Inclusion case has 1 direction (same direction with dPSI) 
                other case has -1 direction (opposite direction with dPSI)
    """
    if df["Simulated_event"] == "EI" or \
       df["Simulated_event"] == "RI" or \
       df["Simulated_event"] == "Alt_A3SS" or \
       df["Simulated_event"] == "Alt_A5SS" or \
       df["Simulated_event"] == "MXE2":
        return 1
    else:   # ES, SI, Ori_A3/5SS, MXE1
        return -1
    

def Draw_Fig5(target_Term):
    ## Get all gene list for individual term
    cont_target = (result[result["Term"].str.startswith(target_Term)]["Genes"].str.split(";").iloc[0])
    print(result[result["Term"].str.startswith(target_Term)])
    ## Query whether those genes in scoring df (Whole)
    sub_data3 = data3[data3["gene"].isin(cont_target)]

    ## Filter out other genes
    for k,gene in enumerate(sub_data3["gene"].unique()):
        temp = sub_data3[(sub_data3["gene"]==gene)].copy()
        temp = temp[temp["Effect_Score"]==temp["Effect_Score"].max()]
        if k == 0:
            uniq_data3 = temp.copy()
        else:
            uniq_data3 = pd.concat([uniq_data3,temp])
    
    ## Remain only ORF1 DS pairs
    sub_data3 = uniq_data3[uniq_data3["ORF"]==1.0]
    sub_data3 = sub_data3.sort_values(by="Effect_Score",
                                      ascending=False)
    # sub_data3 = sub_data3[sub_data3["ORF"]==1.0]
    # sub_data3 = sub_data3[sub_data3["Effect_Score"]>0]
    
    ## Update
    input_df = pd.DataFrame()
    for pair in sub_data3["LongID"].unique():
        pairs = sub_data3[sub_data3["LongID"]==pair]
        max_effect = pairs[pairs["Effect_Score"]==pairs["Effect_Score"].max()]
        max_effect = max_effect.copy()

        if max_effect.shape[0] > 1: # For duplicates check, all duplicates has zero as their maximum
            # One of them should be selected as representative effect. To do this, TU was considered
            max_effect = max_effect[max_effect["Transcript_usage"]==max_effect["Transcript_usage"].max()]
        dir1 = max_effect.apply(Direction, axis=1)
        max_effect["Highly_affected_group"] = (int(dir1.values[0]) * max_effect["Delta_PSI"].values[0])
        max_effect["Highly_affected_group"] = max_effect["Highly_affected_group"].apply(lambda x: "MYC high" if x > 0 else "MYC low")
        input_df = pd.concat([input_df, max_effect])

    ## Make a figure
    plt.figure(figsize=(6,3))
    input_df = input_df[input_df["Effect_Score"]>0]
    # print(input_df.head(10))
    ax = sns.stripplot(x="LongID", 
                       y="Effect_Score",
                       hue="Highly_affected_group",
                       data=input_df,
                       size=8,
                       palette={"MYC low":"#CACACA",
                                "MYC high":"#FCC169"})

    x_ticks = ax.get_xticks()
    x_labels = ax.get_xticklabels()
    x_mapping = {label.get_text(): tick for label, tick in zip(x_labels, x_ticks)}

    sns.despine()
    plt.legend(bbox_to_anchor=(.95,1.0),
               loc="upper left",
               frameon=False)
    ax.set_xticks([])
    ax.set_xlabel("DS event")
    ax.set_title(target_Term)

    first = 0
    for i, row in input_df.iterrows():
        x_position = x_mapping[row['LongID']] + 10
        if x_position < 10:
            x_position += 1
        if target_Term == "Cell Cycle":
            if row["Highly_affected_group"] == "MYC high":
                if row['Effect_Score'] >= 0.07:
                    plt.text(x=x_position+first,
                             y=row['Effect_Score'], 
                             s=row['gene'], 
                             horizontalalignment='right', 
                             fontsize=12, 
                             color='black')
                    first += 1
            else:
                pass
            
        else:
            if row['Effect_Score'] >= 0.1 and first < 4:
                plt.text(x=x_position, 
                         y=row['Effect_Score'], 
                         s=row['gene'], 
                         horizontalalignment='right', 
                         fontsize=12, 
                         color='black')
                first += 1
    
    plt.ylim(0,input_df["Effect_Score"].max()+0.05)
    plt.savefig(OUT+f"{target_Term}_effect_overview.pdf",
                bbox_inches="tight")

# Draw_Fig5("RNA Polymerase II Transcription R-")
Draw_Fig5("Cell Cycle")
# Draw_Fig5("Gene Expression")
Draw_Fig5("mRNA Splicing - Major")
# Draw_Fig5("cellular response to DNA damage stimulus")

# %%
print(gene_score.loc["CHEK2"])

# %%
## For check DS type
for gene in ["REV3L","LRIF1","CHEK2","ATRIP"]:
    max_val = data3[data3["gene"]==gene]["Effect_Score"].max()
    print(data3[data3["Effect_Score"]==max_val])
# %%
