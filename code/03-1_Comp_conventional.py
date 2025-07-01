#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import argparse
from scipy.stats import mannwhitneyu, pearsonr, gaussian_kde, fisher_exact, spearmanr, kruskal, ttest_ind, entropy
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# plt.rcParams["font.family"] = "Arial"

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''

Description
    #########################################################
    This script investigate an association between rMATS value and
    functional changes.
    OUTPUT: FDR/dPSI_asso.pdf, Jointplot_for_FDR-DI.pdf
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains simulated results', 
                        required=True,
                        type=str)
    
    # ## Optional
    # parser.add_argument("--threads", "-t", help="number of threads to use, default: 9", 
    #                     type=int,
    #                     default="9")

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args


def Load(stype):
    tl_data = pd.read_csv(args.input+"result/all_{}_Main_output.txt".format(stype),
                            sep="\t", 
                            skiprows=1)
    tl_data["Ref_CDS"] = (tl_data["Stop"].str.split("-").str[0]).astype(int) - \
                         (tl_data["Start"].str.split("-").str[0]).astype(int)
    tl_data = tl_data[["LongID","Target_TX","Start","3'UTR","dAA","Ref_CDS","ORF_priority","Sim_domain_length","Ref_domain_length","Domain_integrity","Domain_change_ratio","pNMD","occurred_event"]]
    rmat = pd.read_csv(args.input+"rmat.csv",
                       sep=",")
    rmat = rmat[["long_ID","dPSI_2_minus_1","FDR"]]
    merged_df = pd.merge(tl_data, rmat, left_on="LongID", right_on="long_ID")
    merged_df = merged_df[merged_df["Domain_integrity"]!="Frame loss"]  # Filter Frame loss cases
    merged_df = merged_df.dropna(subset="Domain_integrity") # Nan represents ref and sim has no domain
    merged_df["FDR"] = merged_df["FDR"].apply(lambda x : (merged_df[merged_df["FDR"]>0]["FDR"].min()) if x == 0 else x)
    merged_df["FDR"] = -np.log10(merged_df["FDR"])
    
    return merged_df


# %%
def Corr(df,target,target2,pal,axs):
    """ Compare traditional value (tv) and new value (nv)

    Args:
        df (_type_): merged_df (contains tv and nv)
        target (_type_): FDR or dPSI
        target2 (_type_): NMD or domain integ
        axs (_type_): subplots axs number
    """

    df = df[[target,"value"]].astype(float)
    df = df.drop_duplicates()
    if target2 == "NMD":
        df = df[df["value"]!=0]
        # df = df[df["value"]==0]   # To check the neutral category
        df["value"] = df["value"].apply(lambda x : "NMD" if x < 0 else ("PTC remove" if x > 0 else "Neutral"))
    else:
        df["value"] = df["value"].astype(float)
        df = df[df["value"]!=100]
        df = df[df["value"]!=np.inf]
        # df = df[df["value"]==100] # To check the neutral category
        df["value"] = df["value"].apply(lambda x : "LOD" if x < 100 else ("GOD" if x > 100 else "Neutral"))
    df[target] = np.abs(df[target])
    
    bin = []
    if target == "FDR":
        for i in (np.linspace(1,15,num=15)):
            bin.append(-np.log10(1/np.power(10,i)))
    else:
        for i in range(1,11):
            bin.append(i/10)

    df["bin"] = pd.cut(df[target],bin)
    ## Count the total DS in each bin and calculate proportion of events in each bin
    sns.countplot(y="bin", hue="value", 
                  data=df, ax=axs,
                  palette=pal)    # Y axis = mean of "values" in each bin
    axs.legend("",
               frameon=False)
    sns.despine(top=True, right=True, ax=axs)
    
    if target == "FDR":
        axs.set_ylabel("-Log10{}".format(target), fontsize=13)
    else:
        axs.set_ylabel("|{}|".format("dPSI"), fontsize=13)
    axs.set_xlabel("")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/"):
        args.input = args.input+"/"
#    interesting_target = args.interesting_target
    
    ## Run code
    OUT = args.input+"figure/"
    if not os.path.exists(OUT):
        os.mkdir(OUT)

    ## Do analysis
    k = 0
    for splice_type in ["CA", "RI", "A3SS", "A5SS", "MXE"]:
        try:
            sub_df = Load(splice_type)
            if k == 0:
                merged_df = sub_df
            else:
                merged_df = pd.concat([merged_df, sub_df])
            k += 1

        except:
            pass

    merged_df[["pNMD","Domain_integrity"]] = merged_df[["pNMD","Domain_integrity"]].astype(float)
   # %%
    std = {}
    for ds in merged_df["LongID"].unique():
        sub_pair = merged_df[merged_df["LongID"]==ds].copy()
        for stype in sub_pair["occurred_event"].unique():
            sub_class = sub_pair[sub_pair["occurred_event"]==stype]
            if sub_class.shape[0] > 1 and \
            sub_class["Domain_change_ratio"].astype(float).sum() != 0:
                std[ds+"_"+stype] = []
                std[ds+"_"+stype].append(np.std(sub_class["Domain_change_ratio"].astype(float)))
                std[ds+"_"+stype].append(sub_class.shape[0])
    # %%
    std_df = pd.DataFrame.from_dict(std, orient="index")
    for k,over10 in enumerate(std_df[1].value_counts()[std_df[1].value_counts() > 10].index):
        if k == 0:
            robust = std_df[std_df[1]==over10]
        else:
            robust = pd.concat([robust,std_df[std_df[1]==over10]])
    # %%
    plt.figure(figsize=(6,4))
    PROPS = {'boxprops':{'facecolor':'none', 'edgecolor':'grey'},
            'medianprops':{'color':'black', 'linewidth':3},
            'whiskerprops':{'color':'grey'},
            'capprops':{'color':'grey'}}

    ax = sns.boxplot(x=1, y=0, data=robust, showfliers=False,
                **PROPS)
    ylims=ax.get_ylim()
    ax = sns.stripplot(x=1, y=0, data=robust,
                    s=2, edgecolor=None, alpha=.8)
    ax.set(ylim=ylims)
    sns.despine()
    plt.xlabel("# of Matched transcripts")
    plt.ylabel("STD of Delta L")
    plt.savefig(OUT+"Difference_of_change_ratio.pdf",
                bbox_inches="tight")
    # %%
    ## dAA = 1 represents no CDS changes
    ## Check unknown CDS changes within UTR alteration category
    print(merged_df[(merged_df["pNMD"].astype(float)!=-1) & # Not NMD
                    (merged_df["dAA"].astype(float)<1) &    # CDS changes
                    (merged_df["Domain_change_ratio"].astype(float)==0)])    # No domain changes
    # %%
    from scipy import stats
    input = merged_df[(merged_df["Domain_change_ratio"].astype(float)!=1) &
                    (merged_df["Domain_change_ratio"].astype(float)!=0)].copy()
    plt.figure(figsize=(4,4))
    values = np.vstack([input["FDR"].astype(float), input["Domain_change_ratio"].astype(float)])
    kernel = stats.gaussian_kde(values)(values)
    sns.scatterplot(x=input["dPSI_2_minus_1"].astype(float), 
                    y=input["FDR"].astype(float),
                    data=merged_df,
                    c = np.log(input["Domain_change_ratio"].astype(float)),
                    s=5,
                    edgecolor=None)
    
    figfdr,axfdr = plt.subplots(1,2,figsize=(8,3))
    figfdr.subplots_adjust(wspace=0.0)
    pal = {"GOD":"#AABEDB", "LOD":"#487AC5", "Neutral":"#AEAEAE"}
    pal2 = {"PTC remove":"#DBAAAA", "NMD":"#C54848", "Neutral":"#AEAEAE"}
    nmd = merged_df.copy()
    Corr(merged_df[["FDR","Domain_integrity"]].melt(id_vars="FDR"),"FDR","Domain_integrity",pal,axfdr[0])
    Corr(nmd[["FDR","pNMD"]].melt(id_vars="FDR"),"FDR","NMD",pal2,axfdr[1])
    axfdr[0].invert_xaxis()
    axfdr[1].set_yticklabels("")
    axfdr[1].set_ylabel("")
    axfdr[1].tick_params(left=False)
    handle = [matplotlib.patches.Patch(color='#487AC5', label='LOD'),
            matplotlib.patches.Patch(color='#AABEDB', label='GOD'),
            matplotlib.patches.Patch(color='#C54848', label='NMD'), 
            matplotlib.patches.Patch(color='#DBAAAA', label='PTC remove')]
    figfdr.text(0.5, -0.02, '# of Functional event', 
                va='center', ha='center',
                fontsize=13)
    plt.savefig(OUT+"{}_asso.pdf".format("FDR"),
                    bbox_inches="tight")

    figpsi,axpsi= plt.subplots(1,2,figsize=(8,3))
    figpsi.subplots_adjust(wspace=0.0)
    Corr(merged_df[["dPSI_2_minus_1","Domain_integrity"]].melt(id_vars="dPSI_2_minus_1"),"dPSI_2_minus_1","Domain loss",pal,axpsi[0])
    Corr(nmd[["dPSI_2_minus_1","pNMD"]].melt(id_vars="dPSI_2_minus_1"),"dPSI_2_minus_1","NMD",pal2,axpsi[1])
    axpsi[0].invert_xaxis()
    axpsi[1].set_yticklabels("")
    axpsi[1].set_ylabel("")
    axpsi[1].tick_params(left=False)
    axpsi[1].legend(handles=handle,
                    frameon=False,
                    ncol=2,
                    bbox_to_anchor=(1.05,0.2))
    figpsi.text(0.5, -0.02, '# of Functional event', 
                va='center', ha='center',
                fontsize=13)
    plt.savefig(OUT+"{}_asso.pdf".format("PSI"),
                    bbox_inches="tight")
    # %%
    ## Supp figure to show association between FDR and DI 
    plt.figure(figsize=(5,5))
    non_zero = merged_df[(merged_df["Domain_integrity"]!=100) &
                        (merged_df["Domain_integrity"]!=np.inf)].copy()
    non_zero["Domain_integrity"] = np.log10(((non_zero["Domain_integrity"])/100.0))
    non_zero = non_zero[(non_zero["Domain_integrity"] <= 1.5) &
                        (non_zero["Domain_integrity"] >= -1.5)]
    sns.jointplot(x="FDR", 
                y="Domain_integrity",
                data=non_zero,
                kind="reg", scatter_kws={'s':3, "alpha":.3},
                joint_kws={'line_kws': {'linewidth': 3, "color":"red"}},
                height=5, ratio=2)
    plt.xlabel("$-log_{10}FDR$")
    plt.ylabel("$log_{10}DI$")
    plt.savefig(OUT+"Jointplot_for_FDR-DI.pdf",
                bbox_inches="tight")

# %%
def Length_comp():
    ## This plot showed number of matched transcript of certain splicing events
    import datashader as ds
    from datashader.mpl_ext import dsshow
    from scipy.stats import gaussian_kde
    # https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density

    major_length = pd.read_csv("/Users/kangh/Documents/Analysis_note/isoform_translation/MYC_paper/canonical_gencode_length.txt",
                            sep="\t")
    minor_length = pd.read_csv("/Users/kangh/Documents/Analysis_note/isoform_translation/MYC_paper/non-cano_gencode_length.txt",
                            sep="\t")
    # %%
    figure,ax = plt.subplots(3,1, figsize=(4,4), sharex=True)
    sns.kdeplot(np.log10(major_length["1"]),
                label="Canonical",
                color="#48C550",
                linewidth=4, ax=ax[0])
    sns.kdeplot(np.log10(minor_length["1"]),
                label="Non-canonical",
                color="#92C295",
                linewidth=4, ax=ax[0])
    ax[0].legend(frameon=False, bbox_to_anchor=(.8, 1.02), loc="upper left")
    ax[0].set_ylabel('')

    input_df = merged_df.copy()
    input_df["gene"] = input_df["long_ID"].str.split(";").str[1]
    input_df["pNMD"] = input_df["pNMD"].astype(int).apply(lambda x : "NMD" if x == -1 else ("Recovered" if x == 1 else "No NMD"))
    recov_df = input_df[input_df["pNMD"]=="Recovered"] # recov case
    nmd_df = input_df[input_df["pNMD"]=="NMD"] # NMD case
    sns.kdeplot(np.log10(nmd_df["Ref_CDS"]+1),
                label="NMD",
                color="#C54848", ax=ax[1],
                linewidth=4)
    sns.kdeplot(np.log10(recov_df["Ref_CDS"]+1),
                label="PTC remove",
                color="#DBAAAA", ax=ax[1],
                linewidth=4)
    ax[1].set_ylabel('')
    plt.xlabel("$Log_{10}$(Ref TX CDS)")
    sns.despine()
    ax[1].legend(frameon=False, bbox_to_anchor=(.8, 1.02), loc="upper left")
                
    lod = input_df[input_df["Domain_integrity"]<100]
    god = input_df[input_df["Domain_integrity"]>100]
    sns.kdeplot(np.log10(lod["Ref_CDS"]+1),
                label="LOD",
                color="#487AC5",
                ax=ax[2], linewidth=4)
    sns.kdeplot(np.log10(god["Ref_CDS"]+1),
                label="GOD",
                color="#AABEDB",
                ax=ax[2], linewidth=4)
    ax[2].set_ylabel('')
    ax[2].set_xlabel("$Log_{10}$(Ref TX CDS)")
    sns.despine()
    ax[2].legend(frameon=False, bbox_to_anchor=(.8, 1.02),
            loc="upper left")
    print(ttest_ind(np.log10(god["Ref_CDS"]+1), np.log10(lod["Ref_CDS"]+1)))
    print(mannwhitneyu(np.log10(god["Ref_CDS"]+1), np.log10(lod["Ref_CDS"]+1)))
    plt.savefig(OUT+"event_freq_length.pdf",
                bbox_inches="tight")

# %%
