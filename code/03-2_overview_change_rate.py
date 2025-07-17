#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import sys
import os
import argparse
from scipy.stats import pearsonr, mannwhitneyu
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# plt.rcParams["font.family"] = "Arial"

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
 
Description
    #########################################################
    This script draws transcript and domain level differences
    OUTPUT: splicing_categories_stacked_plot.pdf, merged_stacked_plot.pdf
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat result', 
                        required=True,
                        type=str)
    
    ## Optional
    parser.add_argument('--organism', '-or', 
                        help='human or mouse, default: human', 
                        type=str,
                        default="human")
    parser.add_argument("--orf_filter", "-orf", help="Only consider ORF1, default: Y", 
                        type=str,
                        default="Y")

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args


def Load(stype):
    """
    It adds the several types of ID for the given dataframe
    And remove duplicates based on this ID, pNMD, and ORF
    """
    ref_result = pd.read_csv(args.input+f"all_{stype}_Main_output.txt",
                            sep="\t", 
                            skiprows=1)
    ref_result.columns = ["LongID","Gene symbol","ref_tx","event","ORF","Start","Stop","5'UTR","dAA","3'UTR","Domain_integrity","change_rate","Ref_Domain","Sim_Domain","DOA_direction","pNMD"]
    ref_result["pair_ID"] = ref_result["LongID"] + "|" + ref_result["ref_tx"]
    ref_result["ComplexID"] = ref_result["LongID"] + "|" + ref_result["ref_tx"] + "|" + ref_result["event"]
    # ref_result["gene"] = ref_result["LongID"].str.split(";").str[1]
    ref_result_nmd = ref_result[["pNMD","pair_ID","ORF"]].drop_duplicates()

    return ref_result


def Add_UTR(df):
    if df["3'UTR"] != 0 and df["pNMD"] == "0":
        return "UTR"
    elif df["5'UTR"] != 0 and df["pNMD"] == "0":
        return "UTR"
    else:
        return df["pNMD"]
    
def Add_UTR2(df):
    if df["change_rate"] == "inf":  # PTC remove
        return "Gain of domain (GOD)"
    
    elif df["change_rate"] == "Frame loss":
        return "Frame loss"
    
    elif float(df["3'UTR"]) != 0 and float(df["change_rate"]) == 0.0:
        if float(df["dAA"]) == 0:  # No difference in CDS length (Stop - Start)
            return "UTR alt"
        else:
            return "CDS alt"
    
    elif float(df["5'UTR"]) != 0 and float(df["change_rate"]) == 0.0:
        if float(df["dAA"]) == 0:
            return "UTR alt"
        else:
            return "CDS alt"
    
    elif float(df["change_rate"]) > 0.0 and \
         (float(df["Ref_Domain"]) - float(df["Sim_Domain"])) < 0:
        return "Gain of domain (GoD)"
    
    elif float(df["change_rate"]) > 0.0 and \
         (float(df["Ref_Domain"]) - float(df["Sim_Domain"])) > 0:
        return "Loss of domain (LoD)"
    
    else:   # Unknown cases, 100 DI and non-zero change_rate
        return "CDS alt"    # Gain some domain as much as loss some domain
    

# %%
## Make binary class (functional vs non-functional)
def Cal(df):
    if "NMD" == df["pNMD"] and \
       "PTC remove" == df["pNMD"]:
        return 1
    else:
        return 0
    
def Cal2(df):
    if not "UTR alt" in df["change_rate"] and \
       not "Frame loss" in df["change_rate"]:
        return 1
    else:
        return 0

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/result/"):
        args.input = args.input+"/result/"
    OUT = args.input+"figure/"
    if not os.path.exists(OUT):
        os.mkdir(OUT)
    # args.organism = organism
    
    ## Run code
    k = 0
    for splice_type in ["CA", "RI", "A3SS", "A5SS", "MXE"]:
        try:
            sub_df = Load(splice_type)
            if k == 0:
                ref_result = sub_df
            else:
                ref_result = pd.concat([ref_result, sub_df])    # Merged Main_output of all splicing categories
            k += 1

        except:
            pass


    ## ORF filter
    if args.orf_filter == "Y":
        ref_result = ref_result[ref_result["ORF"]=="pORF1"]
    else:
        pass
    
    ref_result.to_csv(args.input+"merged_input_for_eval.txt",
                      sep="\t",
                      index=None)   # Just merged df of Load function
    print(f"\n\nYour Scoring input is ready!\nJust run: python 03-3_scoring_function.py -i {'/'.join(args.input.split('/')[:-2])} -t Y")

    ## Make a general stats figure
    stats = ref_result.copy()
    stats["pNMD"] = stats.apply(Add_UTR, axis=1)
    stats["Key"] = stats["LongID"].str.split(";").str[0]    # DS ID
    stats["pNMD"] = stats["pNMD"].apply(lambda x: "NMD" if x == "-1" else ("PTC remove" if x == "1" else x))
    stats = stats.dropna(subset="Domain_integrity") # Nan represents ref and sim has no domain
    stats["change_rate"] = stats.apply(Add_UTR2, axis=1)
    stats = stats[["LongID","Key","pNMD","Gene symbol","Domain_integrity","change_rate","ref_tx","ORF","event"]]
    stats["NMD"] = stats.apply(Cal, axis=1)
    stats["domain"] = stats.apply(Cal2, axis=1)
    fin = stats[["LongID","NMD","domain"]]

    sp_we = []  # With functional events
    sp_woe = [] # Without functional events
    for i in fin["LongID"].unique():
        temp_fin = fin[fin["LongID"]==i].sum()  # Find at least one functional event
        if temp_fin["NMD"] > 0 or \
        temp_fin["domain"] > 0: # Except for UTR of Neutral
            sp_we.append(i)
        else:   # The splicing event never have at least functional event
            sp_woe.append(i)

    input_df = np.array([len(sp_we), len(sp_woe)]) / len(fin["LongID"].unique())
    # plt.figure(figsize=(0.5,4))
    # sns.barplot(input_df[[1]]+input_df[[0]],
    #             color="#B8D2FF",
    #             label="No Functional change")
    # sns.barplot(input_df[[0]], color="#3E85FF",
    #             label="At least one Functional change")
    # sns.despine()
    # plt.xticks([])
    # plt.ylabel("Proportion of DS")
    # plt.savefig(OUT+"Fig3_non-functional_DS.pdf",
    #             bbox_inches="tight")


    ## Make stats for NMD, Frame loss, and Frame Gain
    new_columns = []
    for k,event in enumerate(stats["event"].unique()):    # Key: each splicing categories
        temp = stats[stats["event"]==event]
        whole_number = temp.shape[0]    # DOA & NMD & UTR
        nmd = pd.DataFrame(temp["pNMD"].value_counts()) # number of NMD
        nmd = nmd[(nmd.index=="NMD") |
                  (nmd.index=="PTC remove")]
        doa = temp[(temp["pNMD"]!="NMD") &
                   (temp["pNMD"]!="PTC remove")]    # Filter out NMD cases
        doa = pd.DataFrame(doa["change_rate"].value_counts())
        temp = pd.concat([nmd, doa])
        temp = temp/whole_number
        
        if "_" in event:
            event = event.replace("_", " ")
        temp.columns = [event]
        if k == 0:
            input = temp.copy()
        else:
            input = pd.merge(input,temp, 
                             left_index=True, 
                             right_index=True,
                             how="outer")
    input = input.fillna(0.0)
   
    figure, axs = plt.subplots(1,len(stats["event"].unique()), 
                               figsize=(6,6), sharey=True)
    figure.subplots_adjust(wspace=.6)
    for k,event in enumerate(input.columns):
        try:
            nmd = input[event].loc["NMD"]
        except:
            nmd = 0
        
        try:
            ptc_remove = input[event].loc["PTC remove"]
        except:
            ptc_remove = 0
        
        try:
            utr = input[event].loc["UTR alt"]
        except:
            utr = 0

        try:
            cds = input[event].loc["CDS alt"]
        except:
            cds = 0

        try:
            lod = input[event].loc["Loss of domain (LoD)"]
        except:
            lod = 0
        
        try:
            god = input[event].loc["Gain of domain (GoD)"]
        except:
            god = 0

        try:
            fl = input[event].loc["Frame loss"]
        except:
            fl = 0

        # Stacked bar chart
        axs[k].bar([0], nmd, color="#C54848")
        axs[k].bar([0], ptc_remove, color="#DBAAAA", bottom = nmd)
        axs[k].bar([0], lod, color="#487AC5", bottom = nmd+ptc_remove)
        axs[k].bar([0], god, color="#AABEDB", bottom = nmd+ptc_remove+lod)
        axs[k].bar([0], cds, color="#7B7B7B", bottom = nmd+ptc_remove+lod+god)
        axs[k].bar([0], utr, color="#B0B0B0", bottom = nmd+ptc_remove+lod+god+cds)
        axs[k].bar([0], fl, color="#DEDEDE", bottom = nmd+ptc_remove+lod+god+cds+utr)
        axs[k].set_xticklabels("", rotation=90)
        sns.despine(left=True)
        labels = ["NMD", "PTC remove", "LoD", "GoD", "CDS alt", "UTR alt", "FL"]
        plt.legend(labels=labels,
                    frameon=False,
                    bbox_to_anchor=(1.02, .99),
                    loc="upper left")
        axs[k].set_xlabel(event, rotation=90)    
    plt.savefig(OUT+"splicing_categories_stacked_plot.pdf",
                bbox_inches="tight")


    ## Merged stats
    temp = stats.copy()
    whole_number = temp.shape[0]    # DOA & NMD & UTR
    nmd = pd.DataFrame(temp["pNMD"].value_counts()) # number of NMD
    nmd = nmd[(nmd.index=="NMD") |
                (nmd.index=="PTC remove")]
    doa = temp[(temp["pNMD"]!="NMD") &
                (temp["pNMD"]!="PTC remove")]    # Filter out NMD cases
    doa = pd.DataFrame(doa["change_rate"].value_counts())
    temp = pd.concat([nmd, doa])
    temp = temp/whole_number
    input = temp.copy()
    input.columns = ["Merged"]
    event = "Merged"
    figure, axs = plt.subplots(1,1, 
                               figsize=(0.5,4))
    figure.subplots_adjust(wspace=.6)
    
    try:
        nmd = input[event].loc["NMD"]
    except:
        nmd = 0
    
    try:
        ptc_remove = input[event].loc["PTC remove"]
    except:
        ptc_remove = 0
    
    try:
        utr = input[event].loc["UTR alt"]
    except:
        utr = 0

    try:
        cds = input[event].loc["CDS alt"]
    except:
        cds = 0

    try:
        lod = input[event].loc["Loss of domain (LoD)"]
    except:
        lod = 0
    
    try:
        god = input[event].loc["Gain of domain (GoD)"]
    except:
        god = 0

    try:
        fl = input[event].loc["Frame loss"]
    except:
        fl = 0

    # Stacked bar chart
    axs.bar([0], nmd, color="#C54848")
    axs.bar([0], ptc_remove, color="#DBAAAA", bottom = nmd)
    axs.bar([0], lod, color="#487AC5", bottom = nmd+ptc_remove)
    axs.bar([0], god, color="#AABEDB", bottom = nmd+ptc_remove+lod)
    axs.bar([0], cds, color="#7B7B7B", bottom = nmd+ptc_remove+lod+god)
    axs.bar([0], utr, color="#B0B0B0", bottom = nmd+ptc_remove+lod+god+cds)
    axs.bar([0], fl, color="#DEDEDE", bottom = nmd+ptc_remove+lod+god+cds+utr)
    axs.set_xticklabels("", rotation=90)
    sns.despine(left=True)
    labels = ["NMD", "PTC remove", "LoD", "GoD", "CDS alt", "UTR alt", "FL"]
    plt.legend(labels=labels,
               frameon=False,
               bbox_to_anchor=(1.02, .99),
               loc="upper left")
    axs.set_xlabel(event)    
    plt.savefig(OUT+"merged_stacked_plot.pdf",
                bbox_inches="tight")
    # %%
