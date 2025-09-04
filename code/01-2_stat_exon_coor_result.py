#!/usr/bin/env python3
# %%
""" Check the mapping for stats between splicing events to the transcriptome

Returns:
    _type_: Some of basic figures and input for 01-4*py
"""
import pandas as pd
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.rcParams["font.family"] = "Arial"

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
Description
    #########################################################
    This script makes summary stats of mapping process and merged input file.
    OUTPUT: mapping.stats, mat_tx_numbers.pdf, mapping_rate.pdf, fig_input.txt
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat file', 
                        required=True,
                        type=str)

    ## Optional
    parser.add_argument('--target', '-t', 
                        help='List of interesting gene IDs (ENSG) file (tsv)', 
                        required=False,
                        default="all",
                        type=str)

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

def Draw_dist(data, axs):
    for col in ["first_diff", "second_diff", "third_diff", "fourth_diff"]:
        data[col] = data[col].replace("SKIP",1e6).replace("RI",1e6).replace("None",1e6).copy()
        data[col] = data[col].astype(float).copy()

    filtered_df = data[data["whole_diff"].astype(int)==0]   # Consider perfectly matched TXs only
    n_tx = []
    if filtered_df.shape[0] > 0:
        for sp in filtered_df["ID"].unique():
            n_tx.append(filtered_df[filtered_df["ID"]==sp].shape[0])    # Number of perfectly matched TXs for each splicing event
    
        sns.kdeplot(n_tx, ax=axs, fill=True,
                    color=col_code[filtered_df["Event_type"].unique()[0]],
                    alpha=1.0)
        sns.despine()
        axs.set_xlim(0,20)

    return filtered_df, n_tx


def Check_Unmapped(data):
    unmapped_df = pd.DataFrame(columns=data.columns)
    for SP in data["ID"].unique():
        pairs = data[data["ID"]==SP]
        if not 0 in pairs["whole_diff"].unique():
            unmapped_df = pd.concat([unmapped_df,pairs])

    unmapped_df.to_csv(args.input+"unmapped.txt",
                       sep="\t",
                       index=None)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/"):
        args.input = args.input+"/"
    rmat = pd.read_csv(args.input+"rmat.csv",
                       sep=",")
    geneID = "geneID"
    
    if args.target != "all":
        target_Genes = pd.read_csv(f"{args.target}",
                                   sep="\t",
                                   header=None)
        rmat = rmat[rmat[geneID].isin(target_Genes[0])]
    
    else:
        pass

        
    ## Run code
    OUT = args.input+"figure/"
    if not os.path.exists(OUT):
        os.mkdir(OUT)

    col_code = {"CA": sns.color_palette("tab10")[0],
                "RI": sns.color_palette("tab10")[1],
                "A3SS": sns.color_palette("tab10")[2],
                "A5SS": sns.color_palette("tab10")[3],
                "MXE": sns.color_palette("tab10")[4]}

    multi_data = pd.read_csv(args.input+"long_output.txt", sep="\t")
    short_data = pd.read_csv(args.input+"short_output.txt", sep="\t")
    merged = pd.concat([multi_data, short_data])
    ######### Temporal #########
#    merged["Event_type"] = "CA"   # To estimate mapping accuracy
    merged["Event_type"] = merged["ID"].str.split(";").str[0]
    ######### Temporal #########
    Check_Unmapped(merged)
    # %%
    mp_array = []
    figure,ax = plt.subplots(5,1,figsize=(4,6),sharex=True)
    for k,SE in enumerate(merged["Event_type"].unique()):
        filtered_merged, n_tx = Draw_dist(merged[merged["Event_type"]==SE], ax[k])
        total_event = len(rmat[rmat["Event_type"]==SE]["long_ID"].unique())
        mapping_rate = len(filtered_merged["ID"].unique()) / total_event
        mp_array.append(mapping_rate*100)

        if k == 0:
            merged_f_data = filtered_merged
        else:
            merged_f_data = pd.concat([merged_f_data,filtered_merged])
    figure.supxlabel("# of matched transcripts")
    plt.savefig(args.input+"figure/mat_tx_numbers.pdf",
                bbox_inches="tight")
    # %%
    mp = pd.DataFrame(mp_array, 
                    index=merged["Event_type"].unique().tolist(),
                    columns=["Mapping rate (%)"])
    plt.figure(figsize=(3,6))
    sns.barplot(x=mp.index, y=mp["Mapping rate (%)"],
                palette=col_code,
                order=merged["Event_type"].unique())
    sns.despine()
    plt.title("Avg mapp rate: {}%".format(round(np.mean(mp)),3))
    plt.savefig(args.input+"figure/mapping_rate.pdf",
                bbox_inches="tight")
    # %%
    ## Make final input AND return total mapping rate
    mr = len(merged_f_data[merged_f_data["whole_diff"]==0]["ID"].unique()) / len(rmat["long_ID"].unique())
    print(f"Total mapping rate: {round(mr*100,3)}%")
    final = merged_f_data[(merged_f_data["whole_diff"] == 0)]
    print(f"Mapped DS events: {len(final['ID'].unique())}")
    final.to_csv(args.input+"fig_input.txt", sep="\t", index=None)
    mp.to_csv(args.input+"mapping.stats",
              sep="\t")
