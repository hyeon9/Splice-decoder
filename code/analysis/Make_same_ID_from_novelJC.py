#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import sys
import os
import sys
import argparse

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
    Description
    #########################################################
    This script makes summary stats of mapping process and merged input file.
    OUTPUT: mapping.stats, mat_tx_numbers.pdf, mapping_rate.pdf, fig_input.txt
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--rMATS_folder', '-r',
                        help='Input directory that contains rmat file',
                        required=True,
                        type=str)
    parser.add_argument('--input', '-i',
                        help='Input directory that contains rmat file',
                        required=True,
                        type=str)
    
    
    ## Optional
    parser.add_argument('--geneid', '-id', nargs='*',
                        help='Column name of gene in GTF', 
                        type=str,
                        default=['Gene'])

    args = parser.parse_args(cmd_args, namespace)

    return args, parse_args

args, parser = parse_args(sys.argv[1:])
dir = args.rMATS_folder
outdir = args.input
files = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f)) if f.startswith("fromGTF.novelJunction")]

for header, matsfile in enumerate(files):
    event_type = matsfile.split(".")[2]
    df = pd.read_csv(f"{dir}{matsfile}",
                    sep="\t")
    if len(df[df.columns[2]].unique()) == 1:    # BUG fixed 250122
        geneid_type = df.columns[1]
    else:
        geneid_type = df.columns[2]

    df["Event_type"] = event_type
    df[[geneid_type]] = df[[geneid_type]].astype(str)   # BUG fixed 241120
    if event_type == "SE" or event_type == "RI":
        try:
            df["exonStart_0base"] += 1
        except:
            df["riExonStart_0base"] += 1
        df["upstreamES"] += 1
        df["downstreamES"] += 1
        pos_df = df[df["strand"]=="+"]
        neg_df = df[df["strand"]!="+"]
        if event_type == "SE":
            pos_df["long_ID"] = pos_df.apply(lambda row: ";".join(["CA", row[geneid_type], row["chr"], row["strand"], str(row["upstreamES"]), str(row["upstreamEE"]), str(row["exonStart_0base"]), str(row["exonEnd"]), str(row["downstreamES"]), str(row["downstreamEE"])]), axis=1)
            pos_df["E1s"] = pos_df["upstreamES"]
            pos_df["E1e"] = pos_df["upstreamEE"]
            pos_df["E2s"] = pos_df["exonStart_0base"]
            pos_df["E2e"] = pos_df["exonEnd"]
            pos_df["E3s"] = pos_df["downstreamES"]
            pos_df["E3e"] = pos_df["downstreamEE"]
            pos_df["E4s"] = None
            pos_df["E4e"] = None

            neg_df["long_ID"] = neg_df.apply(lambda row: ";".join(["CA", row[geneid_type], row["chr"], row["strand"], str(row["downstreamEE"]), str(row["downstreamES"]), str(row["exonEnd"]), str(row["exonStart_0base"]), str(row["upstreamEE"]), str(row["upstreamES"])]), axis=1)
            neg_df["E1e"] = neg_df["downstreamEE"]
            neg_df["E1s"] = neg_df["downstreamES"]
            neg_df["E2e"] = neg_df["exonEnd"]
            neg_df["E2s"] = neg_df["exonStart_0base"]
            neg_df["E3e"] = neg_df["upstreamEE"]
            neg_df["E3s"] = neg_df["upstreamES"]
            neg_df["E4e"] = None
            neg_df["E4s"] = None

            df = pd.concat([pos_df, neg_df])
            df["Event_type"] = "CA"

        else:   # RI
            pos_df["long_ID"] = pos_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["upstreamES"]), str(row["upstreamEE"]), str(row["upstreamEE"]+1), str(row["downstreamES"]-1), str(row["downstreamES"]), str(row["downstreamEE"])]), axis=1)
            pos_df["E1s"] = pos_df["upstreamES"]
            pos_df["E1e"] = pos_df["upstreamEE"]
            pos_df["E2s"] = pos_df["upstreamEE"]+1
            pos_df["E2e"] = pos_df["downstreamES"]-1
            pos_df["E3s"] = pos_df["downstreamES"]
            pos_df["E3e"] = pos_df["downstreamEE"]
            pos_df["E4s"] = None
            pos_df["E4e"] = None

            neg_df["long_ID"] = neg_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["downstreamEE"]), str(row["downstreamES"]), str(row["downstreamES"]-1), str(row["upstreamEE"]+1), str(row["upstreamEE"]), str(row["upstreamES"])]), axis=1)
            neg_df["E1e"] = neg_df["downstreamEE"]
            neg_df["E1s"] = neg_df["downstreamES"]
            neg_df["E2e"] = neg_df["downstreamES"]-1
            neg_df["E2s"] = neg_df["upstreamEE"]+1
            neg_df["E3e"] = neg_df["upstreamEE"]
            neg_df["E3s"] = neg_df["upstreamES"]
            neg_df["E4e"] = None
            neg_df["E4s"] = None

            df = pd.concat([pos_df, neg_df])

    elif event_type == "A3SS" or event_type == "A5SS":
        df["longExonStart_0base"] += 1
        df["shortES"] += 1
        df["flankingES"] += 1
        pos_df = df[df["strand"]=="+"]
        neg_df = df[df["strand"]!="+"]

        if event_type == "A3SS":
            pos_df["long_ID"] = pos_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["flankingES"]), str(row["flankingEE"]), str(row["longExonStart_0base"]), str(row["shortES"]-1), str(row["shortES"]), str(row["shortEE"])]), axis=1)
            pos_df["E1s"] = pos_df["flankingES"]
            pos_df["E1e"] = pos_df["flankingEE"]
            pos_df["E2s"] = pos_df["longExonStart_0base"]
            pos_df["E2e"] = pos_df["shortES"]-1
            pos_df["E3s"] = pos_df["shortES"]
            pos_df["E3e"] = pos_df["shortEE"]
            pos_df["E4s"] = None
            pos_df["E4e"] = None

            neg_df["long_ID"] = neg_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["flankingEE"]), str(row["flankingES"]), str(row["longExonEnd"]), str(row["shortEE"]+1), str(row["shortEE"]), str(row["shortES"])]), axis=1)
            neg_df["E1e"] = neg_df["flankingEE"]
            neg_df["E1s"] = neg_df["flankingES"]
            neg_df["E2e"] = neg_df["longExonEnd"]
            neg_df["E2s"] = neg_df["shortEE"]+1
            neg_df["E3e"] = neg_df["shortEE"]
            neg_df["E3s"] = neg_df["shortES"]
            neg_df["E4e"] = None
            neg_df["E4s"] = None

            df = pd.concat([pos_df, neg_df])
        
        else:
            pos_df["long_ID"] = pos_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["shortES"]), str(row["shortEE"]), str(row["shortEE"]+1), str(row["longExonEnd"]), str(row["flankingES"]), str(row["flankingEE"])]), axis=1)
            pos_df["E1s"] = pos_df["shortES"]
            pos_df["E1e"] = pos_df["shortEE"]
            pos_df["E2s"] = pos_df["shortEE"]+1
            pos_df["E2e"] = pos_df["longExonEnd"]
            pos_df["E3s"] = pos_df["flankingES"]
            pos_df["E3e"] = pos_df["flankingEE"]
            pos_df["E4s"] = None
            pos_df["E4e"] = None

            neg_df["long_ID"] = neg_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["shortEE"]), str(row["shortES"]), str(row["shortES"]-1), str(row["longExonStart_0base"]), str(row["flankingEE"]), str(row["flankingES"])]), axis=1)
            neg_df["E1e"] = neg_df["shortEE"]
            neg_df["E1s"] = neg_df["shortES"]
            neg_df["E2e"] = neg_df["shortES"]-1
            neg_df["E2s"] = neg_df["longExonStart_0base"]
            neg_df["E3e"] = neg_df["flankingEE"]
            neg_df["E3s"] = neg_df["flankingES"]
            neg_df["E4e"] = None
            neg_df["E4s"] = None

            df = pd.concat([pos_df, neg_df])

    elif event_type == "MXE":
        df["1stExonStart_0base"] += 1
        df["2ndExonStart_0base"] += 1
        df["upstreamES"] += 1
        df["downstreamES"] += 1
        pos_df = df[df["strand"]=="+"]
        neg_df = df[df["strand"]!="+"]

        pos_df["long_ID"] = pos_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["upstreamES"]), str(row["upstreamEE"]), str(row["1stExonStart_0base"]), str(row["1stExonEnd"]), str(row["2ndExonStart_0base"]), str(row["2ndExonEnd"]), str(row["downstreamES"]), str(row["downstreamEE"])]), axis=1)
        pos_df["E1s"] = pos_df["upstreamES"]
        pos_df["E1e"] = pos_df["upstreamEE"]
        pos_df["E2s"] = pos_df["1stExonStart_0base"]
        pos_df["E2e"] = pos_df["1stExonEnd"]
        pos_df["E3s"] = pos_df["2ndExonStart_0base"]
        pos_df["E3e"] = pos_df["2ndExonEnd"]
        pos_df["E4s"] = pos_df["downstreamES"]
        pos_df["E4e"] = pos_df["downstreamEE"]

        neg_df["long_ID"] = neg_df.apply(lambda row: ";".join([row["Event_type"], row[geneid_type], row["chr"], row["strand"], str(row["downstreamEE"]), str(row["downstreamES"]), str(row["2ndExonEnd"]), str(row["2ndExonStart_0base"]), str(row["1stExonEnd"]), str(row["1stExonStart_0base"]), str(row["upstreamEE"]), str(row["upstreamES"])]), axis=1)
        neg_df["E1e"] = neg_df["downstreamEE"]
        neg_df["E1s"] = neg_df["downstreamES"]
        neg_df["E2e"] = neg_df["2ndExonEnd"]
        neg_df["E2s"] = neg_df["2ndExonStart_0base"]
        neg_df["E3e"] = neg_df["1stExonEnd"]
        neg_df["E3s"] = neg_df["1stExonStart_0base"]
        neg_df["E4e"] = neg_df["upstreamEE"]
        neg_df["E4s"] = neg_df["upstreamES"]

        df = pd.concat([pos_df, neg_df])

    df["E1"] = df["chr"]+":"+df["E1s"].astype(str)+"-"+df["E1e"].astype(str)
    df["E2"] = df["chr"]+":"+df["E2s"].astype(str)+"-"+df["E2e"].astype(str)
    df["E3"] = df["chr"]+":"+df["E3s"].astype(str)+"-"+df["E3e"].astype(str)
    if event_type == "MXE":
        df["E4"] = df["chr"]+":"+df["E4s"].astype(str)+"-"+df["E4e"].astype(str)
    else:
        df["E4"] = None

    def cal_mean(input_df):
        try:    # For multisample comparison
            lst = input_df.split(",")
        except: # For 1:1 comparison
            lst = [input_df]
        lst = [0 if x == 'NA' else float(x) for x in lst]

        return sum(lst)/len(lst)

    ## Make final output
    out = df[["long_ID","Event_type",geneid_type,"E1","E2","E3","E4","strand"]]
    
    if header > 0:
        out.to_csv(outdir+"novelJC.csv",
                   sep=",",
                   mode="a",
                   index=None,
                   header=None)
    else:
        out.columns = ["long_ID","Event_type",args.geneid[0],"E1","E2","E3","E4","strand"]
        out.to_csv(outdir+"novelJC.csv",
                   sep=",",
                   index=None)
# %%
