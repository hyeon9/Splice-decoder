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
    parser.add_argument('--FDR', '-fdr', nargs='*',
                        help='Filter out DS with higher FDR than cut off', 
                        type=str,
                        default=[0.05])
    parser.add_argument('--PSI', '-psi', nargs='*',
                        help='Filter out DS less PSI than cut off', 
                        type=str,
                        default=[0.1])
    parser.add_argument('--geneid', '-id', nargs='*',
                        help='Column name of gene in GTF', 
                        type=str,
                        default=['Gene'])

    args = parser.parse_args(cmd_args, namespace)

    return args, parse_args

args, parser = parse_args(sys.argv[1:])
dir = args.rMATS_folder
outdir = args.input
FDR_cut = float(args.FDR[0])
PSI_cut = float(args.PSI[0])



df = pd.read_csv(f"{dir}val_PE.txt",
                sep="\t")
geneid_type = "gene_name"
event_type = "SE"   # They only considered alternative exon events
df[[geneid_type]] = df[[geneid_type]].astype(str)   # BUG fixed 241120

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
    neg_df["E4s"] = None
    neg_df["E4e"] = None

    df = pd.concat([pos_df, neg_df])
    df["Event_type"] = "CA"

df["E1"] = df["chr"]+":"+df["E1s"].astype(str)+"-"+df["E1e"].astype(str)
df["E2"] = df["chr"]+":"+df["E2s"].astype(str)+"-"+df["E2e"].astype(str)
df["E3"] = df["chr"]+":"+df["E3s"].astype(str)+"-"+df["E3e"].astype(str)
if event_type == "MXE":
    df["E4"] = df["chr"]+":"+df["E4s"].astype(str)+"-"+df["E4e"].astype(str)
else:
    df["E4"] = None

out = df[["long_ID","Event_type",geneid_type,"E1","E2","E3","E4","strand","NMD_isoforms","gene_name","gene_id"]]
out.columns = ["long_ID","Event_type",args.geneid[0],"E1","E2","E3","E4","strand","NMD_isoforms","gene_name","gene_id"]
out.to_csv(outdir+"PE.csv",
            sep="\t",
            index=None)
# %%
