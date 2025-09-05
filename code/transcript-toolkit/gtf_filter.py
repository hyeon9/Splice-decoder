# %%
import pandas as pd
import numpy as np
import matplotlib
import subprocess
import argparse
import os
import time
import sys

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
start_time = time.time()

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''

Description
    #########################################################
    This script makes simulated transcript (Sim_TX) by using 
    given splicing events on the Ref_TX.
    #########################################################''',

    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat file', 
                        required=True,
                        type=str)
    
    parser.add_argument("--path", "-p", help="Gencode canonical path", 
                        required=True,
                        type=str)
    
    parser.add_argument("--gene_list", "-g", help="Target gene list with full path", 
                        required=True,
                        type=str)
    
    ## Optional
    parser.add_argument("--canonical", "-c", help="Way to define canonical transcripts (default: using Gencode based list)", 
                        type=str,
                        default="default")
    
    
    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args


args, parser = parse_args(sys.argv[1:])


if args.canonical == "default":
    print("Using Gencode canonical")
    gene_list = pd.read_csv(f"{args.gene_list}",
                            sep="\t",
                            header=None)
    
    gencode = pd.read_csv(f"{args.path}/dat/hg38_cano.txt",
                          sep="\t",
                          header=None)
    
    gencode = gencode[gencode[0].isin(gene_list[0])]
    gene_list = gencode.copy()
    
else:
    print("Using user-defined canonical")
    gene_list = pd.read_csv(f"{args.gene_list}",
                            sep="\t",
                            header=None)
    
cano = gene_list[1].tolist()
query_list = []
with open(f"{args.input}main.gtf", 'r') as whole_gtf:
    filtered_gtf = open(f"{args.input}exon_only.gtf", 'w')
    for line in whole_gtf:
        lines = line.strip().split("\t")
        for k,key in enumerate(lines[-1].split("\"")):
            if "gene_id" in key:    # ENSG
                gene = lines[-1].split("\"")[k+1]
            elif "transcript_id" in key:
                tx_id = lines[-1].split("\"")[k+1]
        if gene in gene_list[0].tolist() and lines[2] == "exon":
            filtered_gtf.write(line)
            query_list.append([gene,tx_id])

    filtered_gtf.close()

query_list_df = pd.DataFrame(query_list)
query_list_df = query_list_df.drop_duplicates()
query_table = pd.merge(gene_list, query_list_df,
                       left_on=[0], right_on=[0])
query_table_df = query_table[['1_x','1_y']]
query_table_df.columns = ["Major","query"]
query_table_df = query_table_df.drop_duplicates()
query_table_df.to_csv(f'{args.input}query_list.txt',
                      sep="\t",
                      index=None)

# %%
