#!/usr/bin/env python3
# %%
import argparse
import sys

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(
        description='''
        
Description
    #########################################################
    This script adds exon number tag in your gtf.
    If your gtf file does not have exon number, Splice-decoder does not work properly
    #########################################################
    ''',
        formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains exon gtf file', 
                        required=True,
                        type=str)

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

def Make_gene_dict(DIR):
    import os
    files = [f for f in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, f)) if f.startswith("main") if f.endswith(".gtf") ]
    if len(files) == 0:
        print("Check your main.gtf")
    else:
        f = open(DIR+"{}".format(files[0]))
        o = open(DIR+"tx_gene_dict", "w")
        for line in f:
            read_line = line.strip().split("\t")
            if not line.startswith("#"):
                if read_line[2] == "transcript":
                    if "transcript_id" in read_line[8].split(" "):
                        tx = read_line[8].split(" ")[read_line[8].split(" ").index("transcript_id")+1].split("\"")[1]
                    if "gene_id" in read_line[8].split(" "):
                        gene_id = read_line[8].split(" ")[read_line[8].split(" ").index("gene_id")+1].split("\"")[1]
                        if "gene_name" in read_line[8].split(" "):
                            gene_sym = read_line[8].split(" ")[read_line[8].split(" ").index("gene_name")+1].split("\"")[1]
                        else:
                            gene_sym = ""
                    
                    o.write(tx+"\t"+gene_id+"\t"+gene_sym+"\n")
        f.close()
        o.close()


args, parser = parse_args(sys.argv[1:])
Make_gene_dict(args.input)
