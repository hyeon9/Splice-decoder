#!/usr/bin/env python3
# %%
import pandas as pd
import pybedtools
import os
import argparse
import sys

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(
    description='''
    Description
    #########################################################
    Raw gencode gtf should be processed using below command before run this script
    "awk -F "\t" '{if ($3 == "exon") print }' ${Your.gtf} > exon_only.gtf"
    This script makes proper format of exon gtf
    #########################################################
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Your working directory (input)', 
                        required=True,
                        type=str)
    parser.add_argument('--ref', '-r', 
                        help='reference genome fasta', 
                        required=True,
                        type=str)
    parser.add_argument('--tx', '-t', 
                        help='query transcript', 
                        required=True,
                        type=str)
    parser.add_argument('--event', '-e', 
                        help='simulation event type', 
                        required=True,
                        type=str)
    parser.add_argument('--dat', '-d',
                        help='SpliceDeocder install point',
                        required=True,
                        type=str)

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

def Make_AA():
    ## Make output dir
    out_dir = f"{args.input}AF2/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ## Load output_table
    table = pd.read_csv(f"{args.input}result/Effect_score.tsv",
                        sep="\t")
    sub_table = table[(table["Reference_transcript"]==args.tx) & \
                    (table["Simulated_event"]==args.event)]
    splicing_type = sub_table["LongID"].unique()[0].split(";")[0]
    ref_start = int(float(sub_table["AUG (Ref-Sim)"].unique()[0].split("-")[0]))
    ref_end = int(float(sub_table["Stop_codon (Ref-Sim)"].unique()[0].split("-")[0]))
    sim_start = int(float(sub_table["AUG (Ref-Sim)"].unique()[0].split("-")[1]))
    sim_end = int(float(sub_table["Stop_codon (Ref-Sim)"].unique()[0].split("-")[1]))

    def Make_bed(in_bed):
        ## Load Sim bed
        merged_bed = pd.read_csv(f"{args.input}post_input/{in_bed}",
                                sep="\t",
                                header=None)
        bed = merged_bed[merged_bed[4]==args.tx]
        gene_symbol = bed[3].dropna().unique()[0]
        bed[3] = gene_symbol
        bed_str = bed.to_csv(sep="\t", index=False, header=False)
        bed = pybedtools.BedTool(bed_str, from_string=True)
        
        return bed

    ref_bed = Make_bed(f"merged_wo_{splicing_type}.bed")
    sim_bed = Make_bed(f"merged_w_{splicing_type}.bed")

    ## Specify ref
    if args.ref == "human" or args.ref == "Human":
        ref = f"{args.dat}/dat/reference/hg38.fa"
    else:
        ref = f"{args.dat}/dat/reference/GRCm38.fa"

    ## Extract sequence
    def Get_Fasta(bed,d_type):
        fasta_sequences = bed.sequence(fi=ref, fo=f"{out_dir}{d_type}.fa", s=True)
        
        return fasta_sequences

    Get_Fasta(ref_bed, "ref")
    Get_Fasta(sim_bed, "sim")

    ## CDS specify
    from Bio import SeqIO
    ref_fa = "".join(str(record.seq) for record in SeqIO.parse(f"{out_dir}ref.fa", "fasta"))
    sim_fa = "".join(str(record.seq) for record in SeqIO.parse(f"{out_dir}sim.fa", "fasta"))
    ref_cds = ref_fa[ref_start-1:ref_end]
    sim_cds = sim_fa[sim_start-1:sim_end]

    ## Translate and pairwise align
    from Bio.Seq import Seq
    ref_aa = (Seq(ref_cds).translate())
    sim_aa = (Seq(sim_cds).translate())

    print("Copy and Paste these sequence at alphafold2 web page () to construct your 3D structure")
    print(f"Ref_tx: {ref_aa}","\n")
    print(f"Sim_tx: {sim_aa}")
    # %%
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    else:
        args, parser = parse_args(sys.argv[1:])    
        Make_AA()
