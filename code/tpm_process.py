#!/usr/bin/env python3
# %%
import sys
import argparse
def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(
        description='''
 
 SSSSS   PPPP    L       III  CCC  EEEEE  DDDDD   EEEEE  CCCCC  OOO  DDDDD  EEEEE  RRRRR
S        P   P   L        I  C   C E      D   D  E      C     O   O D   D  E      R    R
 SSSS    PPPP    L        I  C     EEEE   D   D  EEEE   C     O   O D   D  EEEE   RRRRRR
    S    P       L        I  C   C E      D   D  E      C     O   O D   D  E      R   R
 SSSSS   P       LLLLL  III   CCC  EEEEE  DDDDD   EEEEE  CCCCC  OOO  DDDDD  EEEEE  R    R

Description
    #########################################################
    Make merged TPM matrix using individual TPM file from Stringtie
    stringtie -e -B -p 45 -G ${gtf} -o ${id}.tab ${bam}
    #########################################################
    ''',
        formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains TPM file', 
                        required=True,
                        type=str)
    parser.add_argument('--file', '-n', 
                        help='TPM file of individual', 
                        required=True,
                        type=str)

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

args, parser = parse_args(sys.argv[1:])
args.input = args.input+"/"
name = args.file.replace("_","-")
name = "-".join(name.split("-")[:3])

data = open(args.input+args.file, "r")
out = open(args.input+name+".tpm", "w")
out.write("ID"+"\t"+name+"\n")
for line in data:
    readline = line.strip().split("\t")
    if not readline[0].startswith("#") and \
    readline[2] == "transcript":
        if "transcript_id" in readline[8].split(" "):
            tx = readline[8].split(" ")[readline[8].split(" ").index("transcript_id")+1].split("\"")[1]
        if "TPM" in readline[8].split(" "):
            tpm = readline[8].split(" ")[readline[8].split(" ").index("TPM")+1].split("\"")[1]
        else:
            gene = ""
        out.write(tx+"\t"+tpm+"\n")
# %%
