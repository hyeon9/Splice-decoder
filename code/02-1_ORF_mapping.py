#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import matplotlib
import subprocess
import os
import sys
import time
import argparse
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
start_time = time.time()

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
 
 SSSSS   PPPP    L       III  CCC  EEEEE  DDDDD   EEEEE  CCCCC  OOO  DDDDD  EEEEE  RRRRR
S        P   P   L        I  C   C E      D   D  E      C     O   O D   D  E      R    R
 SSSS    PPPP    L        I  C     EEEE   D   D  EEEE   C     O   O D   D  EEEE   RRRRRR
    S    P       L        I  C   C E      D   D  E      C     O   O D   D  E      R   R
 SSSSS   P       LLLLL  III   CCC  EEEEE  DDDDD   EEEEE  CCCCC  OOO  DDDDD  EEEEE  R    R

Description
    #########################################################
    This script makes simulated transcript (Sim_TX) by using 
    given splicing events on the Ref_TX.
    OUTPUT: merged_{splicing_event}_w_event.bed, "{splicing_event}_bestorf_{splitted}.tsv"
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat file', 
                        required=True,
                        type=str)
    parser.add_argument("--cpat", "-cp", help="cpat path", 
                        required=True,
                        type=str)
    parser.add_argument("--cpatdb", "-cpdb", help="cpat db path", 
                        required=True,
                        type=str)
    parser.add_argument("--bedtools", "-b", help="bedtools path", 
                        required=True,
                        type=str)
    parser.add_argument("--genomefa", "-f", help="reference fasta", 
                        required=True,
                        type=str)
    parser.add_argument("--species", "-p", help="species name", 
                        required=True,
                        type=str)
    
    ## Optional
    parser.add_argument("--threads", "-t", help="number of threads to use, default: 9", 
                        type=int,
                        default="9")

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    
# %%
def Estimate(id, final_set, key):
    """ Estimate consequence of splicing change

    Args:
        id (_type_): event ID
        final_set (_type_): pair dataframe
        key (_type_): event type with multiple ID
        if PSI<0 focusing on EI, RI, A5SS, and A3SS >> What about MXE?
    """

    sim_bed = args.input+"sim_bed/"
    fasta_DIR = args.input+"get_Fasta/"
    cpat_DIR = args.input+"cpat/"

#    if not os.path.exists(sim_bed):
#        os.mkdir(sim_bed)

    temp = final_set[final_set["ID"]==id]
    for tx in temp["best match tx"].unique():   # For each matched tx
        mxe_ref = []
        start_ex = int(temp[temp["best match tx"]==tx]["start_exon"].str.replace("exon","").values[0])
        form = temp[temp["best match tx"]==tx]["event"].values[0]
        sub_gtf = exon_dict[exon_dict["ENSTID"]==(tx)]    # Genomic position of start matching exon
        sub_gtf = sub_gtf[sub_gtf["ENSTID"]==tx].sort_values(by="exon_number")
        ups = sub_gtf[sub_gtf["exon_number"].astype(int) <= start_ex][["chr","start","end","gene_symbol","ENSTID","strand"]]  # upstream exons (from matching exon)
        dos = sub_gtf[sub_gtf["exon_number"].astype(int) > start_ex][["chr","start","end","gene_symbol","ENSTID","strand"]]   # downstream exons
        
        ## Make simulated forms
        if not key.startswith("MXE"):
            index_num = temp[temp["best match tx"]==tx]["ID"].index[0]  # Matched TX index number in input
            insert = temp[temp["best match tx"]==tx]["ID"].str.split(";")   # 4-5,6-7,8-9,10-11
    
            if form == "ES":  # Skipped exon is ref > simulate INCLUDED 2nd exon case
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][6],insert[index_num][7],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                else:
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][7],insert[index_num][6],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                dos = pd.concat([new_line.T, dos])
                output = pd.concat([ups, dos])
                sim_form = "EI"   # simulated event
            
            elif form == "SI":  # Skipped intron is ref > simulate INCLUDED intron
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],ups.iloc[-1,1],dos.iloc[0,2],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"]) # Make retained intron
                else:
                    new_line = pd.DataFrame([insert[index_num][2],dos.iloc[0,1],ups.iloc[-1,2],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"]) # Make retained intron
                ups = ups.iloc[:-1,:]   # Remove first exon
                dos = dos.iloc[1:,:]    # Remove third (second) exon
                
                dos = pd.concat([new_line.T, dos])  # Add retained form
                output = pd.concat([ups, dos])
                sim_form = "RI"   # simulated event

            elif form == "Original_A5SS":  # O5SS is ref > simulate first exon to A5SS 
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],ups.iloc[-1,1],insert[index_num][7],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                else:
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][7],ups.iloc[-1,2],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                ups = ups.iloc[:-1,:]   # Remove frist exon
                ups = pd.concat([ups, new_line.T])  # Add first exon with A5SS
                output = pd.concat([ups, dos])
                sim_form = "Alt_A5SS"   # simulated event

            elif form == "Original_A3SS":  # O3SS is ref > simulate third exon to A3SS
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][6],dos.iloc[0,2],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                else:
                    new_line = pd.DataFrame([insert[index_num][2],dos.iloc[0,1],insert[index_num][6],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                dos = dos.iloc[1:,:]    # Remove second exon
                dos = pd.concat([new_line.T, dos])  # Add second exon with A3SS
                output = pd.concat([ups, dos])
                sim_form = "Alt_A3SS"   # simulated event

            elif form == "EI":  # Exon inclusion is ref > simulated Skipped 2nd exon case
                dos = dos.iloc[1:,:]    # remove second exon
                output = pd.concat([ups, dos])
                sim_form = "ES"
            
            elif form == "RI":  # When RI form is ref > simulate excluded intron case
                s = ups.iloc[-1,1]
                e = ups.iloc[-1,2]
                t1 = insert[index_num][5]
                t2 = insert[index_num][8]
                if insert[index_num][3] == "+":
                    new_line1 = pd.DataFrame([insert[index_num][2],s,t1,insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                    new_line2 = pd.DataFrame([insert[index_num][2],t2,e,insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                else:
                    new_line1 = pd.DataFrame([insert[index_num][2],t1,e,insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                    new_line2 = pd.DataFrame([insert[index_num][2],s,t2,insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])        
                ups = ups.iloc[:-1,:]    # Remove RI case
                ups = pd.concat([ups, new_line1.T])  # Insert E1
                ups = pd.concat([ups, new_line2.T])  # Insert E3
                output = pd.concat([ups, dos])
                sim_form = "SI" # simulated event

            elif form == "Alt_A5SS":  # A5SS is ref > simulate first exon to O5SS 
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],ups.iloc[-1,1],insert[index_num][5],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                else:
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][5],ups.iloc[-1,2],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"])
                ups = ups.iloc[:-1,:]   # Remove frist exon
                ups = pd.concat([ups, new_line.T])  # Add first exon with A5SS
                output = pd.concat([ups, dos])
                sim_form = "Ori_A5SS"   # simulated event

            elif form == "Alt_A3SS":  # A3SS is ref > simulate first exon to O3SS 
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][8],dos.iloc[0,2],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"]) # Make retained intron
                else:
                    new_line = pd.DataFrame([insert[index_num][2],dos.iloc[0,1],insert[index_num][8],insert[index_num][1],tx,insert[index_num][3]],
                                            index=["chr","start","end","gene_symbol","ENSTID","strand"],
                                            columns=["INSERT"]) # Make retained intron
                dos = dos.iloc[1:,:]    # Remove second exon
                dos = pd.concat([new_line.T, dos])  # Add second exon with A3SS
                output = pd.concat([ups, dos])
                sim_form = "Ori_A3SS"   # simulated event

        else:   # MXE cases
            insert = temp[temp["best match tx"]==tx]["ID"].str.split(";")
            index_num = temp[temp["best match tx"]==tx]["ID"].index[0]
          
            if form == "MXE1":   # E1,E2,E4
                dos = dos.iloc[1:,:]    # Delete original 2nd exon
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][8],insert[index_num][9],insert[index_num][3],tx],
                                            index=["chr","start","end","strand","ENSTID"],
                                            columns=["INSERT"])
                else:
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][9],insert[index_num][8],insert[index_num][3],tx],
                                            index=["chr","start","end","strand","ENSTID"],
                                             columns=["INSERT"])
                sim_form = "MXE2"
                dos = pd.concat([new_line.T, dos])  # Insert counterpart exon (From other MXE)
                output = pd.concat([ups, dos])

            elif form == "MXE2":   # E1,E3,E4
                dos = dos.iloc[1:,:]    # Delete original 2nd exon
                if insert[index_num][3] == "+":
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][6],insert[index_num][7],insert[index_num][3],tx],
                                            index=["chr","start","end","strand","ENSTID"],
                                            columns=["INSERT"])
                else:
                    new_line = pd.DataFrame([insert[index_num][2],insert[index_num][7],insert[index_num][6],insert[index_num][3],tx],
                                            index=["chr","start","end","strand","ENSTID"],
                                            columns=["INSERT"])
                sim_form = "MXE1"
                dos = pd.concat([new_line.T, dos])  # Insert counterpart exon (From other MXE)
                output = pd.concat([ups, dos])


        if len(mxe_ref) == 0:
            # sub_gtf = sub_gtf[["chr","start","end","ENSGID","gene_symbol","strand","ENSTID","exon_number"]]
            sub_gtf = sub_gtf[["chr","start","end","gene_symbol","ENSTID","strand"]]
            sub_gtf.loc[:,"start"] = sub_gtf.loc[:,"start"]-1   # For getfasta it start at start+1 position
        
        if sim_form != temp[temp["best match tx"]==tx]["event"].values[0]:  # It means Not interest event
            output.loc[:,"start"] = output.loc[:,"start"].astype(int)-1   # 1-based bed to 0-based bed
            if len(mxe_ref) == 0:
                output["sim_event"] = sim_form
                sub_gtf["sim_event"] = sim_form
                output["SID"] = id
                sub_gtf["SID"] = id
                output.to_csv(sim_bed+"{}_w_event.bed".format(key), sep="\t",
                            index=None, header=None)  # With any kinds of event (differ with gtf)
                sub_gtf.to_csv(sim_bed+"{}_wo_event.bed".format(key), sep="\t",
                            index=None, header=None)  # Origianl gtf
                output.to_csv(sim_bed+"merged_{}_w_event.bed".format(key), sep="\t",
                            index=None, header=None, mode="a")  # With any kinds of event (differ with gtf)
                sub_gtf.to_csv(sim_bed+"merged_{}_wo_event.bed".format(key), sep="\t",
                            index=None, header=None, mode="a")  # Origianl gtf

            else:   # for 4 exons match MXE
                output["sim_event"] = sim_form
                mxe_ref["sim_event"] = sim_form
                output["SID"] = id
                mxe_ref["SID"] = id
                mxe_ref.loc[:,"start"] = mxe_ref.loc[:,"start"].astype(int)-1   # 1-based bed to 0-based bed
                output.to_csv(sim_bed+"{}_w_event.bed".format(key), sep="\t",
                            index=None, header=None)  # With any kinds of event (differ with gtf)
                mxe_ref.to_csv(sim_bed+"{}_wo_event.bed".format(key), sep="\t",
                            index=None, header=None)  # Origianl gtf
                output.to_csv(sim_bed+"merged_{}_w_event.bed".format(key), sep="\t",
                            index=None, header=None, mode="a")  # With any kinds of event (differ with gtf)
                mxe_ref.to_csv(sim_bed+"merged_{}_wo_event.bed".format(key), sep="\t",
                            index=None, header=None, mode="a")  # Origianl gtf
                                    

        def Get_seq(input, key):
            """ Extract NT sequence based on gotten bed file (considered strand)

            Args:
                input (_type_): _description_
                key (_type_): _description_
            """
#            if not os.path.exists(fasta_DIR):
#                os.mkdir(fasta_DIR)

            ## Make AA sequence
            cmd = args.bedtools + " getfasta -fi " + args.genomefa + \
                " -bed " + sim_bed + key + "_"+input+".bed" + \
                " -fo " + fasta_DIR + key + "_"+input+".fa" + \
                " -s > /dev/null 2>>" + fasta_DIR +"/getfasta.log"
            subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

            files = open(fasta_DIR + key + "_"+input+".fa", "r")
            cpat_input = open(fasta_DIR + key + "_"+input+"_cpat.fa", "w")
            merge_seq = ""
            for line in files:
                if not line.startswith(">"):
                    merge_seq = merge_seq + line[:-1]
            
            proc_merged_seq = merge_seq.upper()
            cpat_input.write(">"+tx+"_"+input)
            cpat_input.write("\n"+proc_merged_seq)
            cpat_input.close()
                    

        def ORF(input, key):
#            if not os.path.exists(cpat_DIR):
#                os.mkdir(cpat_DIR)

            if input == "wo_event":
                n_orf = 10
            else:
                n_orf = 10
                
            ## Find ORF
            cmd = args.cpat + " -x " + args.cpatdb + args.species + "_Hexamer.tsv" + \
                    " -d " + args.cpatdb + args.species + "*_logit*RData" + \
                    " --min-orf 1" + " --top-orf={}".format(n_orf) + \
                    " -g "+ fasta_DIR+key+"_"+input+"_cpat.fa" + \
                    " -o" + cpat_DIR+input+"_"+key+".orf > /dev/null 2>>" + cpat_DIR +"/cpat.log"
            subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    
            try:    # Best ORF does not exist, fill zeros
                no_orf = pd.read_csv(cpat_DIR+"{}_{}.orf.no_ORF.txt".format(input,key),
                                     sep="\t")  # If no_ORF exists, output will be filled by zeros
                cpat_output = open(cpat_DIR+"{}_bestorf_{}.tsv".format(input,key), "a")
                cpat_output.write(id+"\t"+
                                  tx+"_"+input+"\t"+
                                  form+"\t"+
                                  str(0)+"\t"+
                                  str(0)+"\t"+
                                  str(0)+"\t"+
                                  str(0)+"\t"+
                                  str(0)+"\n")   # Add matched tx, nt length of frame, coding prob, eventid
                cpat_output.close()
                cmd = "rm " + cpat_DIR+"{}_{}.orf.no_ORF.txt".format(input,key)
                subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)	# To avoid re-use this file in next loop


            except: # Best ORF exists
                files = open(cpat_DIR+"{}_{}.orf.ORF_prob.tsv".format(input,key), "r")
                cpat_output = open(cpat_DIR+"{}_bestorf_{}.tsv".format(input,key), "a")
                for line in files:
                    if line.split("\t")[9].startswith("Coding_prob"):
                        pass
                    else:
                        cpat_output.write(id+"\t"+
                                          tx+"_"+input+"\t"+
                                          sim_form+"\t"+
                                          line.split("\t")[1]+"\t"+
                                          line.split("\t")[4]+"\t"+
                                          line.split("\t")[5]+"\t"+
                                          line.split("\t")[6]+"\t"+
                                          line.split("\t")[9])   # 1: tx length, 4,5: ORF start/end, 6: ORF length, 9: P coding
                cpat_output.close()
                cmd = "rm " + cpat_DIR+"{}_{}.orf.ORF_prob.tsv".format(input,key)
                subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)	# To avoid re-use this file in next loop

        if sim_form != temp[temp["best match tx"]==tx]["event"].values[0]:  # It means some of events do not have a RefTX for certain possible case
            Get_seq("w_event", key)
            Get_seq("wo_event", key)
            ORF("w_event", key)
            ORF("wo_event", key)


def Run(final_set, key):
    """ Run get Fasta & ORF prediction

    Args:
        final_set (_type_): query dataset (mapped DS events)
        key (_type_): n_th dataframe
    """
#    try:    # Before make data, remove every old data
#        os.remove(args.input+"sim_bed/merged_{}_w_event.bed".format(key)) # With any kinds of event (differ with gtf)
#        os.remove(args.input+"sim_bed/merged_{}_wo_event.bed".format(key)) # With any kinds of event (differ with gtf)

#    except:
#        pass

    for id in final_set["ID"].unique():
        Estimate(id, final_set, key)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    args.species = args.species.capitalize()
    
    if not args.input.endswith("/"):
        args.input = args.input+"/"
    
    ## Run code
    exon_dict = pd.read_csv(args.input+"exon_only.final.gtf",
                        sep="\t")
    query = pd.read_csv(args.input+"fig_input.txt",
                        sep="\t")
    rmat = pd.read_csv(args.input+"rmat.csv",
                    sep=",")

    import time, os
    from joblib import Parallel, delayed

    for kws in ["CA","RI","A5SS","A3SS","MXE"]:
        df_split = np.array_split(query[query["ID"].str.startswith(kws)],args.threads) # Split into n dataframe
        result = Parallel(n_jobs=args.threads)(delayed(Run)(sub_df, "{}{}".format(kws,k)) for k,sub_df in enumerate(df_split))
