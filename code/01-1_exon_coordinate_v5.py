#!/usr/bin/env python3
#%%
""" 
Make matched transcript list with the given splicing events

Returns:
    _type_: Input of the 01-2_stat_exon_coor_result.py
"""
import pandas as pd
import numpy as np
import argparse
import sys

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
Description
    #########################################################
    This script finds matched transcript (Ref_TX) to simulate splicing 
    events.
    OUTPUT: short_output.txt, long_output.txt
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains exon gtf file and rmat file', 
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


# %%
## Split start, end in foward and reverse genes respectively
def Split(df):
    """ Split genomic position into the start and end

    Args:
        df (_type_): original dataframe

    Returns:
        _type_: start position (end position on reverse)
    """
    try:
        st = df[0].split(":")[1].split("-")[0]
    except:
        return df[0]
    return st

def Split2(df):
    """ Split genomic position into the start and end

    Args:
        df (_type_): original dataframe

    Returns:
        _type_: end position (start position on reverse)
    """
    try:
        ed = df[0].split(":")[1].split("-")[1]
    except:
        return df[0]
    return ed

# %%
def Exon_comp(temp, gene, event_type):
    """ Find best exon structure

    Args:
        temp (_type_): Processed rMATS dataframe
        gene (_type_): gene symbol
        event_type (_type_): Queried splicing event type
    """

    def Comp_with_2exon(eventtype, ts1, te1, *ext):
        # Special case, some annotated isoform has only 2 exon,
        # consider their 2nd exon as third, assume only exon skipping situation
        # ext is additional exon

        best_mat = [1e6]    # default matching value
        best_info = []
        if len(ext) > 0:
            if eventtype == "ES" or eventtype == "SI":
                first = (np.abs(tx_annot.iloc[0,2]-int(temp[te1])))
                second = "SKIP"
                third = np.abs(tx_annot.iloc[1,1]-int(temp[ext[0]]))
                fourth = "NA"
                diff_sum = int(first+third)
            
            elif eventtype.endswith("A3SS"):
                if eventtype.startswith("Original"):
                    first = np.abs(tx_annot.iloc[0,2]-int(temp[te1]))
                    second = "SKIP"
                    third = np.abs(tx_annot.iloc[1,1]-int(temp[ext[0]]))
                    fourth = "NA"
                else:
                    first = np.abs(tx_annot.iloc[0,2]-int(temp[te1]))
                    second = "SKIP"
                    if temp["strand"] == "+":
                        if np.abs(tx_annot.iloc[1,1]-int(temp[ext[0]])) == 0 and \
                            tx_annot.iloc[1,2] > int(temp[ext[1]]):
                            third = 0
                        else:
                            third = 1e6
                    else:
                        if np.abs(tx_annot.iloc[1,1]-int(temp[ext[0]])) == 0 and \
                            tx_annot.iloc[1,2] < int(temp[ext[1]]):
                            third = 0
                        else:
                            third = 1e6
                    fourth = "NA"
                diff_sum = int(first+third)
            
            elif eventtype.endswith("A5SS"):
                if eventtype.startswith("Original"):
                    first = np.abs(tx_annot.iloc[0,2]-int(temp[te1]))
                    second = np.abs(tx_annot.iloc[1,1]-int(temp[ext[0]]))
                    third = "SKIP"
                    fourth = "NA"
                else:
                    if temp["strand"] == "+":
                        if np.abs(tx_annot.iloc[0,2]-int(temp[te1])) == 0 and \
                            tx_annot.iloc[0,1] < int(temp[ts1]):
                            first = 0
                        else:
                            first = 1e6
                    else:
                        if np.abs(tx_annot.iloc[0,2]-int(temp[te1])) == 0 and \
                            tx_annot.iloc[0,1] > int(temp[ts1]):
                            first = 0
                        else:
                            first = 1e6
                    second = np.abs(tx_annot.iloc[1,1]-int(temp[ext[0]]))
                    third = "SKIP"
                    fourth = "NA"
                diff_sum = int(first+second)

        else:   # Intron retention case
            if temp["strand"] == "+":
                if tx_annot.iloc[0,1] <= int(temp[ts1]) and \
                tx_annot.iloc[0,2] >= int(temp[te1]):
                    first = 0
                else:
                    first = 1e6
            else:
                if tx_annot.iloc[0,1] >= int(temp[ts1]) and \
                    tx_annot.iloc[0,2] <= int(temp[te1]):
                    first = 0
                else:
                    first = 1e6
            second = "RI"
            third = "RI"
            fourth = "NA"
            diff_sum = int(first)

        if int(diff_sum) < min(best_mat):   # Report best match form
            best_mat = []
            best_info = []
            best_mat.append(diff_sum)
            best_info.append(tx)
            best_info.append("exon1")    # becasue exon start with 1
            best_info.append(first)
            best_info.append(second)
            best_info.append(third)
            best_info.append(fourth)

        if len(best_mat) > 0 and len(best_info) > 0:
            if len(ext) > 0:
                short_form.write(str(temp["long_ID"])+"\t"+temp[geneID]+"\t"+str(temp[ts1])+"\t"+str(temp[te1])+"\t"+
                                str("SKIP")+"\t"+str("SKIP")+"\t"+str(temp[ext[0]])+"\t"+str(temp[ext[1]])+"\t"+
                                str("NA")+"\t"+str("NA")+"\t"+
                                str(best_mat[0])+"\t"+str(best_info[0])+"\t"+str(best_info[1])+"\t"+str(best_info[2])+"\t"+
                                str(best_info[3])+"\t"+str(best_info[4])+"\t"+str(best_info[5])+"\t"+eventtype+"\n")
            
            else:
                short_form.write(str(temp["long_ID"])+"\t"+temp[geneID]+"\t"+str(temp[ts1])+"\t"+str(temp[te1])+"\t"+
                                str("RI")+"\t"+str("RI")+"\t"+str("RI")+"\t"+str("RI")+"\t"+
                                str("NA")+"\t"+str("NA")+"\t"+
                                str(best_mat[0])+"\t"+str(best_info[0])+"\t"+str(best_info[1])+"\t"+str(best_info[2])+"\t"+
                                str(best_info[3])+"\t"+str(best_info[4])+"\t"+str(best_info[5])+"\t"+eventtype+"\n")


    def Comp_with_multiexon(eventtype, ts1, te1, *ext):
        best_mat = [1e6]    # default matching value
        best_info = []
        try:    # According to the input size, the end point will be set
            end_point = len(ext)/2
            if eventtype == "MXE":
                end_point -= 1  # UPDATE 25.05.07 to avoid early stop in tiling
        except: # For RI, because ext = 0
            end_point = 0
        for case in range(tx_annot.shape[0]):   # Count exon, case is the exon number
            if case+end_point <= (tx_annot.shape[0]-1):
                if len(ext) >= 6:   # MXE cases
                    first = np.abs(tx_annot.iloc[case,2]-int(temp[te1]))
                    second = (np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]])) +
                                np.abs(tx_annot.iloc[case+1,2]-int(temp[ext[1]])))
                    third = "SKIP"
                    fourth = np.abs(tx_annot.iloc[case+2,1]-int(temp[ext[4]]))
                    diff_sum = int(first+second+fourth)
                    eventtype = "MXE1"
                    
                    if diff_sum != 0:
                        first = np.abs(tx_annot.iloc[case,2]-int(temp[te1]))
                        second = (np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[2]])) +
                                    np.abs(tx_annot.iloc[case+1,2]-int(temp[ext[3]])))
                        third = "SKIP"
                        fourth = np.abs(tx_annot.iloc[case+2,1]-int(temp[ext[4]]))
                        diff_sum = int(first+second+fourth)
                        eventtype = "MXE2"
                        
                elif len(ext) >= 4 and len(ext) < 6:   # Exon inclusion
                    first = np.abs(tx_annot.iloc[case,2]-int(temp[te1]))
                    second = (np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]])) +
                              np.abs(tx_annot.iloc[case+1,2]-int(temp[ext[1]])))
                    third = np.abs(tx_annot.iloc[case+2,1]-int(temp[ext[2]]))
                    fourth = "NA"
                    diff_sum = int(first+second+third)
                
                elif len(ext) >= 2 and len(ext) < 4:   # Exon skipped and SI
                    if eventtype == "ES" or eventtype == "SI":
                        first = np.abs(tx_annot.iloc[case,2]-int(temp[te1]))
                        second = "SKIP"
                        third = np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]]))
                        fourth = "NA"
                        diff_sum = int(first+third)
                    
                    elif eventtype.endswith("A3SS"):
                        if eventtype.startswith("Original"):
                            first = np.abs(tx_annot.iloc[case,2]-int(temp[te1]))
                            second = "SKIP"
                            third = np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]]))
                            fourth = "NA"
                        else:
                            first = np.abs(tx_annot.iloc[case,2]-int(temp[te1]))
                            second = "SKIP"
                            if temp["strand"] == "+":
                                if np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]])) == 0 and \
                                    tx_annot.iloc[case+1,2] > int(temp[ext[1]]):
                                    third = 0
                                else:
                                    third = 1e6
                            else:
                                if np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]])) == 0 and \
                                    tx_annot.iloc[case+1,2] < int(temp[ext[1]]):
                                    third = 0
                                else:
                                    third = 1e6
                            fourth = "NA"
                        diff_sum = int(first+third)
                    
                    elif eventtype.endswith("A5SS"):
                        if eventtype.startswith("Original"):
                            first = np.abs(tx_annot.iloc[case,2]-int(temp[te1]))
                            second = np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]]))
                            third = "SKIP"
                            fourth = "NA"
                        else:
                            if temp["strand"] == "+":
                                if np.abs(tx_annot.iloc[case,2]-int(temp[te1])) == 0 and \
                                    tx_annot.iloc[case,1] < int(temp[ts1]):
                                    first = 0
                                else:
                                    first = 1e6
                            else:
                                if np.abs(tx_annot.iloc[case,2]-int(temp[te1])) == 0 and \
                                    tx_annot.iloc[case,1] > int(temp[ts1]):
                                    first = 0
                                else:
                                    first = 1e6
                            second = np.abs(tx_annot.iloc[case+1,1]-int(temp[ext[0]]))
                            third = "SKIP"
                            fourth = "NA"
                        diff_sum = int(first+second)
                
                else:   # Retained intron case, All exon that cover entire intron of retained intron block
                    if temp["strand"] == "+":
                        if tx_annot.iloc[case,1] <= int(temp[ts1]) and \
                        tx_annot.iloc[case,2] >= int(temp[te1]):
                            first = 0
                        else:
                            first = 1e6
                    else:
                        if tx_annot.iloc[case,1] >= int(temp[ts1]) and \
                        tx_annot.iloc[case,2] <= int(temp[te1]):
                            first = 0
                        else:
                            first = 1e6
                    second = "RI"
                    third = "RI"
                    fourth = "NA"
                    diff_sum = int(first) # In case of RI (skipped intron)
                    
                
                if int(diff_sum) < min(best_mat):
                    best_mat = []   # remove previous match info
                    best_info = []  # remove previous match info
                    best_mat.append(diff_sum)
                    best_info.append(tx)
                    best_info.append("exon{}".format(str(case+1)))    # becasue exon start with 1
                    best_info.append(first)
                    best_info.append(second)
                    best_info.append(third)
                    best_info.append(fourth)
                    best_info.append(eventtype)
            
        if len(best_mat) > 0 and len(best_info) > 0:
            if len(ext) >= 6:
                long_form.write(str(temp["long_ID"])+"\t"+temp[geneID]+"\t"+str(temp[ts1])+"\t"+str(temp[te1])+"\t"+
                                str(temp[ext[0]])+"\t"+str(temp[ext[1]])+"\t"+str(temp[ext[2]])+"\t"+str(temp[ext[3]])+"\t"+
                                str(temp[ext[4]])+"\t"+str(temp[ext[5]])+"\t"+
                                str(best_mat[0])+"\t"+str(best_info[0])+"\t"+str(best_info[1])+"\t"+str(best_info[2])+"\t"+
                                str(best_info[3])+"\t"+str(best_info[4])+"\t"+str(best_info[5])+"\t"+str(best_info[6])+"\n")

            elif len(ext) >= 4 and len(ext) < 6:
                long_form.write(str(temp["long_ID"])+"\t"+temp[geneID]+"\t"+str(temp[ts1])+"\t"+str(temp[te1])+"\t"+
                                str(temp[ext[0]])+"\t"+str(temp[ext[1]])+"\t"+str(temp[ext[2]])+"\t"+str(temp[ext[3]])+"\t"+
                                str("NA")+"\t"+str("NA")+"\t"+
                                str(best_mat[0])+"\t"+str(best_info[0])+"\t"+str(best_info[1])+"\t"+str(best_info[2])+"\t"+
                                str(best_info[3])+"\t"+str(best_info[4])+"\t"+str(best_info[5])+"\t"+str(best_info[6])+"\n")
    
            elif len(ext) >= 2 and len(ext) < 4:
                long_form.write(str(temp["long_ID"])+"\t"+temp[geneID]+"\t"+str(temp[ts1])+"\t"+str(temp[te1])+"\t"+
                                str("SKIP")+"\t"+str("SKIP")+"\t"+str(temp[ext[0]])+"\t"+str(temp[ext[1]])+"\t"+
                                str("NA")+"\t"+str("NA")+"\t"+
                                str(best_mat[0])+"\t"+str(best_info[0])+"\t"+str(best_info[1])+"\t"+str(best_info[2])+"\t"+
                                str(best_info[3])+"\t"+str(best_info[4])+"\t"+str(best_info[5])+"\t"+str(best_info[6])+"\n")
    
            else:
                long_form.write(str(temp["long_ID"])+"\t"+temp[geneID]+"\t"+str(temp[ts1])+"\t"+str(temp[te1])+"\t"+
                                str("RI")+"\t"+str("RI")+"\t"+str("RI")+"\t"+str("RI")+"\t"+
                                str("NA")+"\t"+str("NA")+"\t"+
                                str(best_mat[0])+"\t"+str(best_info[0])+"\t"+str(best_info[1])+"\t"+str(best_info[2])+"\t"+
                                str(best_info[3])+"\t"+str(best_info[4])+"\t"+str(best_info[5])+"\t"+str(best_info[6])+"\n")
                

    #annot = gtf[gtf["gene_symbol"]==gene]   # Consider matched gene
    # annot = gtf[(gtf["gene_symbol"]==gene) | (gtf["ENSGID"]==gene)]   # Consider matched gene
    annot = gtf[(gtf["ENSGID"]==gene)]   # Consider matched gene
    for tx in annot["ENSTID"].unique(): # Collect all isoform in each gene from comprehensive gtf
        tx_annot = annot[annot["ENSTID"]==tx].sort_values(by="exon_number")   # Consider one isoform at once, We MUST start from exon 1
        if tx_annot.shape[0]==2:  # A transcript has less than 2 exons
            if event_type == "CA":
                Comp_with_2exon("ES","E1s","E1e","E3s","E3e")
            elif event_type == "A3SS":
                Comp_with_2exon("Original_A3SS","E1s","E1e","E3s","E3e")
                Comp_with_2exon("Alt_A3SS","E1s","E1e","E2s","E3s")
            elif event_type == "A5SS":
                Comp_with_2exon("Original_A5SS","E1s","E1e","E3s","E3e")
                Comp_with_2exon("Alt_A5SS","E1e","E2e","E3s","E3e")
            elif event_type == "RI":
                Comp_with_2exon("RI","E1e","E3s")
                Comp_with_2exon("SI","E1s","E1e","E3s","E3e")
        
        else:
            if event_type == "CA":
                Comp_with_multiexon("EI","E1s","E1e","E2s","E2e","E3s","E3e")
                Comp_with_multiexon("ES","E1s","E1e","E3s","E3e")
            elif event_type == "A3SS":
                Comp_with_multiexon("Original_A3SS","E1s","E1e","E3s","E3e")
                Comp_with_multiexon("Alt_A3SS","E1s","E1e","E2s","E3s")
            elif event_type == "A5SS":
                Comp_with_multiexon("Original_A5SS","E1s","E1e","E3s","E3e")
                Comp_with_multiexon("Alt_A5SS","E1e","E2e","E3s","E3e")
            elif event_type == "RI":
                Comp_with_multiexon("RI","E1e","E3s")
                Comp_with_multiexon("SI","E1s","E1e","E3s","E3e")
            elif event_type == "MXE":
                Comp_with_multiexon("MXE","E1s","E1e","E2s","E2e","E3s","E3e","E4s","E4e")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/"):
        args.input = args.input+"/"
    rmat = pd.read_csv(args.input+"rmat.csv",
                    sep=",")
    print(f"Filter passed DS events: {rmat.shape[0]}")

    ## Run code
    gtf = pd.read_csv(args.input+"exon_only.final.gtf",
                  sep="\t")
    geneID = "Gene"
    
    if args.target != "all":
        target_Genes = pd.read_csv(f"{args.target}",
                                   sep="\t",
                                   header=None)
        rmat = rmat[rmat[geneID].isin(target_Genes[0])]
    
    else:
        pass

    print("Significant DS: {}".format(rmat.shape[0]))
    event_exon = rmat[["long_ID","Event_type",geneID,"E1","E2","E3","E4","strand"]] # input header
#    print(event_exon)
        
    rev_event_exon = event_exon[event_exon["strand"]=="-"]
    fow_event_exon = event_exon[event_exon["strand"]=="+"]

    for i in range(1,5):
        fow_event_exon["E{}s".format(i)] = fow_event_exon[["E{}".format(i)]].apply(Split,axis=1)
        fow_event_exon["E{}e".format(i)] = fow_event_exon[["E{}".format(i)]].apply(Split2,axis=1)
        
        rev_event_exon["E{}s".format(i)] = rev_event_exon[["E{}".format(i)]].apply(Split2,axis=1)
        rev_event_exon["E{}e".format(i)] = rev_event_exon[["E{}".format(i)]].apply(Split,axis=1)

    event_exon = pd.concat([fow_event_exon, rev_event_exon])    # Processed rMATS output

    ## Convert start and end of revers strand genes of gtf
    rev_annot = gtf[gtf["strand"]=="-"]
    fow_annot = gtf[gtf["strand"]!="-"]
    rev_annot.columns = ["chr","end","start","strand","ENSGID","gene_symbol","ENSTID","exon_number"]
    gtf = pd.concat([fow_annot, rev_annot])

    ## Make intermediate file
    short_form = open(args.input+"short_output.txt", "w")
    short_form.write("ID""\t"+"gene_symbol"+"\t""E1s""\t""E1e""\t""E2s""\t""E2e""\t"
                    "E3s""\t""E3e""\t""E4s""\t""E4e""\t""whole_diff""\t""best match tx""\t"
                    "start_exon""\t""first_diff""\t"
                    "second_diff""\t""third_diff""\t"
                    "fourth_diff""\t""event""\n")
    long_form = open(args.input+"long_output.txt", "w")
    long_form.write("ID""\t"+"gene_symbol"+"\t""E1s""\t""E1e""\t""E2s""\t""E2e""\t"
                    "E3s""\t""E3e""\t""E4s""\t""E4e""\t""whole_diff""\t""best match tx""\t"
                    "start_exon""\t""first_diff""\t""second_diff""\t"
                    "third_diff""\t""fourth_diff""\t"
                    "event""\n")

    for gene in event_exon[geneID].unique():
        temp = event_exon[event_exon[geneID]==gene]
        for line in range(temp.shape[0]):   # Per each event
            Exon_comp(temp.iloc[line,:], temp.iloc[line,2], temp.iloc[line,1])

    short_form.close()
    long_form.close()
