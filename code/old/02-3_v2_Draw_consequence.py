#!/usr/bin/env python3
# %%
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import argparse
import os
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({
'axes.titlesize': 15,     # 제목 글꼴 크기
'axes.labelsize': 14,     # x, y축 라벨 글꼴 크기
'xtick.labelsize': 12,    # x축 틱 라벨 글꼴 크기
'ytick.labelsize': 12,    # y축 틱 라벨 글꼴 크기
'legend.fontsize': 12,    # 범례 글꼴 크기
'figure.titlesize': 15    # figure 제목 글꼴 크기
})

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(
        description='''This script discoveries match isoforms with given differential splicing events''')
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains simulated results', 
                        required=True,
                        type=str)
    parser.add_argument('--splicing_event', '-s', 
                        help='Splicing categories e.g., CA, RI, A3SS, A5SS, MXE', 
                        required=True,
                        type=str)
    parser.add_argument('--gene', '-g', 
                        help='Gene symbol', 
                        required=True,
                        type=str)
    parser.add_argument('--sim_splicing_event', '-sim', 
                        help='Simulated splicing e.g., EI, ES, RI, SI, Ori_A3SS, Alt_A3SS MXE1, MXE2', 
                        required=True,
                        type=str)
    parser.add_argument('--transcript', '-t', 
                        help='Reference transcript that contains same exon structure', 
                        required=True,
                        type=str)

    ## Optional
    parser.add_argument('--remove_info', '-ri', nargs='*',
                        help='DO NOT DEPIC these functional category e.g., coiled chain', 
                        type=str,
                        default=[])                   
    
    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

args, parser = parse_args(sys.argv[1:])
args.input = args.input+"/"

def Run(key, input_gene_name, *sub_splicing):
    """ Draw splicing consequence

    Args:
        key (_type_): Splicing categories (CA, RI, A3SS ...)
        input_gene_name (_type_): Gene_symbol
        sub_splicing (_type_): ES, EI, RI, SI ...

    Returns:
        _type_: _description_
    """
    
    def Load_data(key, remove_info):
        """ Load information data
            It could be smaller (perhaps using cmd?)
        Args:
            key (_type_): Splicing category
            remove_info (_type_): Skipped events

        Returns:
            _type_: _description_
        """
        sim_tx = pd.read_csv(args.input+"merged_w_"+key+".bed",
                                sep="\t", header=None)
        ref = pd.read_csv(args.input+"merged_wo_"+key+".bed",
                                sep="\t", header=None)
        query = pd.read_csv(args.input+"result/all_{}_Main_output.txt".format(key),
                            sep="\t", skiprows=1)
        
        query = query[query["5'UTR"]!="Frame loss"] # Filter out frame loss, because they donot have domain
        query = query[query["ORF_priority"]=="pORF1"]
        query["ORF_priority"] = query["ORF_priority"].apply(lambda x : 3 if x =="pORF1" else (2 if x == "pORF2" else 1))
        query = query[["LongID","occurred_event","Start","Stop","dAA","Target_TX","ORF_priority","Domain_change_ratio"]]
        query["simStop"] = query["Stop"].str.split("-").str[1]
        query["Stop"] = query["Stop"].str.split("-").str[0]
        query["simStart"] = query["Start"].str.split("-").str[1]
        query["Start"] = query["Start"].str.split("-").str[0]
        query.columns = ["ID","event_type","Start","STOP","dAA","ENST","ORF","delta L","simSTOP","simStart"]
        query["SID"] = query["ID"]
        query["ID"] = query["ID"]+"|"+query["ENST"]+"|"+query["Start"]+"|"+query["simStart"]
        query[["simSTOP","STOP","Start","simStart"]] = query[["simSTOP","STOP","Start","simStart"]].astype(float)

        w_pfam = pd.read_csv(args.input+"table/{}_w_Pfam.txt".format(key),
                        sep="\t", header=None)
        wo_pfam = pd.read_csv(args.input+"table/{}_wo_Pfam.txt".format(key),
                        sep="\t", header=None)

        w_pfam = w_pfam[~w_pfam[4].isin(remove_info)]
        wo_pfam = wo_pfam[~wo_pfam[4].isin(remove_info)]

        w_pfam[3] = w_pfam[3] +" ("+w_pfam[4]+")"
        wo_pfam[3] = wo_pfam[3] +" ("+wo_pfam[4]+")"

        return sim_tx, ref, query, w_pfam, wo_pfam


    sim_tx, ref, query, w_pfam, wo_pfam = Load_data(key, args.remove_info)

    if len(sub_splicing) > 0 and len(sub_splicing) < 2:
        query = query[query["event_type"]==sub_splicing[0]]
    elif len(sub_splicing) > 1:
        query = query[query["ENST"]==sub_splicing[1]]
    else:
        pass
    query["gene"] = query["ID"].str.split(";").str[1]
    ################## INPUT CHECK
    # print(query)
    ##################
    def Comp_plot(canvas, exonpos, exon_color, *pos_info):
        """ Draw comparison plot

        Args:
            canvas (_type_): subplots axs
            exonpos (list): exon position
            exon_color (_type_): exon block color
            pos_info (list): start, stop_codon, and different exon
        Returns:
            _type_: Figure
        """
        sns.despine(left=True, top=True, bottom=True)
        canvas.tick_params(left=False, bottom=False)
        canvas.set_xticks([])
        canvas.set_yticks([])

        bgcolor="#FFFFFF"  # background
        if len(exonpos) > 0:  # If there are available exon coordinates
            ylims = {'exon_max':2.0, 'exon_min':.5}   # The height of exon box
            canvas.set_facecolor(bgcolor)
            # totalSpan = exonpos[-1][1] - exonpos[0][0]  # End of last exon - start of first exon
            # minExonLen = totalSpan*0.1 # isoform length

            def set_limits():
                ylims['bar_min'] = ylims['exon_max']+0.2
                ylims['bar_max'] = ylims['bar_min']+(ylims['exon_max']-ylims['exon_min'])/5.0

            def draw_exon(span, canvas, *domain):
                """ Draw exon

                Args:
                    span (_type_): [exon start, exon end]
                    canvas (_type_): axs
                    *domain (list): Start, Stop codon, and changed postion, len(domain) >2 means target exon
                """
                ########################################################################
                ## Mark a DS event >> It should be updated for better upderstanding
                if sub_splicing[0] == "RI":
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between([span[1],span[1]], ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor="#FFFFFF",
                                            zorder=4, linestyle=":", linewidth=1.5)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor="none",
                                            zorder=4, linestyle="-", linewidth=1.5)
                
                elif sub_splicing[0] in ["EI", "Alt_A3SS", "Alt_A5SS"]:
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor="#FFFFFF",
                                            zorder=4, linestyle=":", linewidth=1.5)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor="none",
                                            zorder=4, linestyle="-", linewidth=1.5)
                
                elif sub_splicing[0].startswith("MXE"):
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        if domain[3] == "on":
                            fc = "#FFFFFF"
                        else:
                            fc = "none"
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor=fc,
                                            zorder=4, linestyle=":", linewidth=1.5)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        if domain[3] == "on":#0A121D
                            fc = "none"
                        else:
                            fc = "#FFFFFF"
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor=fc,
                                            zorder=4, linestyle=":", linewidth=1.5)

                else:   # All skiping event
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor='none',
                                            zorder=4, linestyle="-", linewidth=1.5)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor='#FFFFFF',
                                            zorder=4, linestyle=":", linewidth=1.5)
                ########################################################################
                
                ## Mark Start and Stop codon
                if len(domain) >= 2:
                    ########################################################################
                    # canvas.scatter(domain[0], 1.6, c="#07FF30", linewidth=1.5, marker="v", zorder=4)
                    # canvas.scatter(domain[1], 1.6, c="#FF1414", linewidth=1.5, marker="v", zorder=4)
                    if span[0] < domain[0] and span[1] < domain[0]: # It is totally 5' UTR 
                        canvas.fill_between(span, ylims['exon_min'],
                                            ylims['exon_max'],
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                    
                    elif span[0] < domain[0] and span[1] >= domain[0]: # It is partially 5' UTR
                        canvas.fill_between([span[0],domain[0]], ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="none", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        canvas.fill_between([domain[0],span[1]], ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="none", facecolor=exon_color,
                                            zorder=3, linewidth=0.5)
                        canvas.plot([span[0], domain[0]], [ylims['exon_min'], ylims['exon_min']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], domain[0]], [ylims['exon_max'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[0], span[1]], [ylims['exon_min'], ylims['exon_min']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[0], span[1]], [ylims['exon_max'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], span[0]], [ylims['exon_min'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[1], span[1]], [ylims['exon_min'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        
                    elif span[0] < domain[1] and span[1] > domain[1]:
                        canvas.fill_between([span[0],domain[1]], ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="none", facecolor=exon_color,
                                            zorder=3, linewidth=0.5)
                        canvas.fill_between([domain[1],span[1]], ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="none", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        canvas.plot([span[0], domain[1]], [ylims['exon_min'], ylims['exon_min']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], domain[1]], [ylims['exon_max'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[1], span[1]], [ylims['exon_min'], ylims['exon_min']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[1], span[1]], [ylims['exon_max'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], span[0]], [ylims['exon_min'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[1], span[1]], [ylims['exon_min'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)

                        
                    elif span[0] > domain[1]:
                        canvas.fill_between(span, ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        
                    else:
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor=exon_color,
                                            zorder=3, linewidth=0.5)
                    ########################################################################
                
                ## Draw domina block
                elif len(domain) == 1:
                    canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                        edgecolor="#FFFFFF", facecolor=exon_color,
                                        zorder=3, linewidth=0.3)
                    canvas.text(0.9,0.4,"{}".format(domain[0]),
                                fontsize=6, zorder=4, transform=canvas.transAxes)   # Add domain name
                else:
                    canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                        edgecolor="#FFFFFF", facecolor=exon_color,
                                        zorder=3, linewidth=0.3)

                return True

            def draw_plot(canvas):
                set_limits()
                domain_num = 0  # Set start point
                mxe_exon = 0    # For MXE exon marking
                for i in range(len(exonpos)):   # For each exon block
                    ########################################################################
                    ## Does it need?
                    # if i > 0 and len(exonpos) > 0:   # for first exon
                        # draw_intron([exonpos[i-1][1], exonpos[i][0]], canvas)
                    ########################################################################
                    ## Drawing the protein domain
                    if len(pos_info) == 1:  # If Sim,Ref-TX has CDS differences
                        if domain_num == 0:
                            draw_exon(exonpos[i], canvas, pos_info[0][domain_num])
                        else:
                            draw_exon(exonpos[i], canvas)
                        domain_num += 1
                    # Drawing the unchanged exon structure
                    elif len(pos_info) > 1 and len(pos_info) <= 2:
                        draw_exon(exonpos[i], canvas, *pos_info)
                    # Drawing the exon structure with changes
                    elif len(pos_info) >= 3:
                        if sub_splicing[0].startswith("MXE"):
                            if i in [pos_info[2], pos_info[2]+1]:    # Changed position
                                new_pos_info = [i for i in pos_info]
                                if mxe_exon == 0:   # For 1st MXE
                                    new_pos_info.append("on")
                                else:   # For 2nd MXE
                                    new_pos_info.append("off")

                                draw_exon(exonpos[i], canvas, *new_pos_info)
                                mxe_exon += 1
                            else:   # If it is not changed exon
                                draw_exon(exonpos[i], canvas, pos_info[0], pos_info[1])
                        
                        else:   # For non MXE cases
                            if i == pos_info[2]:  # It has changed exon
                                draw_exon(exonpos[i], canvas, *pos_info)
                            else:   # If it is not changed exon
                                draw_exon(exonpos[i], canvas, pos_info[0], pos_info[1])

                canvas.fill_between([exonpos[0][0], exonpos[-1][1]],
                                    ylims['bar_min'], ylims['bar_max'],
                                    edgecolor=bgcolor, facecolor=bgcolor)

            draw_plot(canvas)


    def Make_query(final, pfam):
        """ Collect genomic coordinate information from query event

        Args:
            final (dataframe): 01-4 output (merged_bed)

        Returns:
            new_line : Input of the drawing figure function
        """
        diff = []
        new_line = []
        domain_line = []
        domain_name = []
        ##############################
        strand = final[5].values[0]
        for line in range(final.shape[0]):  ## Individual transcript
            if line > 0:    ## From 2nd exon
                if strand == "+":
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                  (pfam[9]==final.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                    intron = np.abs(int(final.iloc[line-1,2] - int(final.iloc[line,1])))
                    diff.append(intron)
                    new_line.append([int(final.iloc[line,1])-np.sum(diff)-zero, 
                                    int(final.iloc[line,2])-np.sum(diff)-zero])
                    if domain.shape[0] > 0: ## There are multi protein domains, it is not only one
                        domain_line.append([int(domain.iloc[0,1])-np.sum(diff)-zero,
                                            int(domain.iloc[0,2])-np.sum(diff)-zero,])
                        if not domain.iloc[0,3] in domain_name:
                            domain_name.append(domain.iloc[0,3])
                        else:
                            domain_name.append(" ")
            
                else:   # Start with -2, because first exon is -1
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                (pfam[9]==final.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                    intron = np.abs(np.abs(int(final.iloc[line-1,1])-zero) - np.abs(int(final.iloc[line,2])-zero))
                    diff.append(intron)
                    new_line.append([np.abs(int(final.iloc[line,2])-zero)-np.sum(diff), 
                                    np.abs(int(final.iloc[line,1])-zero)-np.sum(diff)])    # 2*zero, represents to zero and then direction change
                    
                    if domain.shape[0] > 0: ## There are multi protein domains, it is not only one
                        domain_line.append([np.abs(int(domain.iloc[0,2])-zero)-np.sum(diff),
                                            np.abs(int(domain.iloc[0,1])-zero)-np.sum(diff),])
                        if not domain.iloc[0,3] in domain_name:
                            domain_name.append(domain.iloc[0,3])
                        else:
                            domain_name.append(" ")

            else:   ## For 1st exon
                if strand == "+":
                    zero = int(final.iloc[line,1])  # The most left exon's start site for align in figure
                    new_line.append([int(final.iloc[line,1])-zero,
                                    int(final.iloc[line,2])-zero])
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                (pfam[9]==final.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                    
                    if domain.shape[0] > 0:
                        domain_line.append([np.abs(int(domain.iloc[0,1])-zero),
                                            np.abs(int(domain.iloc[0,2])-zero),])
                        if not domain.iloc[0,3] in domain_name:
                            domain_name.append(domain.iloc[0,3])
                        else:
                            domain_name.append(" ")
                else:
                    zero = int(final.iloc[line,2])  # The most right exon's end site
                    new_line.append([np.abs(int(final.iloc[line,2])-zero),
                                    np.abs(int(final.iloc[line,1])-zero)])
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                (pfam[9]==final.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                    if domain.shape[0] > 0:
                        domain_line.append([np.abs(int(domain.iloc[0,2])-zero),
                                            np.abs(int(domain.iloc[0,1])-zero),])
                        if not domain.iloc[0,3] in domain_name:
                            domain_name.append(domain.iloc[0,3])
                        else:
                            domain_name.append(" ")
        
        return new_line, domain_line, domain_name

    test_list = query.iloc[:,0].unique()
    sid = test_list[0].split("|")[0]
    test_list = query[query["ID"].str.contains(input_gene_name)]["ID"].unique()
    ## Add stop function
    
    ################################ TEMP
    ## Remove Biding site in figures, it spends lots of space
    ## Capture target domains which are contained in target transcript
    w_pfam = w_pfam[(w_pfam[14]==sid) &
                    (w_pfam[11]==sub_splicing[1]) &
                    (w_pfam[13]==sub_splicing[0]) &
                    (w_pfam[4]!="binding")]
    wo_pfam = wo_pfam[(wo_pfam[14]==sid) &
                      (wo_pfam[11]==sub_splicing[1]) &
                      (wo_pfam[13]==sub_splicing[0]) &
                      (wo_pfam[4]!="binding")]
    n_domain = wo_pfam[3].unique()
    ################################

    OUT = args.input+"figure/"
    for id in test_list:    # Certain splicing event with one of the possible ORF
        sub = query[query["ID"]==id]
        deltaL = float(sub["delta L"].values[0])
        sub_sim_tx = sim_tx[(sim_tx[4]==sub["ENST"].values[0]) &
                            (sim_tx[6]==sub["event_type"].values[0]) &
                            (sim_tx[7]==sub["SID"].values[0])].copy() # ASS events can be distinguished by event type, but other event need ID
        sub_ref = ref[(ref[4]==sub["ENST"].values[0]) & 
                      (ref[6]==sub["event_type"].values[0])&
                      (ref[7]==sub["SID"].values[0])].copy()

        def Draw_splicemap(fig_input,ref,sim_tx,fig_on):
            if fig_input.shape[0] == 0: # If input is None
                return
            else:
                n_row = 0
                fig = plt.figure(figsize=(8,(len(n_domain)/3)+1.4))
                for num, kinds_domain in enumerate(range(len(n_domain))):   # Make Figure and dataframe "one domain at once"
                    ## Prepare figure params
                    gene = fig_input["ID"].str.split(";").str[1].values[0]
                    tx = fig_input["ENST"].values[0]
                    eventid = fig_input["event_type"].values[0]

                    ## Make an input (bed) for visualization
                    ref_bed, ref_domain, ref_dname = Make_query(ref, wo_pfam[wo_pfam[3]==n_domain[kinds_domain]])
                    sim_bed, sim_domain, sim_dname = Make_query(sim_tx, w_pfam[w_pfam[3]==n_domain[kinds_domain]])
                    if len(ref_bed) > 0:    # Markup changed regions
                        Make_fig = 0
                        for r,s in zip(ref_bed,sim_bed):
                            if r[0] != s[0] or \
                            r[1] != s[1]:
                                event_pos_ref = ref_bed.index(r)
                                event_pos_query = sim_bed.index(s)  # Posion of sim DS
                                Make_fig += 1   # Make altered regions due to the given DS event
                                break   # there is at least one different CDS alterations)
                        
                        ## Make a dataframe for subsequent analysis
                        def Prunning(domain_bed, dname, st, end, diff_exon):
                            ## Add truncated domain
                            ## Motif could be located on UTR regions
                            l_domain = []   # length
                            p_domain_bed = []   # domain position
                            p_dname = []    # domain name
                            for index in range(len(domain_bed)):
                                if domain_bed[index][0] >= st and \
                                    domain_bed[index][1] <= end:
                                    l_domain.append((domain_bed[index][1] - domain_bed[index][0]))
                                    p_domain_bed.append(domain_bed[index])
                                    p_dname.append(dname[index])

                                elif domain_bed[index][0] < end and \
                                    domain_bed[index][1] > end: # Partial domain
                                    l_domain.append((end - domain_bed[index][0]))
                                    p_domain_bed.append([domain_bed[index][0],end])
                                    p_dname.append(dname[index])

                                elif domain_bed[index][0] < st and \
                                    domain_bed[index][1] > st: # Partial domain
                                    l_domain.append((domain_bed[index][1] - st))
                                    p_domain_bed.append([st,domain_bed[index][1]])
                                    p_dname.append(dname[index])

                                if domain_bed[index][0] >= diff_exon[0]:    # if domain is located after diff exon
                                    space = np.abs(diff_exon[0]-diff_exon[1])
                                    domain_bed[index][0] = domain_bed[index][0]+space
                                    domain_bed[index][1] = domain_bed[index][1]+space

                            return p_domain_bed, p_dname, l_domain

                        if fig_input["event_type"].unique()[0] in ["EI","Alt_A3SS","Alt_A5SS"]:
                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               fig_input["Start"].values[0], 
                                                                               fig_input["STOP"].values[0],
                                                                               sim_bed[event_pos_query])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               fig_input["simStart"].values[0], 
                                                                               fig_input["simSTOP"].values[0],
                                                                               [0,0])
                            
                        elif fig_input["event_type"].unique()[0] in ["ES","Ori_A3SS","Ori_A5SS"]:
                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               fig_input["Start"].values[0], 
                                                                               fig_input["STOP"].values[0],
                                                                               [0,0])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               fig_input["simStart"].values[0], 
                                                                               fig_input["simSTOP"].values[0],
                                                                               ref_bed[event_pos_query])

                        elif fig_input["event_type"].unique()[0] in ["RI","SI"]:
                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               fig_input["Start"].values[0], 
                                                                               fig_input["STOP"].values[0],
                                                                               [0,0])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               fig_input["simStart"].values[0], 
                                                                               fig_input["simSTOP"].values[0],
                                                                               [0,0])

                        else:   # MXE 1 or 2, mxe_bed should be differ
                            ## If Sim MXE2
                            mxe2_leng = np.abs(sim_bed[event_pos_query][0] - sim_bed[event_pos_query][1])
                            mxe_bed = ref_bed[:event_pos_ref+1] + \
                                      [[sim_bed[event_pos_ref][0],sim_bed[event_pos_ref][1]]] + \
                                      [(i+mxe2_leng) for i in ref_bed[event_pos_ref+1:]]
                            
                            ## If Sim MXE1
                            mxe1_leng = np.abs(sim_bed[event_pos_query][0] - sim_bed[event_pos_query][1])
                            mxe_bed = ref_bed[:event_pos_ref] + \
                                      [[sim_bed[event_pos_query][0],sim_bed[event_pos_query][1]]] + \
                                      [(i+mxe1_leng) for i in ref_bed[event_pos_ref:]]
                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               fig_input["Start"].values[0], 
                                                                               fig_input["STOP"].values[0],
                                                                               sim_bed[event_pos_query])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               fig_input["simStart"].values[0], 
                                                                               fig_input["simSTOP"].values[0],
                                                                               ref_bed[event_pos_query])
                            
                        if len(p_ref_domain) > 0 and fig_on == "on":
                            if n_row == 0:
                                ref_axs = plt.subplot(2+(len(n_domain)*2), 1, 1)
                                if Make_fig == 1:
                                    ## If sim splicing is inclusion form, ref tx should have same tx structure
                                    a = np.abs(ref_bed[-1][1]-sim_bed[-1][1])
                                    if fig_input["event_type"].unique()[0] in ["EI","Alt_A3SS","Alt_A5SS"]:
                                        Comp_plot(ref_axs, sim_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0]+a, event_pos_ref)
                                    ## If sim splicing is skipped form, ref tx should have ref tx structure
                                    elif fig_input["event_type"].unique()[0] in ["ES","Ori_A3SS","Ori_A5SS"]:
                                        Comp_plot(ref_axs, ref_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0], event_pos_ref)
                                    elif fig_input["event_type"].unique()[0] in ["RI","SI"]:   # RI and SI
                                        Comp_plot(ref_axs, ref_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0], event_pos_ref)
                                    else:   # MXE >> MXE should have update stop, because it has two exons
                                        Comp_plot(ref_axs, mxe_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0], event_pos_ref)
                                else:
                                    if fig_input["event_type"].unique()[0] in ["EI","Alt_A3SS","Alt_A5SS"]:
                                        Comp_plot(ref_axs, sim_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0]+a)
                                    else:
                                        Comp_plot(ref_axs, ref_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0])
                                plt.text(-0.03,0.2,"{}".format("Ref-Tx"), 
                                        transform=ref_axs.transAxes,
                                        fontsize=10)
                                
                                sim_axs = plt.subplot(2+(len(n_domain)*2), 1, 1+len(n_domain)+1, sharex=ref_axs)
                                if fig_input["event_type"].unique()[0] in ["ES","Ori_A3SS","Ori_A5SS"]:
                                    b = np.abs(ref_bed[-1][1]-sim_bed[-1][1])
                                    Comp_plot(sim_axs, ref_bed, "#73ABFF", fig_input["simStart"].values[0], fig_input["simSTOP"].values[0]+b, event_pos_query)
                                elif fig_input["event_type"].unique()[0] in ["EI","RI","SI","Alt_A3SS","Alt_A5SS"]:
                                    Comp_plot(sim_axs, sim_bed, "#73ABFF", fig_input["simStart"].values[0], fig_input["simSTOP"].values[0], event_pos_query)
                                else:   # MXE
                                    Comp_plot(sim_axs, mxe_bed, "#73ABFF", fig_input["simStart"].values[0], fig_input["simSTOP"].values[0], event_pos_query)
                                plt.text(-0.03,0.2,"{}".format("Sim-Tx"), 
                                        transform=sim_axs.transAxes,
                                        fontsize=10)
                                n_row += 1  # For next domain, tx drawing function will be skipped


                            ## Draw protein domain
                            ref_nmd = ref_bed[-1][0] - fig_input["STOP"].values[0]
                            sim_nmd = sim_bed[-1][0] - fig_input["simSTOP"].values[0]
                            domain_pal = {"diff":"#E72B2B", "NMD":"#C6C6C6", "same":"#69EB89"}
                            
                            def Domain_plot(ref_color_key,sim_color_key):
                                axs = plt.subplot(2+(len(n_domain)*2), 1, n_row+1, sharex=ref_axs)
                                Comp_plot(axs, p_ref_domain, domain_pal[ref_color_key], [n_domain[kinds_domain]])    # Think about dict color map with domain name
                                axs2 = plt.subplot(2+(len(n_domain)*2), 1, n_row+len(n_domain)+2, sharex=ref_axs)
                                Comp_plot(axs2, p_sim_domain, domain_pal[sim_color_key], [n_domain[kinds_domain]])    # Think about dict color map with domain name
 
                            if sim_nmd >= 55 and ref_nmd >= 55:
                                Domain_plot("NMD","NMD")
                            elif ref_nmd >= 55 and sim_nmd < 55:
                                Domain_plot("NMD","diff")
                            elif ref_nmd < 55 and sim_nmd >= 55:
                                Domain_plot("diff","NMD")
                            elif np.abs(np.sum(l_ref_domain) - np.sum(l_sim_domain)) == 0:
                                Domain_plot("same","same")
                            elif np.abs(np.sum(l_ref_domain) - np.sum(l_sim_domain)) != 0: 
                                Domain_plot("diff","diff")
                            n_row += 1


                            fig.suptitle("{}; {}; sim{}; Delta L:{}".format(gene, tx, eventid, round(deltaL,3)), y=1.)

                        else:
                            print("refTX does not have domain")
                    else:
                        plt.close() # Temporal
                
                plt.tight_layout(h_pad=0)
                if not os.path.exists(OUT+"consequence/"):
                    os.mkdir(OUT+"consequence/")

                if len(sub_splicing) > 1:
                    plt.savefig(OUT+"consequence/{}_{}_splicing_map.pdf".format(fig_input["SID"].values[0],sub_splicing[1]),
                                bbox_inches="tight")
                elif len(sub_splicing) == 1:
                    plt.savefig(OUT+"consequence/{}_{}_splicing_map.pdf".format(fig_input["SID"].values[0],sub_splicing[0]),
                                bbox_inches="tight")
                else:
                    plt.savefig(OUT+"consequence/{}_splicing_map.pdf".format(fig_input["ID"].values[0]),
                                bbox_inches="tight")

        if len(n_domain) > 0:
            Draw_splicemap(sub,sub_ref,sub_sim_tx, "on")
        else:
            print("No mapped functional regions")

Run(args.splicing_event, args.gene, args.sim_splicing_event, args.transcript)
