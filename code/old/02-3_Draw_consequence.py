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
plt.rcParams.update({'axes.titlesize': 15,
                     'pdf.fonttype': 42,
                     'ps.fonttype': 42,
                     'font.family': 'Arial'})

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
 
 SSSSS   PPPP    L       III  CCC  EEEEE  DDDDD   EEEEE  CCCCC  OOO  DDDDD  EEEEE  RRRRR
S        P   P   L        I  C   C E      D   D  E      C     O   O D   D  E      R    R
 SSSS    PPPP    L        I  C     EEEE   D   D  EEEE   C     O   O D   D  EEEE   RRRRRR
    S    P       L        I  C   C E      D   D  E      C     O   O D   D  E      R   R
 SSSSS   P       LLLLL  III   CCC  EEEEE  DDDDD   EEEEE  CCCCC  OOO  DDDDD  EEEEE  R    R

Description
    #########################################################
    This script draws transcript and domain level differences 
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains merged bed files', 
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
        ##############################
        ## Remove unwanted information
        w_pfam = w_pfam[~w_pfam[4].isin(remove_info)]
        wo_pfam = wo_pfam[~wo_pfam[4].isin(remove_info)]
        ##############################

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
                    *domain (list): Start, Stop codon, and changed postion
                """
                ########################################################################
                ## Mark a DS event
                if sub_splicing[0] == "RI":
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between([span[1],span[1]-2], ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)
                
                elif sub_splicing[0] in ["EI", "Alt_A3SS", "Alt_A5SS"]:
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between([span[0],span[0]-2], ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)
                
                elif sub_splicing[0].startswith("MXE"):
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)

                else:
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between([span[0],span[0]-2], ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#0000FF", facecolor=exon_color,
                                            zorder=4, linestyle="-", linewidth=2)
                ########################################################################
                
                ## Mark Start and Stop codon
                if len(domain) >= 2:
                    ########################################################################
                    # canvas.scatter(domain[0], 1.6, c="#07FF30", linewidth=1.5, marker="v", zorder=4)
                    # canvas.scatter(domain[1], 1.6, c="#FF1414", linewidth=1.5, marker="v", zorder=4)
                    if span[0] < domain[0] and span[1] < domain[0]:
                        canvas.fill_between(span, ylims['exon_min']+0.5,
                                            ylims['exon_max']-0.5,
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                    
                    elif span[0] < domain[0] and span[1] > domain[0]:
                        canvas.fill_between([span[0],domain[0]], ylims['exon_min']+0.5, 
                                            ylims['exon_max']-0.5,
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        canvas.fill_between([domain[0],span[1]], ylims['exon_min']+0.5, 
                                            ylims['exon_max']-0.5,
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        
                    elif span[0] < domain[1] and span[1] > domain[1]:
                        canvas.fill_between([span[0],domain[1]], ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="#000000", facecolor=exon_color,
                                            zorder=3, linewidth=0.5)
                        canvas.fill_between([domain[1],span[1]], ylims['exon_min']+0.5, 
                                            ylims['exon_max']-0.5,
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        
                    elif span[0] > domain[1]:
                        canvas.fill_between(span, ylims['exon_min']+0.5, 
                                            ylims['exon_max']-0.5,
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
                                fontsize=6, zorder=4, transform=canvas.transAxes)   # Domain name
                else:
                    canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                        edgecolor="#FFFFFF", facecolor=exon_color,
                                        zorder=3, linewidth=0.3)

                return True

            def draw_plot(canvas):
                set_limits()
                domain_num = 0
                for i in range(len(exonpos)):
                    ########################################################################
                    ## Does it need?
                    # if i > 0 and len(exonpos) > 0:   # for first exon
                        # draw_intron([exonpos[i-1][1], exonpos[i][0]], canvas)
                    ########################################################################
                    ## In case of protein domain
                    if len(pos_info) == 1:
                        if domain_num == 0:
                            draw_exon(exonpos[i], canvas, pos_info[0][domain_num])
                        else:
                            draw_exon(exonpos[i], canvas)
                        domain_num += 1
                    elif len(pos_info) > 1 and len(pos_info) <= 2:   # In case of exon structure (Ref)
                        draw_exon(exonpos[i], canvas, *pos_info)
                    elif len(pos_info) >= 3:   # In case of exon structure with changed position (Sim tx)
                        if i == pos_info[2]:
                            draw_exon(exonpos[i], canvas, *pos_info)
                        else:
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
    ## Remove Biding site in figures
    w_pfam = w_pfam[(w_pfam[14]==sid) &
                    (w_pfam[11]==sub_splicing[1]) &
                    (w_pfam[13]==sub_splicing[0]) &
                    (w_pfam[4]!="binding")]
    wo_pfam = wo_pfam[(wo_pfam[14]==sid) &
                      (wo_pfam[11]==sub_splicing[1]) &
                      (wo_pfam[13]==sub_splicing[0]) &
                      (wo_pfam[4]!="binding")]
    n_domain = wo_pfam[3].unique()

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
            if fig_input.shape[0] == 0:
                return
            else:
                n_row = 0
                fig = plt.figure(figsize=(10,(len(n_domain)/3)+1))
                for num, kinds_domain in enumerate(range(len(n_domain))):   # Make Figure and dataframe "one domain at once"
                    ## Prepare figure params
                    gene = fig_input["ID"].str.split(";").str[1].values[0]
                    tx = fig_input["ENST"].values[0]
                    eventid = fig_input["event_type"].values[0]

                    ## Make an input for visualization
                    ref_bed, ref_domain, ref_dname = Make_query(ref, wo_pfam[wo_pfam[3]==n_domain[kinds_domain]])
                    sim_bed, sim_domain, sim_dname = Make_query(sim_tx, w_pfam[w_pfam[3]==n_domain[kinds_domain]])
                    if len(ref_bed) > 0:    # Markup changed regions
                        Make_fig = 0
                        for r,s in zip(ref_bed,sim_bed):
                            if r[0] != s[0] or \
                            r[1] != s[1]:
                                event_pos_ref = ref_bed.index(r)
                                event_pos_query = sim_bed.index(s)
                                Make_fig += 1
                                break
                        
                        ## Make a dataframe for subsequent analysis
                        def Prunning(domain_bed, dname, st, end):
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

                            return p_domain_bed, p_dname, l_domain
                        
                        p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                           fig_input["Start"].values[0], 
                                                                           fig_input["STOP"].values[0])
                        p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                           fig_input["simStart"].values[0], 
                                                                           fig_input["simSTOP"].values[0])
                        
                        if len(p_ref_domain) > 0 and fig_on == "on":    # TO DO it should update, for NMD recov
                            if n_row == 0:
                                ref_axs = plt.subplot(2+(len(n_domain)*2), 1, 1)
                                if Make_fig == 1:
                                    Comp_plot(ref_axs, ref_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0], event_pos_ref)
                                else:
                                    Comp_plot(ref_axs, ref_bed, "#FF7878", fig_input["Start"].values[0], fig_input["STOP"].values[0])
                                plt.text(-0.03,0.2,"{}".format("Ref-Tx"), 
                                        transform=ref_axs.transAxes,
                                        fontsize=10)
                                
                                sim_axs = plt.subplot(2+(len(n_domain)*2), 1, 1+len(n_domain)+1, sharex=ref_axs)
                                Comp_plot(sim_axs, sim_bed, "#73ABFF", fig_input["simStart"].values[0], fig_input["simSTOP"].values[0], event_pos_query)
                                plt.text(-0.03,0.2,"{}".format("Sim-Tx"), 
                                        transform=sim_axs.transAxes,
                                        fontsize=10)
                                n_row += 1


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
        
        Draw_splicemap(sub,sub_ref,sub_sim_tx, "on")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/"):
        args.input = args.input+"/"
    
    ## Run code
    # Run("CA","SRSF2","EI","ENST00000359995.9")
    Run(args.splicing_event, args.gene, args.sim_splicing_event, args.transcript)
# %%
