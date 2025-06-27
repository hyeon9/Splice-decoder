#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"
annot_DIR = "/Users/kangh/Downloads/annot/human/translation_pipeline/"
exon_dict = pd.read_csv(annot_DIR+"v32_exon_only.final.gtf",
                        sep="\t")
# DIR = "/Users/kangh/Documents/Analysis_note/isoform_translation/Young_old/"
DIR = "/Users/kangh/Documents/Analysis_note/isoform_translation/MYC_paper/"
query = pd.read_csv(DIR+"fig_input.txt", 
                    sep="\t")
gene = "CHEK2"
fig_type = "Full" # Zoomed
# DIR = "/Users/kangh/Documents/Analysis_note/isoform_translation/torpor_project/"
# exon_dict = pd.read_csv(DIR+"exon_only.final.gtf",
#                         sep="\t")
# %%
def Comp_plot(canvas, exonIntervals, exon_color, dom_id):
    """ Draw comparison plot

    Args:
        canvas (_type_): subplots axs
        exonIntervals (list): exon position
        exon_color (_type_): exon block color
        dom_id (list): domain id or nt

    Returns:
        _type_: Figure
    """
    ## To do?
    ## Set the same start site, when I conduct AA sequence comparison

    edgecolor="#ADADAD"
    bgcolor="#FFFFFF"  # background

    if len(exonIntervals) > 0:  # If there are available exon coordinates
        ylims = {'exon_max':1, 'exon_min':.5}   # The height of exon box
        canvas.set_facecolor(bgcolor)
        totalSpan = exonIntervals[-1][1] - exonIntervals[0][0]  # End of last exon - start of first exon
        minExonLen = totalSpan*0.1 # isoform length

        def set_limits():
            # ylims['intron_max'] = ylims['exon_max']*0.9    # Using these combi, splicing line can be added
            # ylims['intron_min'] = (ylims['exon_max'] + ylims['exon_min'])/2.0
            ylims["intron_max"] = ylims['exon_max']*.75
            ylims["intron_min"] = ylims['exon_max']*.75
            ylims['bar_min'] = ylims['exon_max']+0.2
            ylims['bar_max'] = ylims['bar_min']+(ylims['exon_max']-ylims['exon_min'])/5.0

        def draw_exon(span, canvas, *domain):
            """ Draw exon

            Args:
                span (_type_): [exon start, exon end]
                canvas (_type_): axs

            """
            if "target" in domain[0]:
                canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                    edgecolor=edgecolor, facecolor=exon_color,
                                    zorder=3, linestyle="-")
            else:
                canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                        edgecolor=edgecolor, facecolor=exon_color,
                                        zorder=3, linestyle="-", linewidth=0)
            # if num == 0:
            #     print(span[0] + span[1])
            #     canvas.text(span[0]+100, .2, 
            #                 "5 prime", fontsize=12, zorder=3)
            
            return True
                
        def draw_intron(span, canvas):
            """ Draw intron

            Args:
                span (_type_): [end of n-1th exon, start of nth exon]
                canvas (_type_): axs

            """
            mid = (span[0]+span[1])/2.0
            canvas.plot([span[0], mid], [ylims['intron_min'], ylims['intron_max']],
                            c="#ADADAD", lw=1.3, ls="-")
            canvas.plot([mid, span[1]], [ylims['intron_max'], ylims['intron_min']],
                            c="#ADADAD", lw=1.3, ls="-")

            return True

        def draw_plot(canvas):
            set_limits()
            # _transform_spans(exonIntervals)   ## line 189
            for i in range(len(exonIntervals)):
                if i > 0 and len(exonIntervals) > 0:   # for first exon
                    draw_intron([exonIntervals[i-1][1], exonIntervals[i][0]], canvas)
                draw_exon(exonIntervals[i], canvas, dom_id)
                # _draw_exon(exonIntervals[i], canvas, dom_id[i])
            canvas.fill_between([exonIntervals[0][0], exonIntervals[-1][1]],
                                        ylims['bar_min'], ylims['bar_max'],
                                        edgecolor=bgcolor, facecolor=bgcolor)

        draw_plot(canvas)
        canvas.tick_params(left=False,
                           bottom=False)
        canvas.set_xticks([])
        canvas.set_yticks([])


def Make_query(final, *scale):
    """ Collect genomic coordinate information from query event

    Args:
        final (dataframe): 01_exon_coordinate.py

    Returns:
        new_line : Input of the drawing figure function
        error : Associate with scaling 1 (error, zoomed instead of scaling), 0 (scaling)
    """
    new_line = []
    final = pd.DataFrame(final).T
    for col in range(2,7,2):
        if final.iloc[0,col] != "SKIP" and \
            final.iloc[0,col] != "RI":
            new_line.append([int(final.iloc[0,col]), 
                             int(final.iloc[0,col+1])])
        
    ## Scaling; same with ref
    if len(scale) > 0:
        scale = scale[0]
        scale = pd.DataFrame.from_dict(scale, orient="index").fillna(0.0)
        scale.loc["sum"] = scale.sum()
        scale.columns = ["first","second","third"]
        new_line[0] = new_line[0] - scale.loc["sum","first"]
        new_line[1] = new_line[1] - scale.loc["sum","second"]
        if len(new_line) > 2:
            new_line[2] = new_line[2] - scale.loc["sum","third"]
            if new_line[2][0] < new_line[1][1]: # the start of third exon is shorter than the end of second exon
                error = 1
            else:
                error = 0
        else:   # If there are only two exon, it doesn't matter?? sure?, same error could be occured in this case
            if new_line[1][0] < new_line[0][1]:
                error=1
            else:
                error = 0
        
        return new_line, error
    
    else:
        return new_line, 0


def Make_ref(target_tx, match_ex, *scale_key):
    """ Make bed using transcript id and exon gtf file

    Args:
        target_tx (_type_): ENST ID
        match_ex (_type_): number of first matched exon

    Returns:
        _type_: bed file and scale factor
    """
    tx = exon_dict[exon_dict["ENSTID"]==target_tx]
    match_ex = int(match_ex.replace("exon",""))
    bed = []
    scale_factor = {}   # with query
    f_scale = 0

    for line in tx.index.tolist():
        st = int(tx.loc[line]["start"])
        ed = int(tx.loc[line]["end"])
        bed.append([st,ed])

    ## Scaling; Find Maximum intron difference and re-scaling
    intron_diff = []
    for pair in range(0,len(bed)-1):
        iss = bed[pair][1]
        ies = bed[pair+1][0]
        intron_length = np.abs(iss-ies)
        intron_diff.append(intron_length)

    if max(intron_diff) > 10000 and len(scale_key) > 0:  # Reducing threshold
        size_factor = max(intron_diff)*0.9  # Reducing size
        pair = intron_diff.index(max(intron_diff))
        
        if (pair+1) - match_ex < 0: # first + second + last
            f_scale += size_factor
            scale_factor['whole'] = []
            scale_factor['whole'].append(f_scale)
            scale_factor['whole'].append(f_scale)
            scale_factor['whole'].append(f_scale)
            bed[pair+1:] = list(bed[pair+1:]-size_factor)
        
        elif (pair+1) - match_ex == 0:  # second + last
            scale_factor['fromsecond'] = []
            scale_factor['fromsecond'].append(0)
            scale_factor['fromsecond'].append(size_factor)
            scale_factor['fromsecond'].append(size_factor)
            bed[pair+1:] = list(bed[pair+1:]-size_factor)

        elif (pair+1) - match_ex == 1:  # last
            scale_factor['last'] = []
            scale_factor['last'].append(0)
            scale_factor['last'].append(0)
            scale_factor['last'].append(size_factor)
            bed[pair+1:] = list(bed[pair+1:]-size_factor)
        
        # else:  # last
        #     scale_factor['whole'] = []
        #     scale_factor['whole'].append(0)
        #     scale_factor['whole'].append(0)
        #     scale_factor['whole'].append(0)
        #     bed = list(bed-size_factor)
    else:
        bed[pair+1:] = list(bed[pair+1:])

    return bed, scale_factor


# %%
# test_list = query[query["ID"].str.contains("ARMCX5-GPRASP2")]["ID"].unique()
# test_list = query[query["ID"].str.contains("KDM6A")]["ID"].unique()
# test_list = query[query["ID"].str.contains("EPB41L1")]["ID"].unique()
test_list = query[query["ID"].str.contains(gene)]["ID"].unique()
# test_list = query.iloc[:,0].unique()
# np.random.seed(2526)  # Example random seed
# rand_num = np.random.randint(0,len(test_list),1)
# test_list = test_list[[rand_num]]
# %%
def Draw_splicemap(fig_input, fig_type):
    """ Make matching results

    Args:
        fig_input (_type_): _description_
        fig_type (_type_): Zoomed or Full
    """
    if fig_input.shape[0] == 0:
        return
    else:
        if fig_input.shape[0] > 5:
            fig = plt.figure(figsize=(8,fig_input.shape[0]*.7))
        else:
            fig = plt.figure(figsize=(8,2))
        outer = fig.add_gridspec(3,1, wspace=0.0, hspace=0.0, bottom=0.05,
                                top=0.95, height_ratios=[1/fig_input.shape[0], 2/fig_input.shape[0], 
                                                            fig_input.shape[0]/fig_input.shape[0]+1])
        
        for k,i in zip(range(0,(fig_input.shape[0]+1)),fig_input.index.tolist()): # For each splicing - tx pair
            ## Prepare figure params
            tx = fig_input.loc[i]["best match tx"]
            eventid = fig_input.loc[i]["ID"]
            st_ex = fig_input.loc[i]["start_exon"]
            event = fig_input.loc[i]["event"]
            zoomed_st = int(st_ex.replace("exon",""))

            ## Make an input for visualization
            if k == 0:  # First figure with splicing event
                inner = outer[0,0].subgridspec(1, 1, wspace=0, hspace=0.0)
                axs = inner.subplots()
                ref_bed, scale = Make_ref(tx,st_ex,"scaling")
                if len(scale) > 0:
                    query_bed, error = Make_query(fig_input.loc[i], scale)
                else:
                    query_bed, error = Make_query(fig_input.loc[i])
                
                Comp_plot(axs, query_bed, "#73ABFF", ["nt"])
                plt.text(0.04,0.8,"{} ({})".format(eventid, event), 
                            transform=axs.transAxes,
                            fontsize=11)
                sns.despine(left=True, top=True, bottom=True, ax=axs)
                
                inner = outer[2,0].subgridspec(fig_input.shape[0], 1, wspace=0, hspace=0.2)
                axs = inner.subplots()

                def Tx_viewer(ax, zoom):
                    """ Draw merged figure using all transcript that has certain splicing event

                    Args:
                        ax (_type_): figure axs automatically set
                        zoom (_type_): zoom level automatically set
                    """
                    tx_block = ref_bed[zoomed_st-1:zoomed_st+4]
                    if zoom == "Zoomed":
                        Comp_plot(ax, tx_block, "#FF7878", ["nt"])
                        Comp_plot(ax, query_bed, "#73ABFF", ["nt"]) # For marking splicing event on the matching tx
                    else:
                        Comp_plot(ax, tx_block, "#FF7878", ["nt"])
                        Comp_plot(ax, ref_bed, "#FF7878", ["nt"])   # Full length transcript
                        Comp_plot(ax, query_bed, "#73ABFF", ["nt"]) # For marking splicing event on the matching tx
                    plt.text(0.04,0.8,"{} ({})".format(tx+"-"+st_ex,zoom), 
                                transform=ax.transAxes,
                                fontsize=11)
                    sns.despine(left=True, top=True, bottom=True, ax=ax)
                
                if fig_input.shape[0] != 1:
                    Tx_viewer(axs[k],fig_type)
                else:
                    Tx_viewer(axs,fig_type)
                
            else:
                ref_bed, scale = Make_ref(tx,st_ex)
                query_bed, error = Make_query(fig_input.loc[i])
                Tx_viewer(axs[k],fig_type)
    plt.savefig(OUT+"{}_{}_mapping.pdf".format(event, gene), bbox_inches="tight")    

# %%
OUT = DIR+"figure/"
for id in test_list:    # For each splicing event
    sub = query[query["ID"]==id]
    sub = sub.sort_values(by="E2s")
    if len(sub["event"].unique()) > 1:
        fig_input = sub[sub["event"]==sub["event"].unique()[0]]
        Draw_splicemap(fig_input,fig_type)
        fig_input = sub[sub["event"]==sub["event"].unique()[1]]
        Draw_splicemap(fig_input,fig_type)
    else:
        Draw_splicemap(sub,fig_type)

    # break
# %%
