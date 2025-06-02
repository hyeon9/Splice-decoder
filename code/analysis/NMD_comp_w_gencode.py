# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import HGfig

# %%
def Make_pie():
    
    data = pd.read_csv("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/pie_input.txt",
                       sep="\t",
                       header=None)
    data = data[data[0]!="lncRNA"]
    data = pd.DataFrame(data[0].value_counts())
    
    plt.figure(figsize=(4,4))
    def func(pct, allvals):
        print(allvals)
        absolute = int(np.round(pct/100.*np.sum(allvals)))
        return f"{pct:.1f}%\n({absolute:d})"


    wedges, texts, autotexts = plt.pie(data["count"], autopct=lambda pct: func(pct, data["count"]),
                                       textprops=dict(color="w"),
                                       colors=plt.cm.Set1.colors)
    plt.legend(wedges, data["count"].index,
               loc="center left",
               bbox_to_anchor=(1, 0, 0.5, 1),
               frameon=False)
    plt.setp(autotexts, size=11, weight="bold")
    plt.savefig("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/figure/pie.pdf",
                bbox_inches="tight")

    plt.figure(figsize=(4,4))
    def func(pct, allvals):
        print(allvals)
        absolute = int(np.round(pct/100.*np.sum(allvals)))
        return f"{pct:.1f}%\n({absolute:d})"


    wedges, texts, autotexts = plt.pie([145,6], autopct=lambda pct: func(pct, [145,6]),
                                       textprops=dict(color="w"),
                                       colors=["#377EB8","#A8A8A8FF"])
    plt.setp(autotexts, size=11, weight="bold")
    plt.title("Retained intron")
    

Make_pie()
# %%
def pre_strategy():
    """
    Make GTF from the targeted bed file
    """
    
    DIR="/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/"
    data = pd.read_csv(f"{DIR}/targeted_tx_250317.bed",
                       sep="\t",
                       header=None) # Check usage.txt (TX w/ NMD causing DS)
    data[8] = data[4]+"-"+data[6]+"-"+data[7]

    with open(f"/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/target_tx_250317.gtf", "w") as out:
        for tx in data[8].unique():
            ind = data[data[8]==tx]
            for i in range(ind.shape[0]):
                tx = tx.replace(";","-")    # To avoid error in gffcompare
                if i == 0:  # +1 to start (it was minus 1 to get fasta sequences)
                    tx_line = f'{ind[0].unique()[0]}\tSpliceDecoder\ttranscript\t{ind[1].min()+1}\t{ind[2].max()}\t.\t{ind[5].unique()[0]}\t.\ttranscript_id "{tx}";'
                    exon_line = f'{ind.iloc[i][0]}\tSpliceDecoder\texon\t{ind.iloc[i][1]+1}\t{ind.iloc[i][2]}\t.\t{ind.iloc[i][5]}\t.\ttranscript_id "{tx}";'
                    out.write(f"{tx_line}\n{exon_line}\n")
                else:
                    exon_line = f'{ind.iloc[i][0]}\tSpliceDecoder\texon\t{ind.iloc[i][1]+1}\t{ind.iloc[i][2]}\t.\t{ind.iloc[i][5]}\t.\ttranscript_id "{tx}";'
                    out.write(f"{exon_line}\n")
                    
pre_strategy()
# %%
def Decoding():
    """ 
    Exchange relative position to genomic position (From SpliceDecoder output)
    """
    
    DIR = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/"
    data = pd.read_csv(f"{DIR}/temp.txt",
                       sep="\t",
                       header=None) # Used Start site within SpliceDecoder
    gtf = pd.read_csv(f"{DIR}/target_tx.gtf",
                      sep="\t",
                      header=None)  # SpliceDecoder input GTF
    gtf = gtf[gtf[2]=="exon"]
    out = []
    for i in data[0].unique():
        temp = gtf[gtf[8]==f'transcript_id "{i.replace(";","")}";'][[3,4,6]]
        temp["len"] = abs(temp[4]-temp[3])+1
        rel = data[data[0]==i][1].values[0] # rel position
        for e in range(temp.shape[0]):  # exon length
            if rel >= temp.iloc[e,3]:   # Compare relative start codon and cumulative exon length
                rel -= temp.iloc[e,3]   # Deduct exon length from rel start codon
            else:
                if temp[6].unique()[0] == "+":
                    g_coor = temp.iloc[e,0]+rel+1 # Start + rest rel length and +1 to adjust rel length
                else:
                    g_coor = temp.iloc[e,1]-rel+1 # Start - rest rel length and +1 to adjust rel length
                break
        out.append([i,g_coor])
    
    pd.DataFrame(out).to_csv(f"{DIR}/pNMD_with_genome_start.txt",
                             sep="\t",
                             header=None,
                             index=None)    # SpliceDecoder simTX w/ genomic pos of start codon
    

Decoding()

# %%
def Make_CDS():
    """ 
    Make CDS bed from SpliceDecoder output
    """
    
    DIR = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/"
    data = pd.read_csv(f"{DIR}/temp_w_stop.txt",
                       sep="\t",
                       header=None) # Used Start site within SpliceDecoder
    gtf = pd.read_csv(f"{DIR}/target_tx.gtf",
                      sep="\t",
                      header=None)  # SpliceDecoder input GTF
    gtf = gtf[gtf[2]=="exon"]
    out = []
    for i in data[0].unique():
        temp = gtf[gtf[8]==f'transcript_id "{i.replace(";","")}";'][[3,4,6]]
        temp["len"] = abs(temp[4]-temp[3])+1
        rel = data[data[0]==i][1].values[0] # rel position
        stop = data[data[0]==i][2].values[0] # rel position
        action = 0
        for e in range(temp.shape[0]):  # exon length
            if rel >= temp.iloc[e,3]:   # Compare relative start codon and cumulative exon length
                rel -= temp.iloc[e,3]   # Deduct exon length from rel start codon
            elif rel < temp.iloc[e,3] and action == 0:
                action += 1
                if temp[6].unique()[0] == "+":
                    g_coor = temp.iloc[e,0]+rel-2 # Smaller
                    g_coor2 = temp.iloc[e,1]    # Larger
                else:
                    g_coor2 = temp.iloc[e,1]-rel+1   # Larger
                    g_coor = temp.iloc[e,0]-1 # Samller
                
                out.append([i.split("-")[4],int(g_coor),int(g_coor2),temp[6].unique()[0],i])  # Add CDS
                rel = 0

            elif rel < temp.iloc[e,3] and action > 0:
                if temp[6].unique()[0] == "+":
                    g_coor = temp.iloc[e,0]-1 # Smaller
                    g_coor2 = temp.iloc[e,1]    # Larger
                else:
                    g_coor2 = temp.iloc[e,1]    # Larger
                    g_coor = temp.iloc[e,0]-1 # Samller
                
                out.append([i.split("-")[4],int(g_coor),int(g_coor2),temp[6].unique()[0],i])  # Add CDS
                rel = 0
            
            if stop >= temp.iloc[e,3]:
                stop -= temp.iloc[e,3]
            else:
                if temp[6].unique()[0] == "+":
                    g_coor = temp.iloc[e,0]-1 # Smaller
                    g_coor2 = temp.iloc[e,0]+stop-2 # Larger
                else:
                    g_coor = temp.iloc[e,1]-stop+1   # Samller
                    g_coor2 = temp.iloc[e,1]-1 # Larger
                
                out.append([i.split("-")[4],int(g_coor),int(g_coor2),temp[6].unique()[0],i])  # Add CDS
                break
    
    pd.DataFrame(out).to_csv(f"{DIR}/pNMD_CDS.bed",
                             sep="\t",
                             header=None,
                             index=None)    # SpliceDecoder simTX w/ genomic pos of start codon

Make_CDS()
# %%
## Compare of translated AA sequences between protine coding tx and matched simulated transcripts
DIR='/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/'
pair = pd.read_csv(f"{DIR}/simtx_w_protein_coding_category.txt",
                   sep="\t",
                   header=None)
pair[5] = pair[5].str.split("|").str[0]
target = pd.read_csv(f"{DIR}/PC_tx/70_txs.txt",
                     sep="\t",
                     header=None)
pair = pair[pair[3].isin(target[0])]
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2
from Bio.Seq import Seq

mat=0
unmat=0
align_score = []
unmat_list = []
for pc,sim in zip(pair[3],pair[5]):
    faa = "".join(str(record.seq) for record in SeqIO.parse(f"{DIR}/PC_tx/{pc}.fa", "fasta"))
    fab = "".join(str(record.seq) for record in SeqIO.parse(f"{DIR}/PC_tx/{sim}.fa", "fasta"))
    
    ## Translate and pairwise align
    from Bio.Seq import Seq
    fa_aa = (Seq(faa).translate())
    fb_aa = (Seq(fab).translate())
    alignments = pairwise2.align.globalxx(fa_aa, fb_aa)
    min_len = min([len(fa_aa), len(fb_aa)])
    for alignment in alignments[:1]:
        align_score.append((float(alignments[0].score)/float(min_len)))
        if (float(alignments[0].score)/float(min_len)) > 0.99:  # According to my observation, 100% cases are against to 50-55 rules
            mat +=1
            pass
            # print(pc,sim)
            # print(pairwise2.format_alignment(*alignment))
        else:
        # elif (float(alignments[0].score)/float(min_len)) > 0.75 and (float(alignments[0].score)/float(min_len))< 1:
            unmat +=1
            unmat_list.append([pc,sim])
            # print(pc,"\n",sim,"\n",[len(fa_aa), len(fb_aa)])
            # print(pairwise2.format_alignment(*alignment))
            break

unmat_list_df = pd.DataFrame(unmat_list)
unmat_list_df.to_csv(DIR+"PC_tx/Diff_frame_55.txt",
                     sep="\t",
                     index=None,
                     header=None)
# %%
plt.figure(figsize=(5,10))
sns.histplot(align_score,
             bins=50)
sns.despine()
plt.xlabel("Alignment score\n(Normalized input length)")
# %%
print(mat,unmat)
# %%
def genome_start(category):
    """
    Extract genomic coordinates and compare gencode tx CDS
    to explain transcripts with protein coding tags
    """
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    

    dict_key = {"PC":["query_tx.gtf","simtx_w_protein_coding_category.txt"],
                "NMD":["NMD_query_tx.gtf","simtx_w_NMD_category.txt"]}
    DIR = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/"
    out = []
    with open(f"{DIR}/{dict_key[category][0]}") as tx_st_file:
        for line in tx_st_file:
            read_line = line.strip().split("\t")
            if read_line[2] == "transcript":
                hit = 0
                for ele in read_line[8].strip().split(";"):
                    if ele.startswith(" transcript_id"):
                        tx_id = ele.split("\"")[1]
            elif read_line[2] == "CDS" and hit == 0:
                if read_line[6] == "+":
                    out.append([tx_id,read_line[3]])
                else:
                    out.append([tx_id,read_line[4]])
                hit += 1
    
    out = pd.DataFrame(out)
    out = out.drop_duplicates()
    out.columns = ["tx","ref_genome_start"]
    
    mat_file = pd.read_csv(f"{DIR}/{dict_key[category][1]}",
                           sep="\t",
                           header=None)
    mat_file[5] = mat_file[5].str.split("|").str[0]
    pnmd = pd.read_csv("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/pNMD_with_genome_start.txt",
                       sep="\t",
                       header=None)
    pnmd.columns = ["simtx","sim_genome_start"]
    pnmd["simtx"] = pnmd["simtx"].str.replace(";","")
    mat_file = pd.merge(mat_file, out,
                        left_on=3,
                        right_on="tx")
    mat_file = pd.merge(mat_file, pnmd,
                        left_on=5,
                        right_on="simtx")
    mat_file[5] = mat_file[5].str.replace("---","-=-")
    mat_file["event"] = mat_file[5].str.split("-").str[1]
    mat_file["strand"] = mat_file[5].str.split("-").str[5].replace("=","-")
    mat_file[5] = mat_file[5].str.split("-").str[6:]
    mat_file = mat_file.drop(["simtx"],axis=1)
    # mat_file = mat_file[mat_file["event"].isin(["EI","ES","Alt_A3SS","Alt_A5SS","Ori_A5SS","Ori_A3SS","RI","SI","MXE1","MXE2"])]
    mat_file["Exon_start"] = mat_file[5].str[2]
    mat_file["Exon_end"] = mat_file[5].str[3]
    p_mat = mat_file[mat_file["strand"]=="+"]
    n_mat = mat_file[(mat_file["strand"]=="-") | (mat_file["strand"]=="chr16")]   # To save BOLA2-SMG1P6

    # mat_file["Start Codon Diff"] = (mat_file["ref_genome_start"].astype(float) - mat_file["sim_genome_start"].astype(float)).apply(lambda x : "Matched" if x == 0 else "Unmatched")
    # plt.figure(figsize=(4,1))
    # sns.countplot(mat_file["Start Codon Diff"])
    mat_file["Start Codon Diff"] = (mat_file["ref_genome_start"].astype(float) - mat_file["sim_genome_start"].astype(float))
    print(mat_file.head(10))
    
    print(f"Number of downstream start (it may avoid the PTC effect)")    # Downstream of NMD asosc DS
    print(f'ref_TX_on + strand: {p_mat[p_mat["Exon_end"].astype(float) < p_mat["ref_genome_start"].astype(float)].shape[0] + n_mat[n_mat["Exon_end"].astype(float) > n_mat["ref_genome_start"].astype(float)].shape[0]}')
    print(f'sim_TX_on + strand: {p_mat[p_mat["Exon_end"].astype(float) < p_mat["sim_genome_start"].astype(float)].shape[0] + n_mat[n_mat["Exon_end"].astype(float) > n_mat["sim_genome_start"].astype(float)].shape[0]}\n')
    # print(f'ref_TX_on - strand: {n_mat[n_mat["Exon_end"].astype(float) > n_mat["ref_genome_start"].astype(float)].shape[0]}')
    # print(f'sim_TX_on + strand: {n_mat[n_mat["Exon_end"].astype(float) > n_mat["sim_genome_start"].astype(float)].shape[0]}\n')
    print(f"Number of upstream start (it should be affected by PTC)")   # Upstream of NMD asosc DS
    print(f'ref_TX_on + strand: {p_mat[p_mat["Exon_end"].astype(float) >= p_mat["ref_genome_start"].astype(float)].shape[0] + n_mat[n_mat["Exon_end"].astype(float) <= n_mat["ref_genome_start"].astype(float)].shape[0]}')
    print(f'sim_TX_on + strand: {p_mat[p_mat["Exon_end"].astype(float) >= p_mat["sim_genome_start"].astype(float)].shape[0] + n_mat[n_mat["Exon_end"].astype(float) <= n_mat["sim_genome_start"].astype(float)].shape[0]}\n')

    ## To investigate upstream PC genes
    query_list = p_mat[p_mat["Exon_end"].astype(float) >= p_mat["ref_genome_start"].astype(float)]
    query_list2 = n_mat[n_mat["Exon_end"].astype(float) <= n_mat["ref_genome_start"].astype(float)]
    pd.DataFrame(query_list[3].tolist() + query_list2[3].tolist()).to_csv('/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/PC_tx/70_txs.txt',
                                                                          sep="\t",
                                                                          index=None,
                                                                          header=None)
    ## Downstream PC genes
    query_list = p_mat[p_mat["Exon_end"].astype(float) < p_mat["ref_genome_start"].astype(float)]
    query_list2 = n_mat[n_mat["Exon_end"].astype(float) > n_mat["ref_genome_start"].astype(float)]
    pd.DataFrame(query_list[3].tolist() + query_list2[3].tolist()).to_csv('/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/PC_tx/47_txs.txt',
                                                                          sep="\t",
                                                                          index=None,
                                                                          header=None)

import matplotlib.pyplot as plt
genome_start("PC")
# genome_start("NMD")
# plt.legend(frameon=False)
# %%
def Apply_frame():
    DIR="/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/"
    query = pd.read_csv(f"{DIR}/NMD_val/PC_tx/47_txs.txt",
                        sep="\t",
                        header=None)
    query2 = pd.read_csv(f"{DIR}/NMD_val/PC_tx/Diff_frame_55.txt",
                         sep="\t",
                         header=None)
    query_tx = query[0].tolist() + query2[0].tolist()
    sim_frame = pd.read_csv(f"{DIR}/pNMD_with_genome_start.txt",
                            sep="\t",
                            header=None)
    pc_tx = pd.read_csv(f"{DIR}/NMD_val/simtx_w_protein_coding_category.txt",
                        sep="\t",
                        header=None)
    pc_tx[5] = pc_tx[5].str.split("|").str[0]
    pc_tx = pc_tx[pc_tx[3].isin(query_tx)]
    pc_gtf = pd.read_csv(f"{DIR}/NMD_val/query_tx.gtf",
                        sep="\t",
                        header=None)
    k = 0
    for pc,sim in zip(pc_tx[3],pc_tx[5]):
        sim_start = sim_frame[sim_frame[0].str.startswith(sim)][1].astype(int).values[0]
        gtf = pc_gtf[(pc_gtf[8].str.contains(pc)) & (pc_gtf[2]!="transcript")].copy()
        cds = gtf[gtf[2]=="CDS"]
        stop = cds.iloc[[-1]]
        start_cds = gtf[(gtf[3] <= sim_start) & (gtf[4] >= sim_start)]
        if start_cds.shape[0] > 0:
            if start_cds[6].unique() == "+":
                start_cds = start_cds[start_cds[2]=="exon"]
                start_cds[3] = sim_start-2  # To get fasta sequence it requires -1 position
                others = gtf.loc[start_cds.index[0]+1:,:]
                others = others[others[2]=="exon"]
                others = others.loc[:stop.index[0],:]   # limit by the exon that has stop codon
                target = pd.concat([start_cds,others])
            else:   # For minus strand
                start_cds = start_cds[start_cds[2]=="exon"]
                start_cds[4] = sim_start
                others = gtf.loc[start_cds.index[0]+1:,:]
                others = others[others[2]=="exon"]
                others = others.loc[:stop.index[0],:]   # limit by the exon that has stop codon
                target = pd.concat([start_cds,others])
            target.iloc[:,:8].to_csv(f"{DIR}/NMD_val/PC_tx/frame_apply/{pc}.bed",
                                    sep="\t",
                                    index=None,
                                    header=None)
        else:
            # print(gtf.iloc[:,:7], sim_start)
            print(pc,sim)
            k += 1
            # break
    print(k)

Apply_frame()
# %%
def Stop_check():
    from Bio import SeqIO
    from Bio.Seq import Seq
    
    out = []
    DIR = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/PC_tx/frame_apply/"
    DIR2 = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/PC_tx/"
    files = [f for f in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, f)) if f.endswith(".fa")]
    hit = 0
    for pc in files:
        original = pc.replace(".bed","")
        # ori_fa = "".join(str(record.seq) for record in SeqIO.parse(f"{DIR2}/{original}", "fasta"))
        # ori_aa = (Seq(ori_fa).translate())
        fa = "".join(str(record.seq) for record in SeqIO.parse(f"{DIR}/{pc}", "fasta"))
        aa = (Seq(fa).translate())
        if "*" in aa:
            # print(f"PTC: {aa.index('*')}, Total CDS length: {len(aa)}")
            out.append([aa.index('*'),len(aa)])
            hit += 1

    df = pd.DataFrame(out)
    df.columns = ["PTC","CDS"]
    df["Coverage of\noriginal CDS (%)"] = (df["PTC"] / df["CDS"]) * 100
    plt.figure(figsize=(1.5,3.5))
    sns.violinplot(df['Coverage of\noriginal CDS (%)'])
    sns.despine()
    plt.savefig("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/NMD_val/figure/PTC_check.pdf",
                bbox_inches="tight")
    # plt.xlim(0,1)
    # print(hit, len(files))
        

Stop_check()
# %%
def rel_start():
    """
    Calculate relative position on the transcript (Ref data)
    """
    import pandas as pd
    DIR = "/flashscratch/kangh/temp_SD/"
    out = []
    with open(f"{DIR}/query_tx.gtf") as tx_st_file:
        for line in tx_st_file:
            read_line = line.strip().split("\t")
            if read_line[2] == "transcript":
                comp = 0
                pre_tx = []
                length = []
                for ele in read_line[8].strip().split(";"):
                    if ele.startswith(" transcript_id"):
                        tx_id = ele.split("\"")[1]
                        pre_tx.append(tx_id)
            elif read_line[2] == "exon" and comp == 0:
                length.append(abs(int(read_line[3]) - int(read_line[4])))
            elif read_line[2] == "CDS" and len(length) >= 1:
                adj_length = abs(length[-1] - abs(int(read_line[3])-int(read_line[4])))
                length[-1] = adj_length
                comp += 1
                out.append([tx_id, sum(length)])
                length = []
            else:
                pass
        out = pd.DataFrame(out)
        out.to_csv(f"{DIR}/Ref_start.txt",
                   sep="\t",
                   index=None,
                   header=None)

        return out


# rel_start()
# %%
def stop_finder():
    """ Find any stop codon for the given sequences (particularly RI sequences)

    Returns:
        result: Stop codon existence (0/1) for all possible reading frame combinations
    """
    from Bio import SeqIO
    import pandas as pd
    def find_stop_codons(seq, start_pos=0):
        stop_codons = {"TAA", "TAG", "TGA"}
        positions = []
        
        for i in range(start_pos, len(seq) - 2, 3):  # Iterate in codon steps
            codon = seq[i:i+3]
            if codon in stop_codons:
                positions.append(i)
        
        return positions
        
            
    DIR="/flashscratch/kangh/temp_SD/"
    with open(f"{DIR}/RI.fa") as ri:
        result = []
        for seq in ri:
            if seq.startswith(">"):
                seq_id = seq.strip()
            else:
                hit = 0
                for frame in range(3):
                    stop_positions = find_stop_codons(seq, start_pos=frame)
                    if len(stop_positions) > 0:
                        hit += 1
                result.append([seq_id,hit])
        
    return pd.DataFrame(result)
            
result = stop_finder()
result.columns = ["seq","n_stop"]
print(result["n_stop"].value_counts())
# %%
import pandas as pd
DIR = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/"
pnmd = pd.read_csv(f"{DIR}result/Whole_DS_score_NMD.txt",
                   sep="\t")
pnmd = pnmd[pnmd["Reference_transcript"].str.contains("ENST")]
pnmd = pnmd[pnmd["Probability_of_NMD"]==-1]   # Collect DS event that cause NMD on sim TX
pnmd = pnmd.drop_duplicates()
print(pnmd.head())

# %%
mapped = pd.read_csv(f"{DIR}/fig_input.txt",
                     sep="\t")
mapped = mapped[(mapped["ID"].isin(pnmd["LongID"].tolist())) & (mapped["whole_diff"]==0)]   # Consider only perfect matching transcripts
print(mapped.head())

# %%
def Make_summ():
    ## Gencode process
    with open("/flashscratch/kangh/temp_SD/gencode.v32.primary_assembly.annotation.gtf", "r") as file:
        out = open("/flashscratch/kangh/temp_SD/tx_info.txt", "w")
        for line in file:
            if not line.startswith("#"):
                read_line = line.strip().split("\t")
                if "transcript" == read_line[2]:
                    ele = read_line[8].split(";")
                    for j in ele:
                        if j.startswith(" transcript_id"):
                            out.write(j.split("\"")[1]+"\t")
                        elif j.startswith(" transcript_type"):
                            out.write(j.split("\"")[1]+"\t")
                        elif j.startswith(" level"):
                            out.write(j.split(" ")[2])
                    out.write("\n")
Make_summ()
# %%
tx_type = pd.read_csv("/flashscratch/kangh/temp_SD/tx_info.txt",
                      sep="\t",
                      header=None)
tx_type = tx_type[tx_type[0].isin(mapped["best match tx"].tolist())]

check = []
for i,j in zip(pnmd["LongID"],pnmd["Simulated_event"]):
    target = mapped[(mapped["ID"]==i) & (mapped["event"]==j) & (mapped["whole_diff"]==0)]
    for tx in target["best match tx"].unique():
        if tx.startswith("ENST"):
            if tx in tx_type[0].tolist():   # Check DS ID
                check.append([i,j,tx,tx_type[tx_type[0]==tx][1].values[0],tx_type[tx_type[0]==tx][2].values[0]])
                # check.append([i,tx]) 
                # break   # At least one match for each DS event
    # print(i,j)
    # print(i,event_dict[j])
    # break
# %%
print(len(pnmd["LongID"].unique()))
# %%
check_df = pd.DataFrame(check)
check_df.to_csv("/flashscratch/kangh/temp_SD/result_gencode_match.txt",
                sep="\t",
                index=None)
print(check_df)
# %%
check_df = check_df[check_df[2]<4]
print(check_df[1].value_counts())
# %%
data = pd.read_csv("/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/SF_NMD.txt",
                   sep="\t",
                   header=None)
data = data[data[2]==-1]
data = data[data[0].str.startswith("CA;")]
print(data)
df = []
for i in data[1].unique():
    temp = data[data[1]==i]
    df.append([i,len(temp[0].unique())])
# %%
out = pd.DataFrame(df)
print(out)
plt.figure(figsize=(2,3))
sns.violinplot(out,
               width=.6)
plt.xticks([])
plt.ylabel("Number of NMD causing DS")
# %%
