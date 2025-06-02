#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import pearsonr, mannwhitneyu
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# plt.rcParams["font.family"] = "Arial"

DIR = "/home/kangh/lab-server/Tool/Splice-decoder_git/example/MYC8h/result/"
OUT = DIR+"figure/"
if not os.path.exists(OUT):
    os.mkdir(OUT)
# doa_data = pd.read_csv(DIR+"Classified_DS_DOA.txt",
#                    sep="\t")
# doa_data["functional_event"] = doa_data["domain_integrity"]
# nmd_data = pd.read_csv(DIR+"Classified_DS_NMD.txt",
#                    sep="\t")
# nmd_data["functional_event"] = nmd_data["pNMD"]
# data = pd.concat([doa_data,nmd_data])
doa_data = pd.read_csv(DIR+"Whole_DS_score_DOA.txt",
                       sep="\t")
doa_data["functional_event"] = doa_data["Domain_integrity"]
nmd_data = pd.read_csv(DIR+"Whole_DS_score_NMD.txt",
                       sep="\t")
nmd_data["functional_event"] = nmd_data["pNMD"]
data = pd.concat([doa_data,nmd_data])

data2 = pd.read_csv("/".join(DIR.split("/")[:-2])+"/rmat.csv",
                    sep=",")
data2 = data2[(data2["FDR"] < 0.05) &
              (np.abs(data2["dPSI_2_minus_1"]) > 0.1)]

data2 = data2[["geneSymbol"]].astype(str)
data2.columns = ["gene"]
print(len(data["gene"].unique()))
print(len(data2["gene"].unique()))
# %%
def GO(target_list, cutoff, term, title):
    key_set = {"mf":"GO_Molecular_Function_2018",
               "bp":"GO_Biological_Process_2018", 
               "react":"Reactome_2022",
               "Msig":"MSigDB_Hallmark_2020",
               "cc":"GO_Cellular_Component_2018"}

    import gseapy as gp
    enr = gp.enrichr(gene_list=target_list,
                    gene_sets=key_set[term],  #GO_Biological_Process_2023, MSigDB_Hallmark_2020, GO_Molecular_Function_2023, GO_Cellular_Component_2023
                    organism="human", # don't forget to set organism to the one you desired! e.g. Yeast
                    outdir=DIR+"GO_{}".format("Merged"), # don't write to disk
                    cutoff=cutoff
                    )
    
    input = enr.results[enr.results["Adjusted P-value"]<cutoff][["Term", "Adjusted P-value","Combined Score","Genes"]]
    input["-log(q)"] = -np.log10(input["Adjusted P-value"])
    input["Term"] = input["Term"].str.split("(").str[0]
    plt.figure(figsize=(5,1.5))
    plt.scatter(y="-log(q)", x="Term", 
                data=input.sort_values(by="-log(q)", ascending=False).iloc[:20,:],
                c="Combined Score",
                s=80)
    plt.ylabel("-Log10(FDR)")
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8, rotation=90)
    sns.despine()
    plt.savefig(DIR+"figure/{}2023.pdf".format(title),
                bbox_inches="tight")

    return input, key_set[term]


terms = "Msig"
result, title = GO(data["gene"], 0.1, terms, "{}_decode".format(terms))
result2, title = GO(data2["gene"], 0.1, terms, "{}_rmat".format(terms))
# %%
def Comp_GO(a,b,title):
    ## Comp GO result between original rMATS vs filtered result
    result["method"] = a
    result2["method"] = b
    shared_term = list(set(result["Term"].tolist()) & set(result2["Term"].tolist()))
    inputa = result[["method","Term","-log(q)"]].set_index("Term")
    inputa = inputa.reindex(shared_term)
    inputa = inputa.sort_values("-log(q)",
                                ascending=False)
    
    inputb = result2[["method","Term","-log(q)"]].set_index("Term")
    inputb = inputb.reindex(shared_term)
    
    merged_df = pd.concat([inputa, inputb])

    plt.figure(figsize=(6,2))
    custom_pal = {a:"#FF7E7E",
                  b:"#9F9F9F"}
    sns.scatterplot(x="Term", y="-log(q)", data=merged_df, hue="method",
                    palette=custom_pal, s=55,
                    edgecolor=None)
    inc_terms = []
    dec_terms = []
    for term in shared_term:
        sc = merged_df[(merged_df.index==term) &
                        (merged_df["method"]==a)]["-log(q)"].values
        mat = merged_df[(merged_df.index==term) &
                        (merged_df["method"]==b)]["-log(q)"].values
        if (sc - mat) > 0.5:
            inc_terms.append(term)
        elif (mat - sc) > 0.5:
            dec_terms.append(term)

    plt.xticks(rotation=90)
    sns.despine()
    plt.legend(frameon=False,
               bbox_to_anchor=(1.02,1.02))
    plt.xlabel(title)
    plt.savefig(DIR+"figure/{}_{}_vs_{}.pdf".format(terms,a,b),
                bbox_inches="tight")
    
    print("Shared: {}, MATS_specific: {}, Decoder_specific: {}".format(len(shared_term),
                                                                       len(list(set(result2["Term"].tolist()) - set(result["Term"].tolist()))),
                                                                       len(list(set(result["Term"].tolist()) - set(result2["Term"].tolist())))
                                                                       )
                                                                       )
    return shared_term, inc_terms, dec_terms
    
shared_term, inc_terms, dec_terms = Comp_GO("Splice-decoder","rMATS", title)
# %%
