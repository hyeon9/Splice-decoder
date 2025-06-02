# %%
import pandas as pd
DIR = "/home/kangh/lab-server/Tool/Splice-decoder/example/MYC_8h/result/PE_val/"
data = pd.read_csv(f"{DIR}8h_unmat_pairs_NMD_info.txt",
                   sep="\t",
                   header=None)
data["pair"] = data[0]+data[1]

k = 0
for i in data["pair"].unique():
    temp = data[data["pair"]==i]
    temp["NMD"] = temp[2].str.split("|").str[1]
    temp["NMD"] = temp["NMD"].copy().apply(lambda x : 1 if x == "NMD" else 0)
    if temp["NMD"].sum() > 1:
        k += 1

print(k/len(data["pair"].unique()))
# %%
