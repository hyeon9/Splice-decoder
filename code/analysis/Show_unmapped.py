#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import sys
import os
import sys
import argparse

DIR = "/home/kangh/lab-server/Tool/Splice-decoder_git/example/Luminal/"
unmap = pd.read_csv(f"{DIR}unmapped.txt",
                    sep="\t")

# %%
nj = pd.read_csv(f"{DIR}/novelJC.csv",
                 sep=",")
# %%
print(len(list(set(unmap["ID"]) & set(nj["long_ID"]))))
print(len(set(unmap["ID"])))
# %%
