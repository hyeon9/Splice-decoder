#!/bin/bash
yml=$1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env create --name splice-decoder -f ${yml} -k &&
conda activate splice-decoder &&
pip install -r requirements_pip.txt
