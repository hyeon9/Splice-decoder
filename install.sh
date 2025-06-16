#!/bin/bash
yml=$1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env create --name testenv -f ${yml} -k &&
conda activate testenv &&
pip install -r requirements_pip.txt
