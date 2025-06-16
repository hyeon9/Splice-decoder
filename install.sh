#!/bin/bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create --name testenv -f yml -k &&
conda activate testenv &&
pip install -r requirements_pip.txt
