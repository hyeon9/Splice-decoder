#!/bin/bash

echo "To make right config file, you should give right answer for each question!"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate test_splice-decoder
conda_path=`conda info --envs | grep '*' | awk '{print $NF}'`
export PATH="`pwd`:$PATH"
config_maker.py ${conda_path}
