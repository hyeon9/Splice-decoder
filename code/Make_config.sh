#!/bin/bash

echo "To make right config file, you should give right answer for each question!"
source /projects/anczukow-lab/kangh/miniforge-pypy3/bin/activate base
conda activate splice-decoder
export PATH="`pwd`:$PATH"
config_maker.py
