#!/bin/bash
yml=$1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env create --name splice-decoder -f ${yml} -k &&
conda activate splice-decoder &&
pip install -r requirements_pip.txt

## Prepara data source
# unzip dat.zip
# cd /dat
# cd reference
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
# gunzip GRCh38.primary_assembly.genome.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
# gunzip GRCm38.primary_assembly.genome.fa.gz
