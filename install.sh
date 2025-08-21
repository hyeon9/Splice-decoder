#!/bin/bash

yml=pkg.yml
source "$(conda info --base)/etc/profile.d/conda.sh"
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --remove channels defaults
conda env create --name splice-decoder -f ${yml} &&
conda activate splice-decoder &&
pip install -r requirements_pip.txt
unzip toy_data.zip

## Prepare data source
unzip dat.zip
cd dat/
cd reference/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
gunzip -c GRCh38.primary_assembly.genome.fa.gz > hg38.fa 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
gunzip -c GRCm38.primary_assembly.genome.fa.gz > GRCm38.fa
