#!/bin/bash
input=$1
cd ${input}

mkdir -p temp
mv *bed temp/
mv *gtf temp/
mv *txt temp/
mv tx_gene_dict temp/
mv get_Fasta/ temp/
mv cpat/ temp/
mv tpm/ temp/
mv rmat.csv temp/

echo "All intermediate files were moved to ${input}temp"
