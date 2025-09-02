#!/bin/bash
input=$1
cd ${input}

## Make temp folder AND post_input folder
mkdir -p temp
mkdir -p post_input

## Clean up working directory
mv *bed temp/
mv *gtf temp/
mv *txt temp/
mv tx_gene_dict temp/
mv get_Fasta/ temp/
mv cpat/ temp/
mv tpm/ temp/
mv rmat.csv temp/

## Move some files to post_input
mv ./temp/merged*bed post_input/

echo "All intermediate files were moved to ${input}temp"
