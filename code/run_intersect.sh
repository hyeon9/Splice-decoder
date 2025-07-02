#!/bin/bash
config_file=$1
splicing_type=$2
source ${config_file}

if [ ${species} == "human" ]
then
	pfam=${cpatdb}"/domain_info/small_modified_uniprot.bed"
else
        pfam=${cpatdb}"/domain_info/mouse_modified_uniprot.bed"
fi

mkdir -p ${input}/table

bedtools intersect -a ${pfam} -b ${input}/merged_w_${splicing_type}.bed -wb > ${input}/table/${splicing_type}_w_Pfam.txt
bedtools intersect -a ${pfam} -b ${input}/merged_wo_${splicing_type}.bed -wb > ${input}/table/${splicing_type}_wo_Pfam.txt
