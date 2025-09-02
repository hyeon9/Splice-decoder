#!/bin/bash
config_file=$1
source ${config_file}

if [[ ${species} == "human" ]]
then
	pfam=${cpatdb}"/domain_info/small_modified_uniprot.bed"
else
        pfam=${cpatdb}"/domain_info/mm10_mouse_modified_uniprot.bed"
fi

mkdir -p ${input}/table

bedtools intersect -a ${pfam} -b ${input}/merged.bed -wb > ${input}/table/Pfam.txt
