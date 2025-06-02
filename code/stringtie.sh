#!/bin/bash
#SBATCH -J Run_stringtie
#SBATCH -c 40
#SBATCH --mem=40G
#SBATCH -o /home/kangh/lab-server/Tool/Splice-decoder/log/%u-%x-%j
#SBATCH -e /home/kangh/lab-server/Tool/Splice-decoder/log/%u-%x-%j.err
#SBATCH --time=72:00:00

## Load conda
source /projects/anczukow-lab/kangh/miniforge-pypy3/bin/activate base
conda activate splice-decoder
gtf=$1
bam_list=$2	# It should have full path
out=$3
mkdir -p ${out}

cat ${bam_list} | while read file
do
	name=`echo $file | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.'`
	stringtie \
		-e -B -p 8 \
		-G ${gtf} \
		-o ${out}/${name}.tab \
		${file}
done
