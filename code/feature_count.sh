#!/bin/bash
#SBATCH -J Run_featureC
#SBATCH -c 40
#SBATCH --mem=60G
#SBATCH -o ./%u-%x-%j
#SBATCH -e ./%u-%x-%j.err
#SBATCH --time=72:00:00

## Load conda
source /projects/anczukow-lab/kangh/miniforge-pypy3/bin/activate base
conda activate splice-decoder

GTF="/flashscratch/horvam/for.Hyeongu/genomic.gtf"
BAM_DIR="/projects/palucka-lab/horvam/Splicing/STAR/"
OUTPUT="transcript_counts.txt"
THREADS=40

featureCounts -T $THREADS -p -t exon -g transcript_id -a $GTF -o $OUTPUT $BAM_DIR/*.bam
