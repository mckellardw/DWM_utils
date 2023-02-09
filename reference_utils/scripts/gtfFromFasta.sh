#!/usr/bin/bash

# Script to build a .gtf file from a .fasta, where each entry is 
# a `transcript` spanning the entire sequence

# Usage:
#   bash gtfFromFasta /path/to/transcriptome.fa /path/to/transcripts.gtf

IN_FA=$1
OUT_GTF=$2

cat ${IN_FA} | \
 awk -v RS="\n>" -v FS="|" '{print $2 "\tCUSTOM\ttranscript\t1\t" length($9) "\t.\t+\t.\tgene_id " $2 "; transcript_id " $1 "; gene_type " $8 "; gene_name " $6 "; transcript_type " $8 "; transcript_name " $5 ";\n"}' > \
 ${OUT_GTF}