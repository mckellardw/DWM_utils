#!/bin/bash

# THis is a simple script that can be used to download .fastq files from SRA. Please note, the files must have originally been uploaded as fastqs
# Dependencies:
#   SRA toolkit: https://hpc.nih.gov/apps/sratoolkit.html (also available via conda)
#   parallel-fastq-dump: https://github.com/rvalieris/parallel-fastq-dump (also available via conda)

# Directory where you want to store the fastqs
DATADIR="/workdir/dwm269/path/to/data/"
cd ${DATADIR}

# SRR ID list- can be directly downloaded from the SRA Run Selector
SAL="SRR_Acc_List.txt"

# Number of threads you want it to use
NTHREADS=6

# load in SRR IDs
readarray -t SRR < ${SAL}

# Loop - prefetch and then convert
for i in "${SRR[@]}"
do
  echo ${i}
  prefetch \
  --verify yes \
  --max-size 9999999999 \
  --output-directory ${DATADIR} \
  ${i}

  # My preferred tool for SRA -> fastq conversion
  parallel-fastq-dump \
  --sra-id ${i}.sra \
  --threads ${NTHREADS} \
  --outdir ${DATADIR} \
  --tmpdir ${DATADIR} \
  --split-files \
  --gzip && \
  rm ${i}.sra
done
