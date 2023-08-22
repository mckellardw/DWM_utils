#!/bin/bash

# Example usage:
# bash script.sh mature.fa "Mus musculus" mmu bt2 2> mmu/bt2.log


# Input parameters
FASTA_FILE=$1
SPECIES=$2
OUTDIR=$3
PREFIX=$4

mkdir -p $OUTDIR

# Replace spaces in species name with underscores
SPECIES_NO_SPACES=$(echo $SPECIES | tr ' ' '_')

# Filter fasta file for specific species, and change Us to Ts for DNA seq data
grep -A 1 "$SPECIES" $FASTA_FILE \
| awk '{if(NR%2==0) gsub("U","T"); print}' \
| awk '/^>/{sub(/.*/, ">"$NF)}1' \
> $OUTDIR/${SPECIES_NO_SPACES}_DNA.fa

# Generate bowtie2 index
bowtie2-build $OUTDIR/${SPECIES_NO_SPACES}_DNA.fa $OUTDIR/$PREFIX
