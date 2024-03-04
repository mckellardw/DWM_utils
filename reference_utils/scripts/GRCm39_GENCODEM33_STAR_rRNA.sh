#!/usr/bin/bash

# Script to extract transcript sequences and annotations for rRNAs, and generate a STAR reference

NCORES=48
# FASTA_GENOME="/gpfs/commons/groups/innovation/dwm/genomes/mus_musculus/gencode/release_M33/GRCm39.genome.fa"
FASTA_cDNA="/gpfs/commons/groups/innovation/dwm/genomes/mus_musculus/gencode/release_M33/gencode.vM33.transcripts.fa.gz"
# GENES="/gpfs/commons/groups/innovation/dwm/genomes/mus_musculus/gencode/release_M33/gencode.vM33.annotation.gtf"

OUTDIR="/gpfs/commons/groups/innovation/dwm/genomes/mus_musculus/STAR_GRCm39_GENCODEM33_rRNA"
FASTA_rRNA=${OUTDIR}/GRCm39_GENCODEM33_rRNA.fa
GENES_rRNA=${OUTDIR}/GRCm39_GENCODEM33_rRNA.gtf

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# Extract rRNA sequences from cDNA/transcript fasta
zcat ${FASTA_cDNA} | \
 awk -v RS="\n>" -v FS="|" '$8=="rRNA" || $8=="Mt_rRNA" { print ">"$2 $9 }' > \
 ${FASTA_rRNA}

# Build custom gtf from the cDNA/transcript fasta
zcat ${FASTA_cDNA} | \
 awk -v RS="\n>" -v FS="|" '$8=="rRNA" || $8=="Mt_rRNA" {print $2 "\tCUSTOM\texon\t1\t" gsub(/A/,"",$9)+gsub(/C/,"",$9)+gsub(/G/,"",$9)+gsub(/T/,"",$9) "\t.\t+\t.\tgene_id " $2 "; transcript_id " $1 "; gene_type " $8 "; gene_name " $6 "; transcript_type " $8 "; transcript_name " $5 ";\n"}' > \
 ${GENES_rRNA}

# Can't build from the genomic .gtf b/c chromoseome names & start/stop positions are all mismatched
# cat ${GENES} | \
#  grep -E 'gene_type "rRNA"|gene_type "Mt_rRNA"' > \
#  awk -v FS="\t" '#TODO' \
#  ${GENES_rRNA}

STAR \
    --runThreadN ${NCORES} \
    --runMode genomeGenerate \
    --genomeDir ${OUTDIR} \
    --genomeFastaFiles ${FASTA_rRNA} \
    --sjdbGTFfile ${GENES_rRNA} \
    --genomeSAindexNbases 6 \
    --sjdbGTFfeatureExon exon
