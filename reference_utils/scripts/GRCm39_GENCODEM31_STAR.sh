#!/usr/bin/bash

NCORES=24
FASTA_GENOME="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/GRCm39.genome.fa"
FASTA_cDNA="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/gencode.vM31.transcripts.fa.gz"
GENES_DIR="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/gencode.vM31.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM31"
FASTA_rRNA=${OUTDIR}/STAR_GRCm39_GENCODEM31_rRNA.fa

mkdir -p ${OUTDIR}
cd ${OUTDIR}

zcat ${FASTA_cDNA} | awk -v RS="\n>" -v FS="|" '$8=="rRNA" {print ">"$0}' > ${FASTA_rRNA}

STAR \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${FASTA_DIR} \
--sjdbGTFfile ${GENES_DIR} \
--sjdbGTFfeatureExon exon
