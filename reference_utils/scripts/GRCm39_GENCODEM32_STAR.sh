#!/usr/bin/bash

NCORES=12
FASTA_GENOME="/workdir/dwm269/genomes/mm39_all/GENCODE_M32/GRCm39.genome.fa"
FASTA_cDNA="/workdir/dwm269/genomes/mm39_all/GENCODE_M32/gencode.vM32.transcripts.fa.gz"
GENES_DIR="/workdir/dwm269/genomes/mm39_all/GENCODE_M32/gencode.vM32.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM32"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${FASTA_GENOME} \
--sjdbGTFfile ${GENES_DIR} \
--sjdbGTFfeatureExon exon
