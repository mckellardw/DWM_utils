#!/usr/bin/bash

NCORES=10
FASTA_DIR="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/GRCm39.genome.fa"
GENES_DIR="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf"
OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM28"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${FASTA_DIR} \
--sjdbGTFfile ${GENES_DIR} \
--sjdbGTFfeatureExon exon
