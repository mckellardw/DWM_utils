#!/usr/bin/bash

NCORES=30

MM_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/GRCm39.genome.fa"
MM_GENES="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf"

TRANS_FASTA="/workdir/dwm269/genomes/custom_sequences/FLExDUX4/FLExDUX4.fa"
TRANS_GTF="/workdir/dwm269/genomes/custom_sequences/FLExDUX4/FLExDUX4.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM28_FLExDUX4"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${MM_FASTA} ${TRANS_FASTA} \
--sjdbGTFfile ${MM_GENES} ${TRANS_GTF} \
--sjdbGTFfeatureExon exon
