#!/usr/bin/bash

NCORES=30

MM_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/GRCm39.genome.fa"
MM_GENES="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf"

Reo_FASTA="/workdir/dwm269/genomes/reo_T1L.fasta"
Reo_GENES="/workdir/dwm269/genomes/reo_T1L_antisense.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM28_REOT1L-as"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${MM_FASTA} ${Reo_FASTA} \
--sjdbGTFfile ${MM_GENES} ${Reo_GENES} \
--sjdbGTFfeatureExon exon

cat ${MM_GENES} ${Reo_GENES} >> ${OUTDIR}/GRCm39_ReoT1L_merged_genes.gtf
