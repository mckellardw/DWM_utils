#!/usr/bin/bash

NCORES=64

MM_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/GRCm39.genome.fa"
MM_GENES="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/gencode.vM31.annotation.gtf"

Reo_FASTA="/workdir/dwm269/genomes/custom_sequences/reo_T1L.fasta"
Reo_GENES="/workdir/dwm269/genomes/custom_sequences/reo_T1L_antisense.gtf"

merged_GENES="GRCm39_M31_ReoT1L_merged_genes.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM31_REOT1L-as"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

#cat ${MM_GENES} ${Reo_GENES} >> ${merged_GENES}

STAR \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${MM_FASTA} ${Reo_FASTA} \
--sjdbGTFfile ${merged_GENES} \
--sjdbGTFfeatureExon exon

# cat /workdir/dwm269/genomes/mm39_all/GENCODE_M31/GRCm39.genome.fa /workdir/dwm269/genomes/custom_sequences/reo_T1L.fasta > /workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM31_REOT1L-as/concatenated_genome.fa
