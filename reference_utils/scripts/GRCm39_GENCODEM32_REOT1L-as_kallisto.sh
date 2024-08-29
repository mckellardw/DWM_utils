#!/usr/bin/bash

# conda activate kallisto1
GENOME_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M32/GRCm39.genome.fa"
FASTA_cDNA="/workdir/dwm269/genomes/mm39_all/GENCODE_M32/gencode.vM32.transcripts.fa.gz"
GENES_GTF="/workdir/dwm269/genomes/mm39_all/GENCODE_M32/gencode.vM32.annotation.gtf"

Reo_FASTA="/workdir/dwm269/genomes/custom_sequences/reo_T1L.fasta"
Reo_GENES="/workdir/dwm269/genomes/custom_sequences/reo_T1L_antisense.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/kallisto_GRCm39_GENCODEM32_REOT1L-as"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

kb ref \
    -i mixed_index.idx \
    -g mixed_t2g.txt \
    -f1 mixed_cdna.fa \
    ${GENOME_FASTA},${Reo_FASTA} \
    ${GENES_GTF},${Reo_GENES}
