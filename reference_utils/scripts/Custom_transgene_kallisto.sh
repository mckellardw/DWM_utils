#!/usr/bin/bash

# conda activate kallisto1
GENOME_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/GRCm39.genome.fa"
GENES_GTF="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf"

TRANS_FASTA="/workdir/dwm269/genomes/custom_sequences/FLExDUX4/FLExDUX4.fa"
TRANS_GTF="/workdir/dwm269/genomes/custom_sequences/FLExDUX4/FLExDUX4.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/kallisto_GRCm39_GENCODEM28_FLExDUX4"
MERGED_FASTA="GRCm39_FLExDUX4.fa"
MERGED_GTF="GRCm39_FLExDUX4.gtf"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# FlexDUX4 transgene source: https://med.unr.edu/jones-lab/resources/flexdux4-mice

# Merge .gtf and .fasta for any future analyses
cat ${GENOME_FASTA} ${TRANS_FASTA} > ${OUTDIR}/${MERGED_FASTA}
cat ${GENOME_GTF} ${TRANS_GTF} > ${OUTDIR}/${MERGED_GTF}

# Build reference with kallisto
kb ref \
-i mixed_index.idx \
-g mixed_t2g.txt \
-f1 mixed_cdna.fa \
-k 31 \
${GENOME_FASTA},${TRANS_FASTA} \
${GENES_GTF},${TRANS_GTF}
