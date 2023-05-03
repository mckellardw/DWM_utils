#!/usr/bin/bash

# conda activate kallisto1
GENOME_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/GRCm39.genome.fa"
GENES_GTF="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf"

Reo_FASTA="/workdir/dwm269/genomes/reo_T1L.fasta"
Reo_GENES="/workdir/dwm269/genomes/reo_T1L_antisense.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/kallisto_GRCm39_GENCODEM28_REOT1L-as_k21"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

kb ref \
-i mixed_index.idx \
-g mixed_t2g.txt \
-f1 mixed_cdna.fa \
-k 21 \
${GENOME_FASTA},${Reo_FASTA} \
${GENES_GTF},${Reo_GENES}
