#!/usr/bin/bash

# conda activate kallisto1
# cDNA_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.transcripts.fa.gz"
GENOME_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/GRCm39.genome.fa"
GENES_GTF="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.miRNAs.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/kallisto_GRCm39_GENCODEM28_miRNAs"
cDNA_FASTA="cdna.fa"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

kb ref \
-i transcriptome.idx \
-g transcripts_to_genes.txt \
-f1 ${cDNA_FASTA} \
${GENOME_FASTA} ${GENES_GTF}
