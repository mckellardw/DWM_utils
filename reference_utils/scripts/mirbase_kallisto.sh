#!/usr/bin/bash

# conda activate kallisto1
# cDNA_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.transcripts.fa.gz"

# GENOME_FASTA="/workdir/dwm269/genomes/mirbase/mature.fa"
# GENES_GTF="/workdir/dwm269/genomes/mirbase/mature.gtf"

GENOME_FASTA="/workdir/dwm269/genomes/mm39_all/STAR_smRNA/smRNA_Mus_musculus.fa"
GENES_GTF="/workdir/dwm269/genomes/mirbase/mature.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/kallisto_mirbase_k21"
cDNA_FASTA="cdna.fa"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

kb ref \
-i transcriptome.idx \
-g transcripts_to_genes.txt \
-k 21 \
-f1 ${cDNA_FASTA} \
${GENOME_FASTA} ${GENES_GTF}
