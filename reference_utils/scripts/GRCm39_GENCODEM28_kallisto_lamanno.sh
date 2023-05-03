#!/usr/bin/bash

# conda activate kallisto1
# cDNA_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.transcripts.fa.gz"
GENOME_FASTA="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/GRCm39.genome.fa"
GENES_GTF="/workdir/dwm269/genomes/mm39_all/GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/kallisto_GRCm39_GENCODEM28_lamanno"
cDNA_FASTA="cdna.fa"
INTRON_FASTA="introns.fa"
cDNA_T2C="cDNA.t2c"
INTRON_T2C="introns.t2c"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

kb ref \
-i transcriptome.idx \
-g transcripts_to_genes.txt \
--workflow lamanno \
-f1 ${cDNA_FASTA} \
-f2 ${INTRON_FASTA} \
-c1 ${cDNA_T2C} \
-c2 ${INTRON_T2C} \
${GENOME_FASTA} ${GENES_GTF}
