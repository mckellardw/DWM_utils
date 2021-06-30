#!/usr/bin/bash

NCORES=15
STAR_EXEC=/programs/STAR/STAR #biohpc current version of STAR

DATADIR=/workdir/dwm269/genomes/mm10_all/GtRNAdb
OUTDIR=/workdir/dwm269/genomes/mm10_all/mm10_GtRNAdb_STAR
mkdir -p ${OUTDIR}

cd ${DATADIR}

# Download fasta files with tRNA sequences

## High confidence tRNAs
# wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.fa

## Mature tRNAs
# wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-mature-tRNAs.fa

FASTA=${DATADIR}"/mm10-tRNAs.fa"

cd ${OUTDIR}

${STAR_EXEC} \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeSAindexNbases 6 \
--genomeFastaFiles ${FASTA}
