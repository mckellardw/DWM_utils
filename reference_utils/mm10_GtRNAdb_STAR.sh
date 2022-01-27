#!/usr/bin/bash

NCORES=30
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

# Build custom annotations file from the .fasta
TRNA_ANNOTATIONS=${DATADIR}"/mm10-tRNAs.gff"
/home/dwm269/miniconda3/envs/STARsolo/bin/python /home/dwm269/DWM_utils/ref_utils/scripts/build_GtRNAdb_gtf.py ${FASTA} ${TRNA_ANNOTATIONS}

ReadLengthMinus1=54

cd ${OUTDIR}

${STAR_EXEC} \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--sjdbGTFtagExonParentTranscript Parent \
--genomeSAindexNbases 6 \
--genomeFastaFiles ${FASTA} \
--sjdbGTFfile ${TRNA_ANNOTATIONS} \
--sjdbOverhang ${ReadLengthMinus1}



#
