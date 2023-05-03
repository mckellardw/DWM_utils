#!/usr/bin/bash

NCORES=30
# FASTA_DIR="/workdir/dwm269/genomes/mirbase/mature.fa"
# GENES_DIR="/workdir/dwm269/genomes/mirbase/mature.gtf"

FASTA_DIR="/home/dwm269/usr/bin/miRge3_Lib/mouse/fasta.Libs/mouse_mature_miRBase.fa"
GENES_DIR="/workdir/dwm269/genomes/mirbase/mature.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_mirbase_mirge"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
--runThreadN ${NCORES} \
--limitGenomeGenerateRAM=64000000000 \
--runMode genomeGenerate \
--genomeSAindexNbases 6 \
--sjdbOverhang 1 \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${FASTA_DIR} \
--sjdbGTFfile ${GENES_DIR} \
--sjdbGTFfeatureExon exon
