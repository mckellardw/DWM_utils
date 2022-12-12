#!/usr/bin/bash
# split_dedup.sh - a bash script that sesuplicates TXG bam files, aligned with STARsolo
#                - deduplicates chr-by-chr, to reduce run time and memory requirements
# Usage:
# bash chr_split_dedup.sh /path/to/file.bam /path/to/whitelist num_cores /path/to/output/directory file.LOG

# bash /home/dwm269/DWM_utils/seq_utils/chr_split_dedup.sh ds.bam /workdir/dwm269/totalRNA/STRS-HD/data/STARsolo/SH2A/tmp/whitelist.txt 24 ./OUTPUT_DEDUP DEDUP.LOG

# chr_split_dedup (v1.0) - input a .bam file, output the deduplicated .bam file in OUTDIR

#TODO: python script instead of bash to make parallelization a bit easier?

INBAM=$1 # path to .bam file (sorted & indexed already!)
BB=$2 # path to barcode whitelist
CORE=$3 # number of cores for parallelization
OUTDIR=$4 # output/deduped bam path
LOG=$5


mkdir -p ${OUTDIR}/tmp
cd ${OUTDIR}/tmp

LOGFULL=${OUTDIR}/${LOG}

echo "Log for chr_split_dedup:" > ${LOGFULL}

echo "Log info can be found in "${LOGFULL}

PREFIX=`echo ${INBAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev`

echo "Using "${PREFIX}" as file prefix..." >> ${LOGFULL}

OUTBAM=${OUTDIR}/${PREFIX}_dedup.bam

echo ".bam file location:    "${INBAM} >> ${LOGFULL}
echo "Max cores:             "${CORE} >> ${LOGFULL}
echo "Output location:       "${OUTDIR} >> ${LOGFULL}
echo >> ${LOGFULL}
echo >> ${LOGFULL}

# Check params...
if [ ! -f ${INBAM} ]
then
    echo "Can't find ${INBAM}"
fi

# Look for .bam index file
if [ ! -f ${INBAM}.bai ]
then
    echo "Indexing input .bam file because you forgot to..."
    samtools index -@ ${CORE} ${INBAM}
fi

# Remove reads that don't have a barcode (CB)
echo "Removing reads without 'CB' tags..." >> ${LOGFULL}
date >> ${LOGFULL}

# Option w/ grep, which doesn't require the barcode whitelist
# samtools view -h Aligned.sortedByCoord.out.bam | grep -v "CB:Z:-" > yestag.sam

samtools view -1 -b \
-@ ${CORE} \
--tag-file CB:${BB} \
${INBAM} \
> ${OUTDIR}/tmp/filter.bam
echo "Done." >> ${LOGFULL}
echo >> ${LOGFULL}

# split bam by chromosome/strand
echo "Splitting by chromosome..." >> ${LOGFULL}
date >> ${LOGFULL}
bamtools split \
-in ${OUTDIR}/tmp/filter.bam \
-reference
echo "Done." >> ${LOGFULL}
echo >> ${LOGFULL}

BAMLIST=${OUTDIR}/tmp/BAMLIST.tmp
ls *REF_*.bam > ${BAMLIST}

# index split .bam's
echo "Indexing split .bam files..." >> ${LOGFULL}
date >> ${LOGFULL}
while read BAM; do
  samtools index \
  -@ ${CORE} \
  ${BAM}
done < ${BAMLIST}
echo "Done." >> ${LOGFULL}
echo >> ${LOGFULL}

# dedup resulting bams one-by-one
echo "Deduplicating split .bam files..." >> ${LOGFULL}
date >> ${LOGFULL}
#TODO: parallelize
while read BAM; do
  samtools index -@ ${CORE} ${BAM}
  umi_tools dedup \
  -I ${BAM} \
  --extract-umi-method=tag \
  --umi-tag=UB \
  --cell-tag=CB \
  --method=unique \
  --per-cell \
  --unmapped-reads=discard \
  --log ${OUTDIR}/umitools.log \
  -S ${OUTDIR}/tmp/dedup_${BAM}
done < ${BAMLIST}
# --output-stats={params.OUTPUT_PREFIX} \
echo "Done." >> ${LOGFULL}
echo >> ${LOGFULL}

# merge, sort, and index dedup'ed .bams
echo "Merging, sorting, and indexing deduplicated .bam files..." >> ${LOGFULL}
date >> ${LOGFULL}
samtools merge \
-f \
-o ${OUTDIR}/tmp/dedup_merge.bam \
-@ ${CORE} \
--write-index \
${OUTDIR}/tmp/dedup_*REF_*.bam

samtools sort ${OUTDIR}/tmp/dedup_merge.bam > ${OUTBAM}
samtools index -@ ${CORE} ${OUTBAM} && rm -r ${OUTDIR}/tmp
echo "Done." >> ${LOGFULL}
echo >> ${LOGFULL}
