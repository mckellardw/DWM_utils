#!/usr/bin/bash

# Usage:
#   kb.sh $OUTDIR $KB_IDX $WHITELIST $CHEMISTRY $LOG $THREADS $MEMLIMIT $R1FQ $R2FQ

# Get params
OUTDIR=$1
KB_IDX=$2
WHITELIST=$3
CHEMISTRY=$4
LOG=$5
THREADS=$6
MEMLIMIT=$7
R1FQ=$8
R2FQ=$9

# Check params
#TODO

# Set up output directory
mkdir -p ${OUTDIR}
cd ${OUTDIR}
touch ${LOG}

# Pseudoalign and generate .bus file
echo "~~~Pseudoaligning with `kallisto bus`...\n\n " >> ${LOG}
kallisto bus \
-i ${KB_IDX} \
-x ${CHEMISTRY} \
--fr-stranded \
-o ${OUTDIR} \
-t ${THREADS} \
--verbose \
${R1FQ} ${R2FQ} 2>> ${LOG}

# Sort .bus file
echo "~~~Sorting output bus...\n\n " >> ${LOG}
bustools sort \
-t ${THREADS} \
-m ${MEMLIMIT} \
-o output.sorted.bus \
output.bus 2>> ${LOG}

# Correct cell/spot barcodes
echo "~~~Correcting barcodes...\n\n " >> ${LOG}
bustools correct \
--whitelist ${WHITELIST} \
-o output.corrected.bus \
output.sorted.bus 2>> ${LOG}

# Convert bus file to text for easier counting
echo "~~~Converting bus to text...\n\n " >> ${LOG}
bustools text \
-o output.corrected.txt \
output.corrected.bus
