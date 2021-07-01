#####################################
## build_GtRNAdb.py
## Written by: David McKellar
## Description: This script takes in a .fasta from GtRNAb and generates a .GFF, which can be used to build a reference for alignment of RNAseq data
#####################################

# usage:
#      python build_GtRNAdb.py GtRNAdb.fasta output_name.gff

import sys
import os.path
# import pandas as pd
from Bio import SeqIO

## .gtf structure
# 1) seqname - name of the chromosome or scaffold
# 2) source - name of the program that generated this feature ("GtRNAdb")
# 3) feature - ("tRNA")
# 4) start - Start position* of the feature, with sequence numbering starting at 1.
# 5) end - End position* of the feature, with sequence numbering starting at 1.
# 6) score - A floating point value.
# 7) strand - defined as + (forward) or - (reverse).
# 8) frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

print("Using the file '" + sys.argv[1] + "' as a tRNA database...")

if os.path.exists(sys.argv[1]):
    f = open(sys.argv[1], "r")

fasta = list(SeqIO.parse(f, "fasta"))

num_tRNAs = len(fasta)
print("     File contains " + str(num_tRNAs) + " records...")

# Build GFF file
## source: https://biopython.org/wiki/GFF_Parsing

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

print("     Building GFF file...")

for i in list(range(0,2)):

    print(re.findall(r":[-+]?\d*\-", fasta[i].description))

    id = fasta[i].id
    rec = SeqRecord(fasta[i].seq, fasta[i].id)

    qualifiers = {
        "source": "GtRNAdb",
        "score": re.findall(r"[-+]?\d*\.\d*\d+", fasta[i].description),
        "other": [fasta[i].description],
        "ID": fasta[i].id,
    }

    chrom = re.findall(r"chr?\S*:", fasta[i].description)[0][0:-1]
    strand = re.findall(r'\((.*?)\)', fasta[i].description)[-1]+'1'

    start = re.findall(r":[-+]?\d*\-", fasta[i].description)[0][1:-1]
    stop = re.findall(r":[-+]?\d*\-+", fasta[i].description)[0][1:-1]

    print(strand)

    # top_feature = SeqFeature(
    #     FeatureLocation(int(start), int(stop)), type="tRNA", strand=strand, qualifiers=qualifiers
    # )
    # rec.features = [top_feature]

out_file = sys.argv[2]

with open(out_file, "w") as out_handle:
    GFF.write([rec], out_handle)
