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

# grab input fasta
if os.path.exists(sys.argv[1]):
    f = open(sys.argv[1], "r")
else:
    print(".fasta file doesn't exist! Try again...")

fasta = list(SeqIO.parse(f, "fasta"))

num_tRNAs = len(fasta)
print("     File contains " + str(num_tRNAs) + " records...")

# Build GFF file
## source: https://biopython.org/wiki/GFF_Parsing
## source: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec36
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

out_file = sys.argv[2]

print("     Building GFF file '" + out_file + "'...")

# initialize .gtf file
tmp_file = open(out_file, "w+")
tmp_file.write("##description: tRNA-only reference\n##provider: GtRNADB\n##format: gtf\n")
tmp_file.close()

# for i in list(range(0,5)):
for i in list(range(0,num_tRNAs)):

    # print(i)
    # print(re.findall(r'\((.*?)\)', fasta[i].description)[-1])

    id = fasta[i].id
    rec = SeqRecord(fasta[i].seq, fasta[i].id)

    qualifiers = {
        "source": "GtRNAdb",
        "score": re.findall(r"[-+]?\d*\.\d*\d+", fasta[i].description),
        "other": [fasta[i].description],
        "ID": fasta[i].id,
    }

    # real genomic start/end positions, chrom, and strand
    # start = re.findall(r":[-+]?\d*\-", fasta[i].description)[0][1:-1]
    # end = re.findall(r":[-+]?\d*\-+", fasta[i].description)[0][1:-1]
    # chrom = re.findall(r"chr?\S*:", fasta[i].description)[0][0:-1]
    # strand = re.findall(r'\((.*?)\)', fasta[i].description)[-1]

    # artificial values to match the tRNA-only reference
    start = 1
    end = len(fasta[i].seq)
    chrom = fasta[i].id
    strand = "+"

    top_feature = SeqFeature(
        FeatureLocation(int(start), int(end)),
        type="tRNA",
        strand=int(strand+'1'),
        qualifiers=qualifiers
    )
    rec.features = [top_feature]

    # buld string to write using rec
    to_write = chrom+"\t"+"GtRNADB"+"\t"+"exon"+"\t"+str(start)+"\t"+str(end)+"\t"+qualifiers['score'][0]+"\t"+str(strand)+"\t"+"0"+"\t"+ "gene_name "+ str(rec.id)+"; transcript_name "+ str(rec.id)+"; gene_biotype 'tRNA'\n"
    # print(to_write)
    tmp_file = open(out_file, "a")  # append mode
    tmp_file.write(to_write)
    tmp_file.close()

print("Done building "+out_file+"!")
