#####################################
## build_GtRNAdb.py
## Written by: David McKellar
## Description: This script takes in a .fasta from mirbase and generates a .gtf, which can be used to build a reference for alignment of RNAseq data
#####################################

# usage:
#      python build_mirbase.py mature.fa.gz output_name.gtf species

import sys
import os.path
from Bio import SeqIO

## .gtf structure
## source: https://useast.ensembl.org/info/website/upload/gff.html

print("Using the file '" + sys.argv[1] + "' as a miRNA database...")

# grab input fasta
if os.path.exists(sys.argv[1]) and re.search(".gz",sys.argv[1]):
    #TODO
    print("gunzip your fasta first please...")
    exit()
    # gzip.decompress(sys.argv[1])
    # f = os.open(str(sys.argv[1])[0:-3], "r")

elif os.path.exists(sys.argv[1])):
    f = os.open(sys.argv[1], "r")

else:
    print(".fasta file doesn't exist! Try again...")

fasta = list(SeqIO.parse(f, "fasta"))

num_miRNAs = len(fasta)
print("     File contains " + str(num_miRNAs) + " records...")

#TODO- make species a passable argument
species = "mus musculus"

#TODO: check for desired species


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
tmp_file.write("##description: miRNA-only reference\n##provider: miRBase\n##format: gtf\n")
tmp_file.close()

for i in list(range(0,5)):
# for i in list(range(0,num_miRNAs)):

    id = fasta[i].id
    if re.search(species,id): #only include mouse!
        rec = SeqRecord(fasta[i].seq, fasta[i].id)

        qualifiers = {
            "source": "miRBase",
            "score": re.findall(r"[-+]?\d*\.\d*\d+", fasta[i].description),
            "other": [fasta[i].description],
            "ID": fasta[i].id,
        }

        # artificial values to match the miRNA-only reference
        start = 1
        end = len(fasta[i].seq)
        chrom = fasta[i].id
        strand = "+"

        top_feature = SeqFeature(
            FeatureLocation(int(start), int(end)),
            type="miRNA",
            strand=int(strand+'1'),
            qualifiers=qualifiers
        )
        rec.features = [top_feature]

        # buld string to write using rec
        to_write = chrom+"\t"+"miRBase"+"\t"+"exon"+"\t"+str(start)+"\t"+str(end)+"\t"+qualifiers['score'][0]+"\t"+str(strand)+"\t"+"0"+"\t"+ "gene_name "+ str(rec.id)+"; transcript_name "+ str(rec.id)+"; gene_biotype 'miRNA'\n"
        # print(to_write)
        tmp_file = open(out_file, "a")  # append mode
        tmp_file.write(to_write)
        tmp_file.close()

print("Done building "+out_file+"!")
