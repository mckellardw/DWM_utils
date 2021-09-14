#####################################
## build_STAR_smRNA_ref.py
## Written by: David McKellar
## Description: This script takes in a .fasta from mirbase and generates a .gtf and STAR reference for a particular species
#####################################

# usage:
#      python build_smRNA_STAR_ref.py mirbase_mature_fasta gtrnadb_fasta species output_directory
#
# example:
#      /home/dwm269/miniconda3/envs/STARsolo/bin/python build_STAR_smRNA_ref.py /workdir/dwm269/genomes/mirbase/mature.fa /workdir/dwm269/genomes/mm39_all/gtrnadb/mm39-mature-tRNAs.fa "Mus musculus" /workdir/dwm269/genomes/mm39_all/mm39_smRNA_STAR 20

import sys
import os.path
import re
from Bio import SeqIO

# Dictionary of species shorthands
#TODO
species_dict = {
    "Mus musculus" : "mmu"
}

## .gtf structure
## source: https://useast.ensembl.org/info/website/upload/gff.html
fasta_mir_path = sys.argv[1]
fasta_tr_path = sys.argv[2]
species = str(sys.argv[3])#"mus musculus"
outdir = sys.argv[4]
num_threads = sys.argv[5]

fasta_out_path = outdir+"/"+"smRNA_"+ species.replace(" ","_")+".fa"
gtf_out_path = outdir+"/"+"smRNA_"+ species.replace(" ","_")+".gtf"

print("Building small RNA reference for "+ species + ":")
print("    fasta output location:  "+ fasta_out_path)
print("    gtf output location:    " + gtf_out_path)
print("    number threads:         " + str(num_threads))

if not os.path.isdir(outdir):
    os.mkdir(outdir)

print("Using the file '" + fasta_mir_path + "' as a miRNA database...")

# grab input fasta for miRNAs
if os.path.exists(fasta_mir_path) and re.search(".gz",fasta_mir_path):
    #TODO
    print("gunzip your fasta first please...")
    exit()
    # gzip.decompress(fasta_mir_path)
    # f = os.open(str(fasta_mir_path)[0:-3], "r")
elif os.path.exists(fasta_mir_path):
    f_mir = open(fasta_mir_path, "r")
else:
    print(".fasta file doesn't exist! Try again...")

fasta_mir = list(SeqIO.parse(f_mir, "fasta"))
num_miRNAs = len(fasta_mir)
print("     Found " + str(num_miRNAs) + " records...")

# grab input fasta for tRNAs
print("Using the file '" + fasta_tr_path + "' as a tRNA database...")

if os.path.exists(fasta_tr_path) and re.search(".gz",fasta_tr_path):
    #TODO
    print("gunzip your fasta first please...")
    exit()
    # gzip.decompress(fasta_mir_path)
    # f = os.open(str(fasta_mir_path)[0:-3], "r")
elif os.path.exists(fasta_tr_path):
    f_tr = open(fasta_tr_path, "r")
else:
    print(".fasta file doesn't exist! Try again...")

fasta_tr = list(SeqIO.parse(f_tr, "fasta"))
num_tRNAs = len(fasta_tr)
print("     Found " + str(num_tRNAs) + " records...")

#TODO: check fasta for desired species


# Add entries for each type of RNA to fasta...
print("Subsetting fasta file for "+species+"...")
print("     Checking for '"+species+"' or '"+species_dict[species]+"'")

# open up new fasta
fasta_out_f = open(fasta_out_path,"w")

print("     Adding miRNAs...")
for i in list(range(0,num_miRNAs)):
    if species_dict[species] in fasta_mir[i].id:
        fasta_out_f.write(">" + fasta_mir[i].id + "\n")
        fasta_out_f.write(str(fasta_mir[i].seq).replace("U", "T") + "\n") #switch U -> T

print("     Adding tRNAs...")
for i in list(range(0,num_tRNAs)):
    if species.replace(" ","_") in fasta_tr[i].id:
        fasta_out_f.write(">" + fasta_tr[i].id + "\n")
        fasta_out_f.write(str(fasta_tr[i].seq) + "\n")
        fasta_out_f.write(str(fasta_mir[i].seq).replace("U", "T") + "\n") #switch U -> T

fasta_out_f.close()

print("     Done!")
print("")

# Build GFF file
## source: https://biopython.org/wiki/GFF_Parsing
## source: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec36
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

print("Building GTF file (" + gtf_out_path + ")...")

fasta_out_f = open(fasta_out_path, "r")
fasta_out = list(SeqIO.parse(fasta_out_f, "fasta"))

# initialize .gtf file
tmp_file = open(gtf_out_path, "w+")
tmp_file.write("##description: small RNA reference\n##provider: miRBase and GtRNAdb\n##format: gtf\n")
# tmp_file.close()

for i in list(range(0,len(fasta_out))):

    id = fasta_out[i].id
    rec = SeqRecord(fasta_out[i].seq, fasta_out[i].id)

    qualifiers = {
        "source": "miRBase and GtRNADB",
        "score": ".", #re.findall(r"[-+]?\d*\.\d*\d+", fasta_out[i].description),
        "other": [fasta_out[i].description],
        "ID": fasta_out[i].id,
    }

    # artificial values to match the reference
    start = 1
    end = len(fasta_out[i].seq)
    chrom = fasta_out[i].id
    strand = "+"

    # Check "A" content for spurious capture
    pct_A = fasta_out[i].seq.count("A") / len(fasta_out[i].seq)

    # top_feature = SeqFeature(
    #     FeatureLocation(int(start), int(end)),
    #     type="miRNA",#TODO
    #     strand=int(strand+'1'),
    #     qualifiers=qualifiers
    # )
    # rec.features = [top_feature]

    # buld string to write using rec
    to_write = \
    chrom + "\t" + \
    "miRBase/GtRNADB" + "\t" + \
    "transcript" + "\t" + \
    str(start) + "\t"+\
    str(end)+ "\t" +\
    qualifiers['score'][0] +"\t" +\
    str(strand) +"\t" +\
    "."+"\t"+ \
    "gene_name \""+ str(rec.id)+"\"; transcript_name \"" + str(rec.id) +"\"; gene_id \""+ fasta_out[i].id +"\";\n"
    #+ "; gene_biotype 'miRNA'\n" # everything else #TODO- add biotype

    tmp_file.write(to_write)

tmp_file.close()

print("    Done!")

print("Building STAR reference...")

cmd = "STAR --runMode genomeGenerate " + \
" --genomeFastaFiles " + fasta_out_path + \
" --genomeDir " + outdir + \
" --runThreadN " + num_threads + \
" --genomeSAindexNbases 7" + \
" --sjdbGTFfile " + gtf_out_path + \
" --sjdbGTFfeatureExon transcript"

os.system(cmd)

print("    Done!")
