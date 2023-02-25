# Parallelized fastq splitting into N files, maintaining read order

# Usage:
#   python splitNfqs.py in.fq.gz 12 12

#TODO: add output options

# imports
import sys
import os
import gzip
import subprocess

# params
fq = sys.argv[1]  # Path to fq file
n_chunks = sys.argv[2] # Number of files to chunk the .fastq into
n_cores = sys.argv[3]  # Number of cores to parallelize with
n_total_lines = sys.argv[4]  # Optional: number of total lines in fastq file (n_reads*4)

# Param checks
## Make sure core count is lower/equal to chunk count
if n_cores > n_chunks:
    n_cores = n_chunks

## Count number of lines if not given
if n_total_lines is None:    
    n_total_lines = int(subprocess.check_output(['zcat', fq, '|', 'wc', '-l']).split()[0])

# Function to extract out a section of the fastq file
def extractReads(fq_in, fq_out, start, stop):
    os.system(
        f"""
        zcat {fq_in} | sed -n '{start},{stop}' > {fq_out}
        gzip {fq_out}
        """
    )

# Compute start/stop pairs for reads 
lines_per_file = n_total_lines // n_chunks
stops = [lines_per_file*n  for n in range(1,n_chunks)]
stops = [stop-stop%4 for stop in stops]            # Ensure that stop positions aren't in the middle of a read (not divisible by 4)
stops.append(n_total_lines)                        # Add final stop (EOF)

starts = [stops[i] - 1 for i in range(1,n_chunks)] # Get start lines
starts.insert(0,1)                                 # Add beginning of file to starts

if "fq" in fq:
    fq_out_list = [f"{fq.replace('fq.gz','')}_{str(n).zfill(3)}.fq" for n in range(0,n_chunks)]
elif "fastq" in fq:
    fq_out_list = [f"{fq.replace('fastq.gz','')}_{str(n).zfill(3)}.fastq" for n in range(0,n_chunks)]

# Extract reads and write new files
if n_cores > 1: # Parallelize with `multiprocessing`
    import multiprocessing
    items = [(fq, fq_out_list[i], starts[i], stops[i]) for i in range(0,n_chunks)]     
    with multiprocessing.Pool(n_cores) as pool:
        multi_out = pool.starmap(extractReads, items)
else: # Single thread
    for i in list(range(0,n_chunks)):
        extractReads(fq, fq_out_list[i], starts[i], stops[i])

# Compress output .fastqs
os.system(
        f"""
        pigz -p{n_cores} {" ".join(fq_out_list)}
        """
    )