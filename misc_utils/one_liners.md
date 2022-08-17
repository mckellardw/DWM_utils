# One-liners: useful code snippets


## Linux utils

- Check disk usage in all sub-directories (of current directory), then sort outputs
```
du -sh /path/to/base/directory/* | sort -hr
```

- Check memory usage acrosss current processes, by ****user***
```
ps aux | awk '{arr[$1]+=$4}; END {for (i in arr) {print i,arr[i]}}' | sort -k2
```

## `awk`, `sed`, etc.

- Convert .csv to .tsv (note that .gtf files are also .tsv's)
```
sed -E 's/("([^"]*)")?,/\2\t/g' file.csv > file.tsv
```

- Add empty string to all entries in a .gtf which are missing 'transcript_id'
```
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }'
```

- Filter .fastq [by read length](https://www.biostars.org/p/66996/)
*Note* that if you use this line in a snakemake shell call, both the `\` and the `{}` have to be escaped (`\`->`\\` and `{}`-> `{{}}`)
*Also note* this only works for a single-end read- use cutadapt for paired end (or reads w/ UMIs/barcodes)
```
zcat your.fastq.gz | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= ${MIN_LENGTH} && length(seq) <= ${MAX_LENGTH}) {print header, seq, qheader, qseq}}' > filtered.fastq
```

- Extract (unique) sequences from a .fastq.gz file
    - Remove `|sort | uniq` if you want all sequences...
```
zcat reads.fastq.gz | awk '(NR%4==2)' | sort | uniq
```

## TXG helpers
- Number of reads per cell/spot barcode (.bam tag `CB`)
 ```
samtools view sub.bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_umi.txt
 ```

 - Number of reads per UMI (.bam tag `CR`). *Note*, UMIs are repeated, so best to use the `sS` tag; even better, concatenate the `CB` and `CR`; best is to actually run UMI/sequence-aware deduplication and count occurrences that way (try `umi_tools`)
 ```
 samtools view sub.bam | grep CR:Z: | sed 's/.*CR:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_umi.txt
 ```

## samtools
- Subset bam by strand
```
samtools view -F 0x10 -b file.bam > file_pos.bam
```

## .fastq finagling

- Get top 1000 most abundant sequences from a file.fastq, output as an out.fa
```
vsearch --sortbysize file.fastq --topn 1000 --output out.fa
```

## .gtf finagling
- Convert .gtf (from GENCODE) to a .bed
```
tail -n +6 /path/to/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf | \
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' | \
gtf2bed --max-mem=16G > gencode.vM28.chr_patch_hapl_scaff.annotation.bed
```

- Filter a .gtf by gene biotype
```
cat gencode.vM28.chr_patch_hapl_scaff.annotation.gtf | grep -E 'transcript_type "miRNA"' > gencode.vM28.miRNAs.gtf
```
