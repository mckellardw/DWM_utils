# One-liners: useful code snippets
## Table of Contents
<!---toc start-->
  * [Other one-liner resources](#Other one-liner resources)
  * [General utils](#General utils)
  * [`awk`, `sed`, etc.](#`awk`, `sed`, etc.)
  * [bam wrangling](#bam wrangling)
  * [fastq finagling](#fastq finagling)
  * [gtf finagling](#gtf finagling)
<!---toc end-->

## Other one-liner resources
Other people's one-liners that are more useful...
- Stephen Turner [[link](https://github.com/stephenturner/oneliners#etc)]


## General utils

- Check disk usage in all sub-directories (of current directory), then sort outputs
```
du -sh /path/to/base/directory/* | sort -hr
```

- Check memory usage across current processes, by ****user***
```
ps aux | awk '{arr[$1]+=$4}; END {for (i in arr) {print i,arr[i]}}' | sort -k2
```

- Broad memory usage stats
```
cat /proc/meminfo
```
[Useful guide from redhat for interpreting output](https://access.redhat.com/solutions/406773)

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
    - Remove `| sort | uniq` if you want all sequences...
```
zcat reads.fastq.gz | awk '(NR%4==2)' | sort | uniq
```

- Remove a fixed length segment from all reads in a .fastq file
```
zcat file.fastq.gz | \
awk -v s=6 -v S=26 '(FNR-1) % 2 == 0 { name=$1; chr=$2; len=$3; next }
                    (FNR-2) % 4 == 0 { seq=substr($0,s,S) }
                                     { print name "." seq, chr, len
                                       print substr($0,n+1) }'
```

- Extract sequences in a .fasta based on a pattern in the header
*Note:* `FS` here is the "field separator", which may need to be changed for your .fasta
```
zcat gencode.vM31.transcripts.fa.gz | awk -v RS="\n>" -v FS="|" '$8=="rRNA" {print ">"$0}'
```
*Example output:*
```
>ENSMUST00000240339.1|ENSMUSG00000119695.1|-|-|n-R5s209-201|n-R5s209|119|rRNA|
GTCTAGGGCCATACCACCCTGAACGCGCCCAATCTCGTCTGTTCTCAGAAGCTAAGCAGG
GTTGGGCCTGGTTAGTACTTGGATGGGAGACTGTCCAGGATTACCGGGTGCTGTAGGAT
>ENSMUST00020182636.1|ENSMUSG00002076113.1|-|-|Gm55778-201|Gm55778|138|rRNA|
ATCTATGGCCATAGTGTATACGCACCTGATCTCATCTGATCTCAGATGAAATATCTGTGT
ATAACACAGTACAATTTACAATGCACACATACAGTGCTGGTCAGGAGAAAGGTTAAAGTA
GGATGTAGCCTGCAGGTC
>ENSMUST00000177688.3|ENSMUSG00000119706.1|-|-|n-R5s210-201|n-R5s210|119|rRNA|
GTCTATGGCCATACCACCCTGAACATGCCTGATCTCATCTGATCTCGGAAGCTAAGCAGG
GTCGGGCCTTGTTAGTACTTGGATGGGAGACCGCCTGAAAATACCAGGTTCTGTAGGCT
```

## .bam-boozling
- Number of reads per cell/spot barcode (.bam tag `CB`)
 ```
samtools view sub.bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_umi.txt
 ```

 - Number of reads per UMI (.bam tag `CR`). *Note*, UMIs are repeated, so best to use the `sS` tag; even better, concatenate the `CB` and `CR`; best is to actually run UMI/sequence-aware deduplication and count occurrences that way (try `umi_tools`)
 ```
 samtools view sub.bam | grep CR:Z: | sed 's/.*CR:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_umi.txt
 ```

- Subset bam by strand
```
samtools view -F 0x10 -b file.bam > file_pos.bam
```

- Remove reads that don't have a cell barcode tag
```
samtools view -h Aligned.sortedByCoord.out.bam | grep -v "CB:Z:-" > yestag.sam
```

## .fastq finagling

- Get top 1000 most abundant sequences from a file.fastq, output as an out.fa
```
vsearch --sortbysize file.fastq --topn 1000 --output out.fa
```

## .gtf jamboree
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
