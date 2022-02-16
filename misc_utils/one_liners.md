# One-liners: useful code snippets


## Linux utils

- Check memory usage in all sub-directories (of current directory), then sort outputs
```
du -sh /path/to/base/directory/* | sort -hr
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

## TXG helpers
- Number of reads per cell/spot barcode (.bam tag `CB`)
 ```
samtools view sub.bam | grep CB:Z: | sed 's/.*CR:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_umi.txt
 ```

 - Number of reads per UMI (.bam tag `CR`). *Note*, UMIs are repeated, so best to use the `sS` tag; even better, concatenate the `CB` and `CR`; best is to actually run UMI/sequence-aware deduplication and count occurrences that way (try `umi_tools`)
 ```
 samtools view sub.bam | grep CR:Z: | sed 's/.*CR:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_umi.txt
 ```

## samtools

```

```
##
- Convert .gtf (from GENCODE) to a .bed
```
tail -n +6 /path/to/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf | \
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' | \
gtf2bed --max-mem=16G > gencode.vM28.chr_patch_hapl_scaff.annotation.bed
```
