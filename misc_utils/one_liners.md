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
