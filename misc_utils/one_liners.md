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

## samtools

```

```
