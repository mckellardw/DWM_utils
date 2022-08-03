## seq_utils
Multipurpose scripts/one-liners for handling sequencing data

#### bam2splitBigWig
Convert .bam files to stranded .bigWigs! Outputs linear and log10 scaled .bw files for each strand, plus a merged (both strands with stranded sign, +/-, for coverage).

Usage:
```
bash /path/to/bam2splitBigWig.sh /path/to/file.bam num_cores /path/to/outdir /path/to/genome.chrom.sizes
```
