

## Install:
Build mamba/conda environment:
```
mamba create --name bulk_STAR_ChROseq  bioconda::snakemake bioconda::STAR bioconda::cutadapt bioconda::fastqc bioconda::qualimap bioconda::samtools conda-forge::pigz
mamba activate bulk_STAR_ChROseq
```

## Run:
```
snakemake -j 56
```