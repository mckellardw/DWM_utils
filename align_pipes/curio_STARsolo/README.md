# `curio_STARsolo`
Preprocessing, alignment, QC, and quantification workflow for Curio Bioscience data (SlideSeq/Seeker)
**David W. McKellar**


## Dependencies:
- `cutadapt` [v3.4](https://cutadapt.readthedocs.io/en/stable/)
- `fastqc` [vv0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `STAR` [v2.7.10a](https://github.com/alexdobin/STAR)
- `qualimap` [v2.2.2a](http://qualimap.conesalab.org/)
- `vsearch` [v2.17.0_linux_x86_64](https://github.com/torognes/vsearch)
- `BLAST`

## Outputs:
```
SAMPLE_ID/
├── cutadapt_polyA_report.txt
├── cutadapt_polyG_report.txt
├── postTrim_fastqc_R2_out
│   ├── SAMPLE_ID_R2_final_fastqc.html
│   └── SAMPLE_ID_R2_final_fastqc.zip
├── preTrim_fastqc_R2_out
│   ├── SAMPLE_ID_R2_fastqc.html
│   └── SAMPLE_ID_R2_fastqc.zip
├── STARsolo
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.bam.bai
│   ├── Log.final.out
│   ├── Log.out
│   ├── Log.progress.out
│   ├── SJ.out.tab
│   ├── Solo.out
│   │   ├── Barcodes.stats
│   │   ├── Gene
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   ├── GeneFull
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   └── Velocyto
│   │       ├── Features.stats
│   │       ├── filtered
│   │       │   ├── ambiguous.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   ├── features.tsv.gz
│   │       │   ├── spliced.mtx.gz
│   │       │   └── unspliced.mtx.gz
│   │       ├── raw
│   │       │   ├── ambiguous.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   ├── features.tsv.gz
│   │       │   ├── spliced.mtx.gz
│   │       │   └── unspliced.mtx.gz
│   │       └── Summary.csv
│   ├── Unmapped.out.mate1.fastq.gz
│   └── Unmapped.out.mate2.fastq.gz
├── tmp
│   ├── whitelist_1.txt
│   ├── whitelist_2.txt
│   └── whitelist.txt
├── Unmapped_fastqc_out
│   ├── Unmapped.out.mate1_fastqc.html
│   ├── Unmapped.out.mate1_fastqc.zip
│   ├── Unmapped.out.mate2_fastqc.html
│   └── Unmapped.out.mate2_fastqc.zip
└── Unmapped.out.mate2_blastResults.txt
```


## Helpful links
- [Barcode download](http://3.222.73.200/)
