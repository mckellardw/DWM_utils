
# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STAR command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
rule STAR_align:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/postTrim_R1.fq.gz", 
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/postTrim_R2.fq.gz",
    output:  #TODO- add more output files?
        BAM="{OUTDIR}/{SAMPLE}/star/Aligned.sortedByCoord.out.bam",
        UNMAPPED1="{OUTDIR}/{SAMPLE}/star/Unmapped.out.mate1",
        UNMAPPED2="{OUTDIR}/{SAMPLE}/star/Unmapped.out.mate2",
        COUNTS="{OUTDIR}/{SAMPLE}/star/ReadsPerGene.out.tab",
    params:
        R1_FQ=lambda wildcards: R1_FQS[wildcards.SAMPLE],
        R2_FQ=lambda wildcards: R2_FQS[wildcards.SAMPLE],
        STAR_REF=config["STAR_REF"],
    resources:
        mem=config["MEMLIMIT"],
    threads: config["CORES"]
    shell:
        """
        mkdir -p $(dirname {output.BAM})

        STAR \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --genomeDir {params.STAR_REF} \
            --limitBAMsortRAM={resources.mem} \
            --readFilesIn {input.R1_FQ} {input.R2_FQ} \
            --outFileNamePrefix $(dirname {output.BAM})/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI NM MD AS nM jM jI XS \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterMatchNmin 16 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMultimapNmax 50
        """
            # --genomeLoad LoadAndKeep \

# Run fastqc on unmapped reads
rule rename_compress_unmapped:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/star/Unmapped.out.mate1",
        R2_FQ="{OUTDIR}/{SAMPLE}/star/Unmapped.out.mate2",
    output:
        R1_FQ="{OUTDIR}/{SAMPLE}/star/unmapped_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/star/unmapped_R2.fq.gz",
    threads: config["CORES"]
    shell:
        """
        mv {input.R1_FQ} {wildcards.OUTDIR}/{wildcards.SAMPLE}/star/unmapped_R1.fq
        mv {input.R2_FQ} {wildcards.OUTDIR}/{wildcards.SAMPLE}/star/unmapped_R2.fq
        pigz -p {threads} {wildcards.OUTDIR}/{wildcards.SAMPLE}/star/unmapped_R1.fq
        pigz -p {threads} {wildcards.OUTDIR}/{wildcards.SAMPLE}/star/unmapped_R2.fq
        """
