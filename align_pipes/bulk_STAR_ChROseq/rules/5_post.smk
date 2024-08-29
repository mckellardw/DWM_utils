# generate stranded bigWigs for plotting/visualization
rule bamToSplitBigWig:
    input:
        BAM="{OUTDIR}/{SAMPLE}/star/Aligned.sortedByCoord.out.bam",
        BAI="{OUTDIR}/{SAMPLE}/star/Aligned.sortedByCoord.out.bam.bai",
    output:
        POS_BW="{OUTDIR}/{SAMPLE}/star/Aligned.sortedByCoord.out_plus.bw",
    params:
        STAR_REF=config["STAR_REF"],
    threads: config["CORES"]
    shell:
        """
        scripts/bam2splitBigWig.sh \
            {input.BAM} \
            {threads} \
            $(dirname {input.BAM}) \
            {STAR_REF}/chrNameLength.txt
        """

#TODO: rule to generate some qc plots on TSS's, etc