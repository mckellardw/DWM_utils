
# Merge .fastq files (in case more than one sesquencing run was performed)
rule merge_fastqs:
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/merged_R1.fq.gz"),
        R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/merged_R2.fq.gz"),
    params:
        TMP_DIR="{OUTDIR}/{SAMPLE}/tmp",
        R1_FQ=lambda wildcards: R1_FQS[wildcards.SAMPLE],
        R2_FQ=lambda wildcards: R2_FQS[wildcards.SAMPLE],
    threads: config["CORES"]
    run:
        if len(params.R1_FQ.split(" ")) == 1:  # shell for single fastq input
            shell("cp {params.R1_FQ} {output.R1_FQ}")
            shell("cp {params.R2_FQ} {output.R2_FQ}")
        else:  # shell enablinging multi-fast input; concatenate inputs
            print(
                "Concatenating",
                len(params.R1_FQ.split(" ")),
                ".fastq's for",
                wildcards.SAMPLE,
            )
            shell("mkdir -p {params.TMP_DIR}")
            shell("zcat {params.R1_FQ} > {params.TMP_DIR}/merged_R1.fq")
            shell("zcat {params.R2_FQ} > {params.TMP_DIR}/merged_R2.fq")
            shell("pigz -p {threads} {params.TMP_DIR}/*.fq")



# Index .bam file
rule index_BAM:
    input:
        BAM="{BAM}",
    output:
        BAI="{BAM}.bai",
    # wildcard_constraints:
    #     BAM=".*\.(bam)$"
    # resources:
    threads: config["CORES"]
    shell:
        """
        samtools index -@ {threads} {input.BAM}
        """