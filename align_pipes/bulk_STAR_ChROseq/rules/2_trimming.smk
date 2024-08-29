

# rule trimgalore:
#     input:
#         MERGED_R1_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R1.fq.gz",
#         MERGED_R2_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R2.fq.gz",
#     output:
#         # OUTDIR = directory('{OUTDIR}/{SAMPLE}/tmp/'),
#         A_TRIMMED_R1_FQ="{OUTDIR}/{SAMPLE}/tmp/trimmed_R1.fq.gz",
#         A_TRIMMED_R2_FQ="{OUTDIR}/{SAMPLE}/tmp/trimmed_R2.fq.gz",
#     params:
#         OUTDIR="{OUTDIR}/{SAMPLE}/tmp/",
#     # config['CORES']
#     threads: min([config["CORES"], 8])  
#     shell:
#         """
#         trim_galore \
#             --length 10 \
#             --2colour \
#             --paired \
#             -o {params.OUTDIR} \
#             --cores {threads} \
#             {input.MERGED_R1_FQ} {input.MERGED_R2_FQ}
#         """

rule cutadapt:
    input:
        R1_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R1.fq.gz",
        R2_FQ="{OUTDIR}/{SAMPLE}/tmp/merged_R2.fq.gz",
    output:
        R1_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/postTrim_R1.fq.gz"),
        R2_FQ=temp("{OUTDIR}/{SAMPLE}/tmp/postTrim_R2.fq.gz"),
        JSON="{OUTDIR}/{SAMPLE}/misc_logs/cutadapt.json",
    params:
        LENGTH_MIN = 12,
        QUALITY_MIN=20,
        OVERLAP=5,
        HOMOPOLYMER_ERROR_RATE=0.2,  # default error rate is 0.1
        POLYA="A" * 100,
        POLYT="T" * 100,
        NEXTERA="CTGTCTCTTATA",  # Nextera sequence
        rcNEXTERA="TATAAGAGACAG",  # Rev Comp of Nextera sequence        
        TSO="AAGCTGGTATCAACGCAGAGTGAATGGG",  # TSO
        rcTSO="CCCATTCACTCTGCGTTGATACCAGCTT",  # rev comp of TSO
        ILMN_UNIVERSAL="AGATCGGAAGAG",  # Illumina Universal
    log:
        log="{OUTDIR}/{SAMPLE}/misc_logs/cutadapt.log",
    resources:
        mem="16G",
    threads: config["CORES"]
    # conda:
    #     f"{workflow.basedir}/envs/cutadapt.yml"
    shell:
        """
        mkdir -p $(dirname {output.JSON})

        cutadapt \
            --minimum-length {params.LENGTH_MIN}:{params.LENGTH_MIN} \
            --quality-cutoff {params.QUALITY_MIN} \
            --overlap {params.OVERLAP} \
            --match-read-wildcards \
            --nextseq-trim=20 \
            -A POLYA_3p="{params.POLYA}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -A POLYT_3p="{params.POLYT}X;max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -B TSO={params.TSO} \
            -B NEXTERA="{params.NEXTERA}" \
            -A ILMN_UNIVERSAL_3p="{params.ILMN_UNIVERSAL}" \
            --pair-filter=any \
            -o {output.R1_FQ} \
            -p {output.R2_FQ} \
            --cores {threads} \
            --json {output.JSON} \
            {input.R1_FQ} {input.R2_FQ} \
        1> {log.log}
        """