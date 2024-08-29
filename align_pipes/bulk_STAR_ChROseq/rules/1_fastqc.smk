


rule fastqc_all:
    input:
        FQ="{OUTDIR}/{SAMPLE}/tmp/merged_{READ}.fq.gz",
    output:
        fastqcDir=directory("{OUTDIR}/{SAMPLE}/fastqc/{STEP}_{READ}"),
    threads: config["CORES"]
    shell:
        """
        mkdir -p {output.fastqcDir}

        fastqc \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            {input.FQ}
        """

