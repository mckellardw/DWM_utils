
#############################################
## qualimap on aligned reads
#############################################
rule qualimapQC:
    input:
        BAM="{OUTDIR}/{SAMPLE}/star/Aligned.sortedByCoord.out.bam",
    output:
        REPORT="{OUTDIR}/{SAMPLE}/qualimap/qualimapReport.html",
    params:
        GENES_GTF=config["GENES_GTF"],
    threads: config["CORES"]
    resources:
        mem="4G"
    shell:
        """
        mkdir -p $(dirname {output.REPORT})

        qualimap rnaseq \
            -bam {input.BAM} \
            -gtf {params.GENES_GTF} \
            --sequencing-protocol strand-specific-forward \
            --sorted \
            --java-mem-size={resources.mem} \
            -outdir $(dirname {output.REPORT}) \
            -outformat html
        """
