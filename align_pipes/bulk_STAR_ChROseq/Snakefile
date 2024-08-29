########################################################################################################
# bulk_STAR_ChROseq
#   Snakemake workflow to use STAR to align and quantify bulk nascent RNAseq datasets
#   v1.0
#   Written by David McKellar
########################################################################################################

import pdb
import pandas as pd
import glob


########################################################################################################
# Config file
########################################################################################################
configfile: "config.yaml"


########################################################################################################
# Directories and locations
########################################################################################################
TMPDIR = config["TMPDIR"]
OUTDIR = config["OUTDIR"]

########################################################################################################
# Variables and references
########################################################################################################
SAMPLES = list(pd.read_csv(config["SAMPLE_SHEET"])["sampleID"])

R1_FQS = dict(zip(SAMPLES, list(pd.read_csv(config["SAMPLE_SHEET"])["fastq_R1"])))

R2_FQS = dict(zip(SAMPLES, list(pd.read_csv(config["SAMPLE_SHEET"])["fastq_R2"])))
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET']

# SRR = dict(zip(SAMPLES, list(pd.read_csv(config['SAMPLE_SHEET'])['SRR'])))

STAR_REF = config["STAR_REF"]

### Wildcard constraints ###############################################################
wildcard_constraints:
    OUTDIR = config["OUTDIR"],
    SAMPLE = "[A-Za-z0-9_-]+"

########################################################################################################
rule all:
    input:
        expand(
            "{OUTDIR}/{SAMPLE}/fastqc/{STEP}_{READ}",
            OUTDIR=config["OUTDIR"],
            SAMPLE=SAMPLES,
            STEP=["preTrim","postTrim","unmapped"],
            READ=["R1", "R2"]
        ),  # fastQC results
        expand(
            "{OUTDIR}/{SAMPLE}/misc_logs/{FILE}",
            OUTDIR=config["OUTDIR"],
            SAMPLE=SAMPLES,
            FILE=["cutadapt.json"]
        ),  #
        expand(
            "{OUTDIR}/{SAMPLE}/star/{FILE}",
            OUTDIR=config["OUTDIR"],
            SAMPLE=SAMPLES,
            FILE=["Aligned.sortedByCoord.out.bam", "Aligned.sortedByCoord.out.bam.bai", "ReadsPerGene.out.tab"]
        ),  # STAR output(s)
        expand(
            "{OUTDIR}/{SAMPLE}/qualimap/qualimapReport.html",
            OUTDIR=config["OUTDIR"],
            SAMPLE=SAMPLES,
        ),  # alignment QC qith qualimap
        expand(
            "{OUTDIR}/{SAMPLE}/star/Aligned.sortedByCoord.out_plus.bw",
            OUTDIR=config["OUTDIR"],
            SAMPLE=SAMPLES,
        ),  # strand-split bigWigs


#############################################
## Import rules:
#############################################

include: "rules/0_utils.smk"
include: "rules/1_fastqc.smk"
include: "rules/2_trimming.smk"
include: "rules/3a_star.smk"
# include: "rules/3b_bowtie.smk" #TODO- write
include: "rules/4_qc.smk"
include: "rules/5_post.smk"
# include: "rules/X_blast.smk" #TODO- update