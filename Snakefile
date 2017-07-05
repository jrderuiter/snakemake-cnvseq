import pandas as pd


configfile: "config.yaml"

################################################################################
# Globals                                                                      #
################################################################################

samples = pd.read_csv("samples.tsv", sep="\t")


################################################################################
# Functions                                                                    #
################################################################################

def get_samples():
    """Returns list of all samples."""
    return list(samples["sample"].unique())


def get_normal_samples():
    """Returns list of normal samples."""
    subset = samples.loc[samples["type"] == "normal"]
    return list(subset["sample"].unique())


def get_samples_with_lane():
    """Returns list of all combined lane/sample identifiers."""
    return list((samples["sample"] + "." + samples["lane"]).unique())


def get_sample_lanes(sample):
    """Returns lanes for given sample."""
    subset = samples.loc[samples["sample"] == sample]
    return list(subset["lane"].unique())


################################################################################
# Rules                                                                        #
################################################################################

rule all:
    input:
        expand("bam/final/{sample}.bam", sample=get_samples()),
        expand("bam/final/{sample}.bam.bai", sample=get_samples()),
        "qdnaseq/logratios.txt",
        "qc/multiqc_report.html"

include: "rules/input.smk"
include: "rules/fastq.smk"
include: "rules/alignment.smk"
include: "rules/qdnaseq.smk"
include: "rules/qc.smk"
