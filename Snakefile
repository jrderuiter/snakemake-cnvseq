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

def all_inputs(wildcards):
    inputs = ["qdnaseq/logratios.txt", "qc/multiqc_report.html"]

    samples = get_samples()
    inputs += expand("bam/final/{sample}.bam", sample=samples)
    inputs += expand("bam/final/{sample}.bam.bai", sample=samples)

    if config["options"]["annotate_qdnaseq"]:
        datatypes = ["calls", "logratios", "probs", "segmented"]
        inputs += expand("qdnaseq/{datatype}.ann.txt", datatype=datatypes)

    return inputs

rule all:
    input: all_inputs

include: "rules/input.smk"
include: "rules/fastq.smk"
include: "rules/alignment.smk"
include: "rules/qdnaseq.smk"
include: "rules/qc.smk"
