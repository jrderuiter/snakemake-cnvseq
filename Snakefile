import pandas as pd


################################################################################
# Functions                                                                    #
################################################################################

def get_samples():
    """Returns list of all samples."""
    return list(config["samples"].keys())

def get_units():
    """Returns list of units."""
    return list(config["units"].keys())

def get_sample_units(sample):
    """Returns lanes for given sample."""
    return config["samples"][sample]

def get_normals():
    """Returns list of normal samples."""
    return config["normals"]


################################################################################
# Rules                                                                        #
################################################################################

def all_inputs(wildcards):
    inputs = ["qdnaseq/logratios.txt", "qc/multiqc_report.html"]

    samples = get_samples()
    inputs += expand("bam/final/{sample}.bam", sample=samples)
    inputs += expand("bam/final/{sample}.bam.bai", sample=samples)

    if config["options"]["qdnaseq"]["annotate"]:
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
