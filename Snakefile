import pandas as pd

if not config:
    raise ValueError("A config file must be provided using --configfile")

def _invert_dict(d):
    return dict( (v,k) for k in d for v in d[k] )

_unit_sample_lookup = _invert_dict(config['samples'])


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

def get_sample_for_unit(unit):
    """Returns sample for given unit."""
    return _unit_sample_lookup[unit]

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
    output: touch(".all")

include: "rules/input.smk"
include: "rules/fastq.smk"
include: "rules/alignment.smk"
include: "rules/qdnaseq.smk"
include: "rules/qc.smk"
