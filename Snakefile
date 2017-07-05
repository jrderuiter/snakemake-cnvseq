import pandas as pd


configfile: 'config.yaml'

################################################################################
# Globals                                                                      #
################################################################################

samples = pd.read_csv('samples.tsv', sep='\t')


################################################################################
# Functions                                                                    #
################################################################################

def get_samples():
    return list(samples['sample'].unique())

def get_spleen_samples():
    subset = samples.loc[samples['type'] == 'spleen']
    return list(subset['sample'].unique())

def get_samples_with_lane():
    return list((samples['sample'] + '.' + samples['lane']).unique())

def get_sample_lanes(sample):
    subset = samples.loc[samples['sample'] == sample]
    return list(subset['lane'].unique())


################################################################################
# Rules                                                                        #
################################################################################

rule all:
    input:
        'qdnaseq/logratios.txt',
        'qc/multiqc_report.html'

include: "rules/input.smk"
include: "rules/fastq.smk"
include: "rules/alignment.smk"
include: "rules/qdnaseq.smk"
include: "rules/qc.smk"
