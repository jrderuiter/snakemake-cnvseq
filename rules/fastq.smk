from os import path

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider


HTTP = HTTPRemoteProvider()


def cutadapt_input(wildcards):
    # Lookup input paths.
    key = (wildcards.sample, wildcards.lane)
    row = samples.set_index(['sample', 'lane']).loc[key]

    input_ = row['fastq']

    # Wrap as URL if needed.
    if input_.startswith('http'):
        input_ = HTTP.remote(input_, keep_local=True)

    return input_


rule cutadapt:
    input:
        cutadapt_input
    output:
        fastq='fastq/trimmed/{sample}.{lane}.R1.fastq.gz',
        qc='qc/cutadapt/{sample}.{lane}.qc.txt'
    params:
        config['cutadapt']['extra']
    log:
        'logs/cutadapt/{sample}.{lane}.log'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/cutadapt/se')
