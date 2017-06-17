from os import path


rule bwa_aln:
    input:
        "fastq/trimmed/{sample}.{lane}.{pair}.fastq.gz"
    output:
        "bam/aligned/{sample}.{lane}.{pair}.sai"
    params:
        index=config['bwa_aln']['index'],
        extra=config['bwa_aln']['extra']
    log:
        "logs/bwa_aln/{sample}.{lane}.{pair}.log"
    threads:
        config['bwa_aln']['threads']
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/bwa/aln')
        #"master/bio/bwa/aln"

rule bwa_samse:
    input:
        fastq="fastq/trimmed/{sample}.{lane}.R1.fastq.gz",
        sai="bam/aligned/{sample}.{lane}.R1.sai"
    output:
        "bam/aligned/{sample}.{lane}.bam"
    params:
        index=config['bwa_samse']['index'],
        extra=config['bwa_samse']['extra'],
        sort="samtools",
        sort_order="coordinate",
        sort_extra=config['bwa_samse']['sort_extra']
    log:
        "logs/bwa_samse/{sample}.{lane}.log"
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/bwa/samse')
        #"master/bio/bwa/aln"


def merge_inputs(wildcards):
    lanes = get_sample_lanes(wildcards.sample)

    file_paths = ['bam/aligned/{}.{}.bam'.format(wildcards.sample, lane)
                  for lane in lanes]

    return file_paths


rule picard_merge_bam:
    input:
        merge_inputs
    output:
        'bam/merged/{sample}.bam'
    params:
        config['picard_merge_bam']['extra']
    log:
        'logs/picard_merge_bam/{sample}.log'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/picard/mergesamfiles')


rule picard_mark_duplicates:
    input:
        'bam/merged/{sample}.bam'
    output:
        bam='bam/deduped/{sample}.bam',
        metrics='bam/deduped/{sample}.metrics'
    params:
        config['picard_mark_duplicates']['extra']
    log:
        'logs/picard_mark_duplicates/{sample}.log'
    wrapper:
        '0.15.4/bio/picard/markduplicates'


rule samtools_index:
    input:
        'bam/deduped/{sample}.bam'
    output:
        'bam/deduped/{sample}.bam.bai'
    wrapper:
        "0.15.4/bio/samtools/index"
