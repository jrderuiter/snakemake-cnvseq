from os import path


rule bwa_aln:
    input:
        "fastq/trimmed/{unit}.{pair}.fastq.gz"
    output:
        temp("bam/aligned/{unit}.{pair}.sai")
    params:
        index=config["references"]["bwa_index"],
        extra=" ".join(config["rules"]["bwa_aln"]["extra"])
    log:
        "logs/bwa_aln/{unit}.{pair}.log"
    threads:
        config["rules"]["bwa_aln"]["threads"]
    wrapper:
        "0.17.0/bio/bwa/aln"


def samse_extra(wildcards):
    """Generates bwa samse extra arguments."""

    extra = list(config["rules"]["bwa_samse"]["extra"])

    readgroup_str = ('\"@RG\\tID:{unit}\\tSM:{sample}\\t'
                     'LB:{sample}\\tPU:{unit}\\t'
                     'PL:{platform}\\tCN:{centre}\"')

    readgroup_str = readgroup_str.format(
        sample=get_sample_for_unit(wildcards.unit),
        unit=wildcards.unit,
        platform=config["options"]["readgroup"]["platform"],
        centre=config["options"]["readgroup"]["centre"])

    extra += ['-r ' + readgroup_str]

    return " ".join(extra)


rule bwa_samse:
    input:
        fastq="fastq/trimmed/{unit}.R1.fastq.gz",
        sai="bam/aligned/{unit}.R1.sai"
    output:
        temp("bam/aligned/{unit}.bam")
    params:
        index=config["references"]["bwa_index"],
        extra=lambda wc: samse_extra(wc),
        sort="samtools",
        sort_order="coordinate",
        sort_extra=" ".join(config["rules"]["bwa_samse"]["sort_extra"])
    log:
        "logs/bwa_samse/{unit}.log"
    wrapper:
        "0.17.0/bio/bwa/samse"


def merge_inputs(wildcards):
    units = get_sample_units(wildcards.sample)
    return ["bam/aligned/{unit}.bam".format(unit=unit)
            for unit in units]


rule samtools_merge:
    input:
        merge_inputs
    output:
        temp("bam/merged/{sample}.bam")
    params:
        " ".join(config["rules"]["samtools_merge"]["extra"])
    threads:
        config["rules"]["samtools_merge"]["threads"]
    wrapper:
        "0.17.0/bio/samtools/merge"


rule picard_mark_duplicates:
    input:
        "bam/merged/{sample}.bam"
    output:
        bam="bam/final/{sample}.bam",
        metrics="qc/picard_mark_duplicates/{sample}.metrics"
    params:
        " ".join(config["rules"]["picard_mark_duplicates"]["extra"])
    log:
        "logs/picard_mark_duplicates/{sample}.log"
    wrapper:
        "0.17.0/bio/picard/markduplicates"


rule samtools_index:
    input:
        "bam/final/{sample}.bam"
    output:
        "bam/final/{sample}.bam.bai"
    wrapper:
        "0.17.0/bio/samtools/index"
