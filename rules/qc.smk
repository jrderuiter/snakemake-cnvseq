from os import path


rule multiqc:
    input:
        expand("qc/fastqc/{sample_lane}.{pair}_fastqc.html",
               sample_lane=get_samples_with_lane(), pair=["R1"]),
        expand("qc/samtools_stats/{sample}.txt",
               sample=get_samples()),
        expand("qc/picard_mark_duplicates/{sample}.metrics",
               sample=get_samples())
    output:
        "qc/multiqc_report.html"
    params:
        config["multiqc"]["extra"]
    log:
        "logs/multiqc.log"
    wrapper:
        "0.17.0/bio/multiqc"


rule fastqc:
    input:
        "fastq/trimmed/{sample}.{lane}.{pair}.fastq.gz"
    output:
        html="qc/fastqc/{sample}.{lane}.{pair}_fastqc.html",
        zip="qc/fastqc/{sample}.{lane}.{pair}_fastqc.zip"
    params:
        config["fastqc"]["extra"]
    wrapper:
        "0.17.0/bio/fastqc"


rule samtools_stats:
    input:
        "bam/final/{sample}.bam"
    output:
        "qc/samtools_stats/{sample}.txt"
    wrapper:
        "0.17.0/bio/samtools/stats"
