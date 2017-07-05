rule cutadapt:
    input:
        "fastq/raw/{sample}.{lane}.R1.fastq.gz"
    output:
        fastq=temp("fastq/trimmed/{sample}.{lane}.R1.fastq.gz"),
        qc="qc/cutadapt/{sample}.{lane}.txt"
    params:
        config["cutadapt"]["extra"]
    log:
        "logs/cutadapt/{sample}.{lane}.log"
    wrapper:
        "0.17.0/bio/cutadapt/se"
