def cutadapt_extra(wildcards):

    extra = config["rules"]["cutadapt"]["extra"]
    extra += ['--length {}'.format(
        config["options"]["qdnaseq"]["read_length"])]

    return " ".join(extra)


rule cutadapt:
    input:
        "fastq/raw/{unit}.R1.fastq.gz"
    output:
        fastq=temp("fastq/trimmed/{unit}.R1.fastq.gz"),
        qc="qc/cutadapt/{unit}.txt"
    params:
        cutadapt_extra
    threads:
        config["rules"]["cutadapt"]["threads"]
    log:
        "logs/cutadapt/{unit}.log"
    wrapper:
        path.join(workflow.basedir, "wrappers", "cutadapt", "se")
