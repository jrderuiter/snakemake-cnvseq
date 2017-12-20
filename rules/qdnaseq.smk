from os import path


rule qdnaseq:
    input:
        bams=expand("bam/final/{sample}.bam", sample=get_samples()),
        blacklists=config["references"]["qdnaseq"]["blacklists"]
    output:
        rds="qdnaseq/cgh.rds",
        logratios="qdnaseq/logratios.txt",
        segments="qdnaseq/segmented.txt",
        calls="qdnaseq/calls.txt",
        probs="qdnaseq/probs.txt",
        plots="qdnaseq/plots"
    params:
        genome=config["references"]["qdnaseq"]["genome"],
        organism=config["references"]["qdnaseq"]["organism"],
        bin_size=config["options"]["qdnaseq"]["bin_size"],
        read_length=config["options"]["qdnaseq"]["read_length"],
        normals=get_normals()
    conda:
        path.join(workflow.basedir, "envs/qdnaseq.yaml")
    log:
        "logs/qdnaseq.log"
    script:
        path.join(workflow.basedir, "scripts/qdnaseq.R")


rule qdnaseq_annotate:
    input:
        data="qdnaseq/{datatype}.txt",
        gtf=config["references"]["gtf"]
    output:
        "qdnaseq/{datatype}.ann.txt"
    params:
        script=path.join(workflow.basedir, "scripts/qdnaseq_annotate.py"),
        extra=" ".join(config["rules"]["qdnaseq_annotate"]["extra"])
    conda:
        path.join(workflow.basedir, "envs/genopandas.yaml")
    shell:
        "python {params.script} --input {input.data}"
        " --output {output[0]} --gtf {input.gtf} {params.extra}"
