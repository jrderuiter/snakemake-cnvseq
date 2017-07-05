from os import path


rule qdnaseq:
    input:
        expand("bam/final/{sample}.bam", sample=get_samples())
    output:
        rds="qdnaseq/cgh.rds",
        logratios="qdnaseq/logratios.txt",
        segments="qdnaseq/segmented.txt",
        calls="qdnaseq/calls.txt",
        probs="qdnaseq/probs.txt",
        plots="qdnaseq/plots"
    params:
        genome=config["qdnaseq"]["genome"],
        organism=config["qdnaseq"]["organism"],
        bin_size=config["qdnaseq"]["bin_size"],
        read_length=config["qdnaseq"]["read_length"],
        blacklists=config["qdnaseq"]["blacklists"],
        normals=get_normal_samples()
    conda:
        path.join(workflow.basedir, "envs/qdnaseq.yaml")
    log:
        "logs/qdnaseq.log"
    script:
        path.join(workflow.basedir, "scripts/qdnaseq.R")
