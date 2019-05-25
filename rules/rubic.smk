from os import path


rule rubic_prepare_markers:
    input:
        "qdnaseq/segmented.txt"
    params:
        chromosomes=config["references"]["rubic"]["chromosomes"],
        chromosome_map=config["references"]["rubic"]["chromosome_map"]
    output:
        "rubic/inputs/markers.txt"
    conda:
        path.join(workflow.basedir, "envs", "genopandas.yaml")
    script:
        path.join(workflow.basedir, "scripts", "rubic", "prepare_markers.py")


rule rubic_prepare_segments:
    input:
        "qdnaseq/segmented.txt"
    params:
        chromosomes=config["references"]["rubic"]["chromosomes"],
        chromosome_map=config["references"]["rubic"]["chromosome_map"]
    output:
        "rubic/inputs/segments.txt"
    conda:
        path.join(workflow.basedir, "envs", "genopandas.yaml")
    script:
        path.join(workflow.basedir, "scripts", "rubic", "prepare_segments.py")


rule rubic_prepare_genes:
    output:
        "rubic/inputs/genes.txt"
    params:
        host=config["references"]["rubic"]["biomart_host"],
        dataset=config["references"]["rubic"]["biomart_dataset"],
        chromosomes=config["references"]["rubic"]["chromosomes"]
    conda:
        path.join(workflow.basedir, "envs", "pybiomart.yaml")
    script:
        path.join(workflow.basedir, "scripts", "rubic", "fetch_genes.py")


rule rubic:
    input:
        segments="rubic/inputs/segments.txt",
        markers="rubic/inputs/markers.txt",
        genes="rubic/inputs/genes.txt"
    params:
        samples=get_samples_for_group,
        focal_threshold=config["rules"]["rubic"]["focal_threshold"],
        min_probes=config["rules"]["rubic"]["min_probes"],
        fdr=config["rules"]["rubic"]["fdr"]
    output:
        rds="rubic/results/{group}/rubic.rds",
        gains="rubic/results/{group}/focal_gains.tsv",
        losses="rubic/results/{group}/focal_losses.tsv",
        plots="rubic/results/{group}/plots"
    conda:
        path.join(workflow.basedir, "envs/rubic.yaml")
    log:
        'logs/rubic/{group}.log'
    script:
        path.join(workflow.basedir, "scripts/rubic/rubic.R")
