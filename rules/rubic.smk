from os import path

import pandas as pd
import pybiomart

from genopandas.ngs.cnv import CnvValueMatrix


rule rubic_prepare_markers:
    input:
        "qdnaseq/segmented.txt"
    params:
        chromosomes=config["references"]["rubic"]["chromosomes"],
        chromosome_map=config["references"]["rubic"]["chromosome_map"]
    output:
        "rubic/inputs/markers.txt"
    run:
        # Read values.
        segmented = CnvValueMatrix.from_csv_condensed(input[0], sep='\t')

        # Map chromosomes (if needed) and subset.
        if params.chromosome_map:
            segmented = segmented.rename_chromosomes(params.chromosome_map)

        segmented = segmented.gloc[params.chromosomes]

        # Convert to markers.
        markers = pd.DataFrame(
            {
                'Name': ['P{}'.format(i + 1) for i in
                         range(segmented.shape[0])],
                'Chromosome': segmented.gloc.chromosome,
                'Position': segmented.gloc.position
            },
            columns=['Name', 'Chromosome', 'Position'])

        markers.to_csv(output[0], sep="\t", index=False)


rule rubic_prepare_segments:
    input:
        "qdnaseq/segmented.txt"
    params:
        chromosomes=config["references"]["rubic"]["chromosomes"],
        chromosome_map=config["references"]["rubic"]["chromosome_map"]
    output:
        "rubic/inputs/segments.txt"
    run:
        # Read segmented values.
        segmented = CnvValueMatrix.from_csv_condensed(input[0], sep='\t')

        # Map chromosomes (if needed).
        if params.chromosome_map:
            segmented = segmented.rename_chromosomes(params.chromosome_map)

        # Convert to segments.
        segments = segmented.as_segments().reset_index()

        # Order by sample/genomic position.
        segments = (
            segments
            .assign(chromosome=lambda df: pd.Categorical(
                df['chromosome'], categories=params.chromosomes))
            .dropna(subset=['chromosome'])
            .sort_values(by=["sample", "chromosome", "start", "end"]))

        # Rename and re-order columns.
        segments = segments.rename(columns={
            'chromosome': 'Chromosome',
            'start': 'Start',
            'end': 'End',
            'value': 'LogRatio',
            'sample': 'Sample'
        })

        segments = segments.reindex(columns=[
            'Sample', 'Chromosome', 'Start', 'End', 'LogRatio'])

        # Write output.
        segments.to_csv(output[0], sep="\t", index=False)


def fetch_genes(host, dataset, chromosomes=None):
    """Fetches gene definition from biomart."""

    dataset = pybiomart.Dataset(name=dataset, host=host)

    if chromosomes is not None:
        filters = {'chromosome_name': chromosomes}
    else:
        filters = {}

    result = dataset.query(
        attributes=[
            'ensembl_gene_id', 'external_gene_name', 'chromosome_name',
            'start_position', 'end_position'
        ],
        filters=filters,
        use_attr_names=True)

    genes = result.rename(columns={
        'ensembl_gene_id': 'ID',
        'external_gene_name': 'Name',
        'chromosome_name': 'Chromosome',
        'start_position': 'Start',
        'end_position': 'End'
    })

    return genes


rule rubic_prepare_genes:
    output:
        "rubic/inputs/genes.txt"
    params:
        host=config["references"]["rubic"]["biomart_host"],
        dataset=config["references"]["rubic"]["biomart_dataset"],
        chromosomes=config["references"]["rubic"]["chromosomes"]
    run:
        genes = fetch_genes(
            host=params.host,
            dataset=params.dataset,
            chromosomes=params.chromosomes)
        genes.to_csv(output[0], sep="\t", index=False)


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
        path.join(workflow.basedir, "scripts/rubic.R")
