import pybiomart

dataset = pybiomart.Dataset(
    name=snakemake.params.dataset, host=snakemake.params.host)

if snakemake.params.chromosomes is not None:
    filters = {'chromosome_name': snakemake.params.chromosomes}
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

genes.to_csv(snakemake.output[0], sep="\t", index=False)
