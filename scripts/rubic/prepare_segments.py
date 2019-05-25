import pandas as pd

from genopandas.ngs.cnv import CnvValueMatrix

# Read segmented values.
segmented = CnvValueMatrix.from_csv_condensed(snakemake.input[0], sep='\t')

# Map chromosomes (if needed).
if snakemake.params.chromosome_map:
    segmented = segmented.rename_chromosomes(snakemake.params.chromosome_map)

# Convert to segments.
segments = segmented.as_segments().reset_index()

# Order by sample/genomic position.
segments = (
    segments.assign(
        chromosome=lambda df: pd.Categorical(
            df['chromosome'], categories=snakemake.params.chromosomes))
    .dropna(subset=['chromosome']).sort_values(
            by=["sample", "chromosome", "start", "end"]))  # yapf: disable

# Rename and re-order columns.
segments = segments.rename(columns={
    'chromosome': 'Chromosome',
    'start': 'Start',
    'end': 'End',
    'value': 'LogRatio',
    'sample': 'Sample'
})

segments = segments.reindex(
    columns=['Sample', 'Chromosome', 'Start', 'End', 'LogRatio'])

# Write output.
segments.to_csv(snakemake.output[0], sep="\t", index=False)
