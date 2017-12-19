import pandas as pd

from genopandas.ngs.cnv import CnvValueMatrix

# Read values.
segmented = CnvValueMatrix.from_csv_condensed(snakemake.input[0], sep='\t')

# Map chromosomes (if needed) and subset.
if snakemake.params.chromosome_map:
    segmented = segmented.rename_chromosomes(snakemake.params.chromosome_map)

segmented = segmented.gloc[snakemake.params.chromosomes]

# Convert to markers.
markers = pd.DataFrame(
    {
        'Name': ['P{}'.format(i + 1) for i in range(segmented.shape[0])],
        'Chromosome': segmented.gloc.chromosome,
        'Position': segmented.gloc.position
    },
    columns=['Name', 'Chromosome', 'Position'])

markers.to_csv(snakemake.output[0], sep="\t", index=False)
