#!/usr/bin/env python

import argparse

from genopandas import GenomicDataFrame, RegionMatrix


def main():
    """Main function."""

    args = parse_args()

    # Read matrix.
    region_values = RegionMatrix.from_csv(
        args.input, sep='\t', expand_index=True)

    # Annotate with genes.
    genes = GenomicDataFrame.from_gtf(
        args.gtf, filter=lambda rec: rec.feature == 'gene')
    gene_values = region_values.annotate(genes, id_='gene_name')

    # Write output.
    gene_values.values.to_csv(args.output, sep='\t', index=True)


def parse_args():
    """Parses commandline args."""

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--gtf', required=True)

    return parser.parse_args()
