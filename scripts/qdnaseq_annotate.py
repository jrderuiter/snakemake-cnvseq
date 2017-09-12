#!/usr/bin/env python

import argparse
import logging

from genopandas import GenomicDataFrame, RegionMatrix

logging.basicConfig(
    format='[%(asctime)-15s]  %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def main():
    """Main function."""

    args = parse_args()

    # Read matrix.
    logging.info('Reading matrix values')
    region_values = RegionMatrix.from_csv(
        args.input, sep='\t', expand_index=True)

    # Annotate with genes.
    logging.info('Reading gene annotation')
    genes = GenomicDataFrame.from_gtf(
        args.gtf, filter_=lambda rec: rec.feature == 'gene')

    logging.info('Annotating regions (may take a long time)')
    gene_values = region_values.annotate(genes, id_='gene_name')

    # Write output.
    logging.info('Writing outputs')
    gene_values.values.to_csv(args.output, sep='\t', index=True)


def parse_args():
    """Parses commandline args."""

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--gtf', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    main()
