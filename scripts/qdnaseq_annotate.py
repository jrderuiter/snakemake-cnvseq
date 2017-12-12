#!/usr/bin/env python

import argparse
import logging

from genopandas import GenomicDataFrame, GenomicMatrix

logging.basicConfig(
    format='[%(asctime)-15s]  %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def main():
    """Main function."""

    args = parse_args()

    # Read matrix.
    logging.info('Reading matrix values')
    region_values = GenomicMatrix.from_csv_condensed(
        args.input,
        sep='\t',
        is_one_based=args.one_based,
        is_inclusive=args.inclusive)

    # Impute if needed.
    if args.impute:
        logging.info('Imputing values using neighboring bins')
        region_values = region_values.impute(
            window=args.impute_window, min_probes=args.impute_probes)

    # Annotate with genes.
    logging.info('Reading gene annotation')
    genes = GenomicDataFrame.from_gtf(
        args.gtf, filter_=lambda rec: rec.feature == args.feature_type)

    logging.info('Annotating regions (may take a long time)')
    gene_values = region_values.annotate(genes, feature_id=args.feature_id)

    # Write output.
    logging.info('Writing outputs')
    gene_values.values.to_csv(args.output, sep='\t', index=True)


def parse_args():
    """Parses commandline args."""

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=True)
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--output', required=True)

    parser.add_argument('--impute', default=False, action='store_true')
    parser.add_argument('--impute_window', default=11, type=int)
    parser.add_argument('--impute_probes', default=5, type=int)

    parser.add_argument('--feature_type', default='gene')
    parser.add_argument('--feature_id', default='gene_name')

    parser.add_argument('--one_based', default=False, action='store_true')
    parser.add_argument('--inclusive', default=False, action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()
