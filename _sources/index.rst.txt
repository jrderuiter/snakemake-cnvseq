Snakemake-cnvseq-qdnaseq
========================

|Snakemake| |Wercker|

Snakemake-cnvseq-qdnaseq is a snakemake workflow that generates copy number
profiles/calls from CNV-seq (shallow whole genome sequencing) data. The
workflow is designed to handle single-end and (optionally) multi-lane sequencing
data. The actual copy number profiles and calls are generated using the
QDNAseq R package.

If you use this workflow in a paper, don't forget to give credits
to the authors by citing the URL of this repository and, if available, its
DOI (see above).

.. toctree::
   :maxdepth: 2
   :hidden:

   overview
   installation
   configuration
   usage
   contributing
   authors
   history

.. |Snakemake| image:: https://img.shields.io/badge/snakemake-≥3.13.3-brightgreen.svg
   :target: https://snakemake.bitbucket.io

.. |Wercker| image:: https://app.wercker.com/status/ace261cedf02ae669a31189a1363e61d/s/develop
   :target: https://app.wercker.com/project/byKey/ace261cedf02ae669a31189a1363e61d
