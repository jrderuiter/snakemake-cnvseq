# Snakemake workflow: cnvseq-qdnaseq

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/cnvseq-qdnaseq.svg?branch=master)](https://travis-ci.org/snakemake-workflows/cnvseq-qdnaseq)

This is a Snakemake workflow for generating CNV profiles/calls from CNV-seq
(shallow WGS) sequencing data using QDNAseq. The workflow is designed to handle
single-end (and optionally multi-lane) sequencing data.

The workflow essentially performs the following steps:

* The input reads are trimmed to remove adapters and/or poor quality base calls
  using cutadapt. Optionally, reads are also trimmed to the read length
  used by the QDNAseq reference.
* The trimmed reads are aligned to the reference genome using bwa aln/samse
  and sorted using samtools. (Bwa samse is used because we generally have read
  lengths shorter than 70bp.)
* Bam files from multiple lanes are merged using picard MergeSamFiles.
* QDNAseq is used (via an R script) to analyse the bam files and generate
  the CNV profiles and calls.

QC statistics are generated using fastqc and samtools stats. The statistics are
summarized in a single overview using multiqc.

**Note that this workflow is still under active development.**

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the
[latest release](https://github.com/snakemake-workflows/cnvseq-qdnaseq/releases).
If you intend to modify and further develop this workflow, fork this
repository. Please consider providing any generally applicable modifications
via a pull request.

In any case, if you use this workflow in a paper, don't forget to give
credits to the authors by citing the URL of this repository and, if available,
its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs by editing the config file
`config.yaml` and the samples file `samples.tsv`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for
further details.

## Authors

* Julian de Ruiter (@jrderuiter)
