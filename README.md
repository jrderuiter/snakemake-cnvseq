# Snakemake workflow: cnvseq-qdnaseq

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![wercker status](https://app.wercker.com/status/ace261cedf02ae669a31189a1363e61d/s/master "wercker status")](https://app.wercker.com/project/byKey/ace261cedf02ae669a31189a1363e61d)

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
summarized into a single report using multiqc.

**Note that this workflow is still under active development.**

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the
[latest release](https://github.com/jrderuiter/snakemake-cnvseq-qdnaseq/releases).
If you intend to modify and further develop this workflow, fork this
repository. Please consider providing any generally applicable modifications
via a pull request.

In any case, if you use this workflow in a paper, don't forget to give
credits to the authors by citing the URL of this repository and, if available,
its DOI (see above).

### Step 2: Install dependencies

To be able to run the workflow, you need to have snakemake and pandas
installed. The various tools (e.g. bwa, samtools) also need to be installed
or can be managed via snakemake using conda (with the --use-conda flag).

### Step 3: Configure workflow

Configure the workflow according to your needs by editing the files
`config.yaml` and `samples.tsv`. Note that fastq file paths can be specified
as local file paths or remote http-based urls (other options can be added
on request).

Optionally, spleen samples can be supplied to normalize for background
variation that is not corrected for by QDNAseq's adjustment for QC-content and
mappability. Spleen samples are be distinguished by labeling them as "spleen"
in the `type` column of `samples.tsv` file. Tumor samples should simply be
labeled "tumor".

### Step 4: Execute workflow

Test your configuration by performing a dry-run using

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment using

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

The workflow can be executed in a different directory using

    snakemake --directory ~/scratch/exome

Note that this directory should contain the appropriate sample and config files.

See the [Snakemake documentation](https://snakemake.readthedocs.io) for
further details.

## Authors

* Julian de Ruiter (@jrderuiter)
