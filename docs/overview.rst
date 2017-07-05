Overview
========

The workflow essentially performs the following steps:

* The input reads are trimmed to remove adapters and/or poor quality base
  calls using cutadapt.
* The trimmed reads are aligned to the reference genome using bwa aln.
  (Bwa aln is used because we typically have data with short read lengths.
  Bwa mem support can be added on request.)
* The alignments are sorted and indexed using samtools.
* Bam files from multiple lanes are merged using samtools.
* Picard MarkDuplicates is used to remove optical/PCR duplicates.
* The final alignments are indexed using samtools index.
* QDNAseq is used to generate copy number profiles/calls from the alignments.
  This step also supports normalization using normal samples and
  blacklisting of specific genomic regions.

QC statistics are generated using fastqc, samtools stats and picard
MarkDuplicates. The stats are summarized into a single report using multiqc.

Altogether, this results in the following dependency graph:

.. figure:: images/dag.svg
  :align: center

Note that the QDNAseq step should be configured carefully to reflect the used
genome and read length. Make sure that a bin annotation is available for your
reference genome and read length, otherwise you may need to generate one.
(See the `QDNAseq documentation`_ for more details.) Alternatively, if no
bin annotation is available for your specific read length, you can also trim
your reads to a supported length in the ``cutadapt`` step.

.. _QDNAseq documentation: https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf
