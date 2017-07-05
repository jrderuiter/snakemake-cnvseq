Configuration
=============

The pipeline is configured using two config files, ``samples.tsv`` and
``config.yaml``. The ``samples.tsv`` file defines the input samples,
together with the paths to the respective source files. The ``config.yaml``
file provides detailed configuration options for the different
steps of the pipeline.

Sample definition
-----------------

The sample definition file is a tab-separated file that lists the samples
that are to be processed by the pipeline. Each row represents a single set of
(paired-end) fastq files for a given sample on a given lane. As such, samples
that have been sequenced on multiple lanes with typically span multiple rows
that share the same sample ID.

For example, a single sample sequenced over two lanes would be described
as follows:

.. literalinclude:: ../samples.tsv

The ``fastq1`` column should contain paths to the input files, which are
expected to be fastq files. These paths can either be provided as local
relative/absolute paths, or as remote http/ftp urls. Note that relative file
paths are taken relative to the input directory defined in the configuration
file (see below for more details), if specified.

The ``lane`` column is used to distinguish sequencing data from the same sample
that has been sequenced in different lanes. This column can be filled with
dummy values (i.e. L999) if lane information is not available and samples
were sequenced on a single lane.

The ``type`` column is used to indicate whether a given sample represents a
tumor or a normal samples. If provided, normal samples are used to (soft)
correct for biases that are not corrected by QDNAseq's correction for
mappability and QC content.

Pipeline configuration
----------------------

The individual steps of the pipeline are configured using the ``config.yaml``
file. This config file contains two different sections, which define
configurations for the inputs and for each specific rule, respectively.

Input options
~~~~~~~~~~~~~

The input section defines several options regarding the handling of the
input files:

.. literalinclude:: ../config.yaml
    :lines: 1-12

Here, ``dir`` is an optional value that defines the directory containing
the input files. If given, file paths provided in ``samples.tsv`` are
sought relative to this directory. Its value is ignored if http/ftp urls
are used for the inputs.

The ``ftp`` section defines the username/password to use when downloading
samples over an ftp connection. These values can be omitted when downloading
files from an anonymous ftp server.

Rule options
~~~~~~~~~~~~

The rule section provides detailed configuration options for the different
rule of the workflow. In general, each rule has a set of configurable options
under the same name as the step itself. The options themselves are specific for
each step and the corresponding tool, but each step typically has an ``extra``
option, which allows you to pass arbitrary arguments to the underlying tool.

.. literalinclude:: ../config.yaml
    :lines: 15-

The final rule, ``qdnaseq``, executes an R script which runs the various
QDNAseq functions. This script requires a suitable QDNAseq bin annotation to be
installed, which matches the given genome and read length. See the
`QDNAseq documentation`_ for more details on how to obtain or generate such a
reference. If no bin annotation is available for your given read length, you
can also trim your reads to the required length in the ``cutadapt`` step (using
the extra parameter ``--length``).

.. _QDNAseq documentation: https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf
