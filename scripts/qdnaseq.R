#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  require(docopt)
  require(methods)
})

"
Usage:
    qdnaseq.R [options] <bam_paths>...

Description: Wrapper script for running QDNAseq.

Options:
    --out_dir=<path>              Output directory.
    --bin_size=<int>              Bin size for generating counts. [Default: 50]
    --read_length=<int>           Read length of input data. [Default: 50]
    --genome=<str>                Genome to use. [Default: mm10]
    --organism=<str>              Organism from which data is derived. [Default: Other]
    --blacklists=<path,path...>   BED files describing blacklisted regions.
    --spleens=<str,str...>        Samples to use for spleen correction.

" -> doc

# Retrieve the command-line arguments
opts <- docopt(doc)

# Load remaining libraries. Note that these libraries are not imported until
# after the docopt call to increase the responsiveness of the argument parsing.
suppressMessages({
  library(matrixStats)  # quick fix for missing colMedians import
  library(QDNAseq)
  library(QDNAseq.mm10)
  library(CGHcall)
})

# -------------------------- Helper functions ----------------------------------


setGeneric("normalizeBinsBySpleens",
  function(object, spleens, method="median", force=FALSE)
  standardGeneric("normalizeBinsBySpleens"))


setMethod("normalizeBinsBySpleens", signature=c(object="QDNAseqCopyNumbers"),
          definition=function(object, spleens,
                              method=c("median", "mean", "mode"),
                              force=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument "object":
  if (!force && ("segmented" %in% assayDataElementNames(object)))
    stop("Data has already been segmented. Re-normalizing will ",
          "remove segmentation (and possible calling) results. ",
          "Please specify force=TRUE, if you want this.")

  # Argument "method":
  method <- match.arg(method);

  if ("segmented" %in% assayDataElementNames(object))
    assayDataElement(object, "segmented") <- NULL
  if ("calls" %in% assayDataElementNames(object)) {
    assayDataElement(object, "calls") <- NULL
    assayDataElement(object, "probloss") <- NULL
    assayDataElement(object, "probnorm") <- NULL
    assayDataElement(object, "probgain") <- NULL
    if ("probdloss" %in% assayDataElementNames(object))
      assayDataElement(object, "probdloss") <- NULL
    if ("probamp" %in% assayDataElementNames(object))
      assayDataElement(object, "probamp") <- NULL
  }

  # Extract corrected counts
  copynumber <- assayDataElement(object, "copynumber")

  # Extract annotation data
  fData <- fData(object)

  # Sanity check
  stopifnot(is.matrix(copynumber))

  # Filter
  # condition <- QDNAseq:::binsToUse(object)

  cn_spleens = copynumber[, spleens]

  QDNAseq:::vmsg("Applying ", method, " normalization (spleens)...",
                  appendLF=FALSE)
  if (method == "mean") {
    values <- rowMeans(cn_spleens, na.rm=TRUE)
  } else if (method == "median") {
    values <- rowMedians(cn_spleens, na.rm=TRUE)
  } else if (method == "mode") {
    values <- apply(cn_spleens, MARGIN=1L,
      FUN=function(x) {
        d <- density(x, na.rm=TRUE); d$x[which.max(d$y)]
      })
  }

  # Replace zeroes (for division).
  values[values == 0] <- 1

  copynumber2 <- t(scale(t(copynumber), center=FALSE, scale=values))

  # Assign
  assayDataElement(object, "copynumber") <- copynumber2

  # Not needed anymore
  copynumber2 <- NULL

  object
})

# ----------------------------- Main script ------------------------------------

# Get bin annotations.
binsize <- as.integer(opts$bin_size)
bins <- getBinAnnotations(binSize=binsize, genome=opts$genome)

# Use external blacklists if given.
if (!is.null(opts$blacklists)) {
  blacklists = strsplit(opts$blacklists, ',')[[1]]
  bins$blacklist = calculateBlacklist(bins@data, blacklists)
}

# Perform analysis
readCounts <- binReadCounts(bins, bamfiles=opts$bam_paths)

readCountsFiltered <- applyFilters(
  readCounts, residual=TRUE, blacklist=TRUE, mappability=20, chromosomes="Y")
readCountsFiltered <- estimateCorrection(readCountsFiltered)

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)

if (!is.null(opts$spleens)) {
  spleens = strsplit(opts$spleens, ',')[[1]]
  copyNumbersNormalized = normalizeBinsBySpleens(
    copyNumbersNormalized, spleens=spleens)
}

copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

# Note: standard transformFun is log2.
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

copyNumbersCalled <- callBins(copyNumbersSegmented, organism=opts$organism)
cgh <- makeCgh(copyNumbersCalled)

# Create probability matrix.
probs = matrix(0, dim(cgh)[1], dim(cgh)[2])
colnames(probs) = colnames(cgh)
rownames(probs) = rownames(cgh)

for (i in 1:dim(cgh)[1]){
  for (j in 1:dim(cgh)[2]){
    if (calls(cgh)[i,j] == 2){
      probs[i,j] = probamp(cgh)[i,j]
    } else if (calls(cgh)[i,j] == 1){
      probs[i,j] = probgain(cgh)[i,j]
    } else if (calls(cgh)[i,j] == -1){
      probs[i,j] = probloss(cgh)[i,j]
    } else if (calls(cgh)[i,j] == -2){
      probs[i,j] = probdloss(cgh)[i,j]
      }
  }
}

# Define/create output directories.
calls_dir = file.path(opts$out_dir)
object_dir = file.path(opts$out_dir, 'images')
plots_dir = file.path(opts$out_dir, 'plots')

if (!dir.exists(object_dir)) { dir.create(object_dir, recursive=T) }
if (!dir.exists(calls_dir)) { dir.create(calls_dir, recursive=T) }
if (!dir.exists(plots_dir)) { dir.create(plots_dir, recursive=T) }

# Dump images.
image_path = file.path(object_dir, 'image.rdata')
save.image(image_path)

cgh_obj_path = file.path(object_dir, 'cgh.rds')
saveRDS(cgh, cgh_obj_path)

# Save results.
probs_path = file.path(calls_dir, 'probs.txt')
write.table(probs, probs_path, sep='\t',
            quote=FALSE, col.names=NA)

calls_path = file.path(calls_dir, 'calls.txt')
write.table(calls(cgh), calls_path, sep='\t',
            quote=FALSE, col.names=NA)

segments_path = file.path(calls_dir, 'segments.txt')
write.table(segmented(cgh), segments_path, sep='\t',
            quote=FALSE, col.names=NA)

cgh_path = file.path(calls_dir, 'logratios.txt')
write.table(copynumber(cgh), cgh_path, sep='\t',
            quote=FALSE, col.names=NA)

# Create call plots (pdf).
for (i in 1:length(colnames(copyNumbersCalled))){
  name = colnames(copyNumbersCalled)[i]
  pdf_path = file.path(plots_dir, paste0(name, '.pdf'))

  pdf(pdf_path, width=21, height=8)
  plot(copyNumbersCalled[, i])
  dev.off()
}
