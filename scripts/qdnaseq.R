suppressMessages({
    library(matrixStats)  # quick fix for missing colMedians import
    library(QDNAseq)
    library(QDNAseq.mm10)
    library(CGHcall)
})

# Redirect output to log file.
if (!is.null(snakemake@log)) {
    log_file = con <- file(snakemake@log[[1]])
    sink(log_file)
    sink(log_file, type="message")
}


################################################################################
# Functions                                                                    #
################################################################################

setGeneric("normalizeBinsBySpleens",
    function(object, spleens, method="median", force=FALSE)
    standardGeneric("normalizeBinsBySpleens"))


setMethod("normalizeBinsBySpleens",
            signature=c(object="QDNAseqCopyNumbers"),
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
            }
        )
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


probability_matrix <- function(cgh) {
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

    probs
}


################################################################################
# Script                                                                       #
################################################################################

bins <- getBinAnnotations(
    binSize=as.integer(snakemake@params$bin_size),
    genome=snakemake@params$genome)

# Use external blacklists if given.
if (!is.null(snakemake@params$blacklists)) {
    bins$blacklist <- calculateBlacklist(
        bins@data, snakemake@params$blacklists)
}

# Perform analysis
readCounts <- binReadCounts(bins, bamfiles=as.character(snakemake@input))

readCountsFiltered <- applyFilters(
    readCounts, residual=TRUE, blacklist=TRUE,
    mappability=20, chromosomes="Y")
readCountsFiltered <- estimateCorrection(readCountsFiltered)

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)

if (!is.null(snakemake@params$spleens)) {
    copyNumbersNormalized = normalizeBinsBySpleens(
        copyNumbersNormalized, spleens=snakemake@params$spleens)
}

copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

# Note: standard transformFun is log2.
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

copyNumbersCalled <- callBins(
    copyNumbersSegmented, organism=snakemake@params$organism)
cgh <- makeCgh(copyNumbersCalled)

# Save results.
if (!is.null(snakemake@output$rds)) {
    saveRDS(cgh, snakemake@output$rds)
}

if (!is.null(snakemake@output$logratios)) {
    write.table(copynumber(cgh), snakemake@output$logratios,
                sep='\t', quote=FALSE, col.names=NA)
}

if (!is.null(snakemake@output$segments)) {
    write.table(segmented(cgh), snakemake@output$segments,
                sep='\t', quote=FALSE, col.names=NA)
}


if (!is.null(snakemake@output$calls)) {
    write.table(calls(cgh), snakemake@output$calls,
                sep='\t', quote=FALSE, col.names=NA)
}

if (!is.null(snakemake@output$probs)) {
    write.table(probability_matrix(cgh), snakemake@output$probs,
                sep='\t', quote=FALSE, col.names=NA)
}

# Create call plots (pdf).
if (!is.null(snakemake@output$plots)) {
    if (!dir.exists(snakemake@output$plots)) {
        dir.create(snakemake@output$plots, recursive=T)
    }

    for (i in 1:length(colnames(copyNumbersCalled))) {
        name = colnames(copyNumbersCalled)[i]
        pdf_path = file.path(snakemake@output$plots, paste0(name, '.pdf'))

        pdf(pdf_path, width=21, height=8)
        plot(copyNumbersCalled[, i])
        dev.off()
    }
}
