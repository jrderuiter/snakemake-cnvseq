suppressMessages({
    require(methods)
    require(data.table)
    require(RUBIC)
})

# Redirect output to log file.
if (length(snakemake@log) > 0) {
    log_file = con <- file(snakemake@log[[1]])
    sink(log_file)
    sink(log_file, type="message")
}

# Read input files.
segments = fread(snakemake@input$segments)
markers = fread(snakemake@input$markers)
genes = fread(snakemake@input$genes)

# Extract extra options.
if (!is.null(snakemake@params$samples)) {
    samples = snakemake@params$samples
} else {
    samples = unique(segments$Sample)
}

if (!is.null(snakemake@params$focal_threshold)) {
    focal_threshold = as.numeric(snakemake@params$focal_threshold)
} else {
    focal_threshold = 1e+08
}

if (!is.null(snakemake@params$min_probes)) {
    min_probes = as.numeric(snakemake@params$min_probes)
} else {
    min_probes = 4
}

if (!is.null(snakemake@params$fdr)) {
    fdr = as.numeric(snakemake@params$fdr)
} else {
    fdr = 0.25
}

# Run RUBIC.
rbc = rubic(
    fdr=fdr,
    seg.cna=segments,
    markers=markers,
    genes=genes,
    samples=samples,
    focal.threshold=focal_threshold,
    min.probes=min_probes)

rbc$estimate.parameters()
rbc$call.events()
rbc$call.focal.events()

# Write outputs.
if (!is.null(snakemake@output$rds)) {
    rbc$save(snakemake@output$rds)
}

if (!is.null(snakemake@output$gains)) {
    rbc$save.focal.gains(snakemake@output$gains)
}

if (!is.null(snakemake@output$losses)) {
    rbc$save.focal.losses(snakemake@output$losses)
}

if (!is.null(snakemake@output$plots)) {
    rbc$save.plots(snakemake@output$plots)
}
