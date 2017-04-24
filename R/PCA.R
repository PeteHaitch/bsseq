# TODO: Standardise argument naming style (arg.name, arg_name, argName, ArgName)

# Helper function to construct bins containing a minimum number of loci with a
# minimum width in bp
.makeBins <- function(gr, min_loci, min_width) {
    stopifnot(!is.unsorted(gr))
    val <- .Call("makeBinEndpoints",
                 ss = start(gr),
                 max.bins = ceiling(end(gr[length(gr)]) / min_width),
                 min_loci = min_loci,
                 min_width = min_width)

    # Sanity check
    stopifnot(val[[1]][1] != 0)

    # NOTE: Need to remove trailing zeros (which are padding created by
    #       makeBinEndpoints
    val <- lapply(val, function(v) v[v != 0])
    GRanges(seqnames = unique(seqnames(gr)),
            ranges = IRanges(start = val[[1]],
                             end = val[[2]]))
}

# Helper function to compute average methylation level in bins along genome
.methPerBin <- function(BSseq, min_width, min_loci, type, mc.cores, verbose) {
    if (is.unsorted(BSseq)) {
        # TODO: Error or message() + sort()?
        message("Sorting BSseq")
        BSseq <- sort(BSseq)
    }

    # Make bins containing at least min_loci and at least min_width wide
    if (verbose) {
        message("Constructing bins with 'min_loci' = ", min_loci,
                " and 'min_width' = ", min_width)
    }
    grl <- split(rowRanges(BSseq), seqnames(BSseq))
    bins <- unlist(endoapply(grl, .makeBins,
                             min_loci = min_loci,
                             min_width = min_width))

    # Compute avereage methylation level in each bin
    if (verbose) {
        message("Computing average methylation in bins")
    }
    # TODO: This bit is going to be very slow, ~22 mins for BSseq with 7.8
    #       million rows and 24 columns (CH loci on chr22); perhaps it can be
    #       parallelised by seqlevel or sample
    # TODO: How big will the output be? ~18 MB (35000 elements) for BSseq with
    #       7.8 million rows and 24 columns (CH loci on chr22)
    getMeth(BSseq, regions = bins, type = type, what = "perRegion")
}

# Internal helper function for plotPCA,BSseq-method
# TODO: Worth having a max binwidth (e.g., to avoid spanning centromere)?
.PCA <- function(BSseq, min_width, min_loci, min_cov, type, ntop, mc.cores,
                 verbose) {
    type <- match.arg(type)
    # TODO: This means the min_cov requirement must hold for all samples; is
    #       this what we want?
    if (verbose) {
        message("Subsetting by 'min_cov' = ", min_cov)
    }
    if (min_cov > 0) {
        BSseq <- BSseq[rowSums(getCoverage(BSseq) >= min_cov) == ncol(BSseq), ]
    }

    # NOTE: prcomp requires matrix to have samples as rows and features as
    #       columns
    binned_meth <- do.call(cbind,
                           .methPerBin(BSseq = BSseq,
                                       min_width = min_width,
                                       min_loci = min_loci,
                                       type = type,
                                       mc.cores = mc.cores,
                                       verbose = verbose))

    message("Computing PCA ")
    if (!is.infinite(ntop)) {
        rv <- colVars(binned_meth)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        binned_meth <- binned_meth[, select, drop = FALSE]
    }
    prcomp(binned_meth)
}

# Internal helper function for plotPCA,BSseq-method
# TODO: How to specify what plot is annotated by (e.g., colour, shape, etc.)
.plotPCA <- function(pca, ...) {

}

# TODO: How to specify what plot is annotated by (e.g., colour, shape, etc.)
setMethod("plotPCA", "BSseq",
          function(object, min_cov = 0L, min_width = 1000, min_loci = 10L,
                   type = c("smooth", "raw"), ntop = Inf, labels = NULL,
                   pch = NULL, cex = 1, dim.plot = c(1, 2),
                   ndim = max(dim.plot), xlab = NULL, ylab = NULL, plot = TRUE,
                   mc.cores = 1, verbose = FALSE) {
              pca <- .PCA(BSseq = BSseq,
                          min_width = min_width,
                          min_loci = min_loci,
                          min_cov = min_cov,
                          type = match.arg(type),
                          ntop = ntop,
                          mc.cores = mc.cores,
                          verbose = verbose)
              # TODO: Write .plotPCA() based on limma::plotMDS() and
              #       DESeq2::plotPCA()
              .plotPCA(pca, ...)
          })
