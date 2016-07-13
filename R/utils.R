data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
    stopifnot(class(df) == "data.frame")
    stopifnot(all(c("start", "end") %in% names(df)))
    stopifnot(any(c("chr", "seqnames") %in% names(df)))
    if("seqnames" %in% names(df))
        names(df)[names(df) == "seqnames"] <- "chr"
    if(!ignoreStrand && "strand" %in% names(df)) {
        if(is.numeric(df$strand)) {
            strand <- ifelse(df$strand == 1, "+", "*")
            strand[df$strand == -1] <- "-"
            df$strand <- strand
        }
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end),
                      strand = df$strand)
    } else {
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end))
    }
    if(keepColumns) {
        dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
                     "DataFrame")
        mcols(gr) <- dt
    }
    names(gr) <- rownames(df)
    gr
}

.checkAssayNames <- function(object, names) {
    nms <- assayNames(object)
    if(!all(names %in% nms))
        return(sprintf("object of class '%s' needs to have assay slots with names '%s'",
                       class(object), paste0(names, collapse = ", ")))
    else
        NULL
}

# TODO: No longer required
setMethod("assays", "BSseq",
          function(x, ..., withDimnames = TRUE) {
              x@assays$field("data")
          })

# TODO: No longer required
setMethod("assayNames", "BSseq",
          function(x, ...) {
              names(x@assays$field("data"))
          })
# stats::plogis() generalised to handle DelayedArray input
# TODO: Herve has added plogis to HDF5Array, but needs a version bump to
#       propagate to the build system
#       Once working, remove this hack
.plogis <- function(x) {
    # TODO: Have asked Hervé what is the correct way to register a delayed
    # op, i.e. is there an officially supported (and exported) method?
    if (is(x, "DelayedArray")) {
        y <- HDF5Array:::register_delayed_op(x, plogis)
    } else {
        y <- plogis(x)
    }
    y
}

# NOTE: This used to be defined in both BSmooth.tstat() and smoothSds(); pulled
#       out to avoid code duplication and renamed to indicate this is an
#       internal helper function
#' Smooth standard deviations using a running mean
#'
#' @param Sds A vector or column matrix (not checked) of standard deviations
#' to be smoothed
#' @param k The window size used by runmean()
#' @param qSd The quantile to be used in thresholding the standard deviations
#'
#' @note This is really a general running mean smoother but is named with
#' 'standard deviation' because that is what it is used on in bsseq
#'
#' @return A vector of length equal to the length of the Sds
.smoothSd <- function(Sds, k, qSd) {
    k0 <- floor(k/2)
    if(all(is.na(Sds))) return(Sds)
    thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
    addSD <- rep(median(Sds, na.rm = TRUE), k0)
    sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
    sSds
}
