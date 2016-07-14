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

# Internal helper used by .combineMatrixLike (TODO: and .combineListMatrixList?)
# While designed for DelayedMatrix objects, it will also work for matrix
# objects (although it is not the most efficient way to do so in this case
# since it (A) cbinds() vectors to form the final matrix-like object, and
# (B) returns a HDF5-backed DelayedMatrix object rather than a matrix object)
.combineDelayedMatrix <- function(m, idx, nrow, fill) {
    # TODO: Could consider using mclapply(), but (A) is it worth it? and
    #       (B) would need to be careful in case .combineMatrixList() is
    #       iteself called from a function using mclapply().
    M <- lapply(seq_len(ncol(m)), function(j) {
        mm <- matrix(fill, nrow = nrow, ncol = 1L)
        mm[idx, 1L] <- as.array(m[, j, drop = FALSE])
        HDF5Array(mm)
    })
    if (length(M) > 1) {
        M <- do.call(cbind, M)
    } else {
        M <- M[[1]]
    }
    M
}

#' Combine two matrix-like objects
#'
#' A helper function used by combine (TODO: And combineList?).
#' Combines two matrix-like objects (x, y) into a new matrix-like object with
#' given dimensions (nrow, nccol) by merging according to row indices
#' (idx_x, idx_y) with columns of x to the left of columns of y.
#' The matrix-like object can be a matrix or a DelayedMatrix (from HDF5Array
#' package). The return value is a matrix if both x and y are' matrix objects
#' (i.e. realised in memory), and is otherwise a HDF5Array object (i.e.
#' realised on disk) [note that this means that if the x or y is a
#' DelayedMatrix with an in-memory @seed the result is still written to disk as
#' a HDF5-backed DelayedMatrix].
#' @param x,y A matrix-like object
#' @param idx_x,idx_y A vector of row indices for elements of x (resp. y) in
#' the returned matrix-like object
#' @param nrow,ncol The dimensions of the returned matrix-like object
#' @param fill The value to be used for filling in elements of the returned
#' matrix-like object where a row is not found in one of x or y
.combineMatrixLike <- function(x, y, idx_x, idx_y, nrow, ncol, fill = 0L) {
    if (is.matrix(x) && is.matrix(y)) {
        z <- matrix(fill, nrow = nrow, ncol = ncol)
        z[idx_x, seq_len(ncol(x))] <- x
        z[idx_y, ncol(x) + seq_len(ncol(y))] <- y
    } else if (is(x, "DelayedMatrix") || is(y, "DelayedMatrix")) {
        # NOTE: The strategy for combining x and y when one or both is a
        #       DelayedMatrix object is designed to avoid realising in memory
        #       either object in its entirety at any one time. Instead, we
        #       realise columns of x and y, perform the combine, write the
        #       combined column to disk, and finally cbind() these columns and
        #       write as a new HDF5Matrix.
        z_x <- .combineDelayedMatrix(x, idx_x, nrow, fill)
        z_y <- .combineDelayedMatrix(y, idx_y, nrow, fill)
        # NOTE: z is a DelayedMatrix with two HDF5-backed seeds
        z <- cbind(z_x, z_y)

    } else {
        stop("Cannot combine objects with classes '", class(x),
             "' and '", class(y), "'")
    }
    z
}

.combineListMatrixLike <- function(matrix_list, idx_list, nrow, ncol,
                                   fill = NA_real_, hdf5 = FALSE) {
    cl <- vapply(matrix_list, class, character(1L))
    if (all(cl == "matrix") && !hdf5) {
        z <- matrix(fill, nrow, ncol)
        j0 <- 0
        for (j in seq_along(matrix_list)) {
            i <- idx_list[[j]]
            jj <- j0 + seq_len(ncol(matrix_list[[j]]))
            z[i, jj] <- matrix_list[[j]]
            j0 <- j0 + ncol(matrix_list[[j]])
        }
        if (hdf5) {
            z <- HDF5Array(z)
        }
    } else if (any(cl == "DelayedMatrix") || hdf5) {
        z <- do.call(cbind, mapply(function(m, idx) {
            .combineDelayedMatrix(m, idx, nrow, fill)
        }, m = matrix_list, idx = idx_list))
    } else {
        stop("Cannot combine list of objects with classes: ", paste0(cl, sep = ", "))
    }
    z
}
