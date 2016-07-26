# Internal helper used by .combineListMatrixList
# While designed for DelayedMatrix objects, it will also work for matrix
# objects (although it is not the most efficient way to do so in this case
# since it (A) cbinds() vectors to form the final matrix-like object, and
# (B) returns a HDF5-backed DelayedMatrix object rather than a matrix object)
.combineDelayedMatrix <- function(m, idx, nrow, fill) {
    # TODO: Could consider using mclapply(), but (A) is it worth it? and
    #       (B) would need to be careful in case .combineMatrixList() is
    #       itself called from a function using mclapply()
    M <- lapply(seq_len(ncol(m)), function(j) {
        mm <- matrix(fill, nrow = nrow, ncol = 1L)
        mm[idx, 1L] <- as.array(m[, j, drop = FALSE])
        .safeHDF5Array(mm, "BSseq.", paste0("combineDelayedMatrix.", j))
    })
    if (length(M) > 1) {
        M <- do.call(cbind, M)
    } else {
        M <- M[[1]]
    }
    M
}

# Internal helper used by combineList()
.combineListMatrixLike <- function(matrix_list, idx_list, nrow, ncol,
                                   fill = NA_real_, hdf5 = FALSE) {
    is_matrix <- vapply(matrix_list, is.matrix, logical(1L))
    is_DelayedMatrix <- vapply(matrix_list, function(x) is(x, "DelayedMatrix"),
                               logical(1L))
    if (all(is_matrix) && !hdf5) {
        z <- matrix(fill, nrow, ncol)
        j0 <- 0
        for (j in seq_along(matrix_list)) {
            i <- idx_list[[j]]
            jj <- j0 + seq_len(ncol(matrix_list[[j]]))
            z[i, jj] <- matrix_list[[j]]
            j0 <- j0 + ncol(matrix_list[[j]])
        }
        if (hdf5) {
            z <- .safeHDF5Array(z, "BSseq.", "z")
        }
    } else if (any(is_DelayedMatrix) || hdf5) {
        # TODO: Could consider using mcmapply(), but (A) is it worth it? and
        #       (B) would need to be careful in case .combineMatrixList() is
        #       itself called from a function using mclapply(). Would also
        #       need to ensure that each object was written to a new .h5
        #       (`file`)
        # TODO: Consolidate into a single .h5 file if running serially
        z <- do.call(cbind, Map(function(matrix, idx) {
            .combineDelayedMatrix(matrix, idx, nrow, fill)
        }, matrix = matrix_list, idx = idx_list))
    } else {
        stop("Cannot combine list of objects with classes: ",
             paste0(vapply(matrix_list, class, character(length(matrix_list))),
                    collapse = ", "))
    }
    z
}

# NOTE: If some assay elements are matrix and some are DelayedMatrix then all
#       matrix objects are first converted to HDF5-backed DelayedMatrix objects
#       before combining (regardless of the value of `hdf5`)
# NOTE: hdf5 argument may be useful even if all inputs are memory-backed
#       because the combined object may be so big that we want to stick it on
#       disk rather than keep it in memory
combineList <- function(x, ..., hdf5 = FALSE) {
    if (class(x) == "BSseq") {
        x <- list(x, ...)
    }
    stopifnot(all(vapply(x, class, character(1L)) == "BSseq"))
    gr <- getBSseq(x[[1]], "gr")
    trans <- getBSseq(x[[1]], "trans")
    sameTrans <- vapply(x[-1], function(xx) {
        identical(trans, getBSseq(xx, "trans"))
    }, logical(1L))
    if (!all(sameTrans)) {
        stop("all elements of '...' in combineList needs to have the same trans")
    }
    sameGr <- vapply(x[-1], function(xx) {
        identical(gr, getBSseq(xx, "gr"))
    }, logical(1L))
    sampleNames <- do.call(c, unname(lapply(x, sampleNames)))
    if (length(sampleNames) > 1 &&
        length(Reduce(intersect, lapply(x, sampleNames))) != 0L) {
        stop("sampleNames of inputs must not overlap")
    }
    if (all(sameGr)) {
        # NOTE: Can shortcut and cbind() the BSseq objects since they all have
        #       the same GRanges. But first need to:
        #       1. Remove smoothing information (coef, se.coef, trans,
        #          parameters) unless all objects have been smoothed
        #       2. Ensure the assays have compatible classes
        #       3. Act on the hdf5 argument
        #
        #       However, calling do.call(cbind, x) defers to
        #       cbind,SummarizedExperiment-method which checks that the
        #       @rowRanges slot in each object is identical **but we've already
        #       done that**. So, to avoid duplicated effort, we simply combine
        #       the remaining slots in the appropriate fashion and construct a
        #       new BSseq object (whilst also avoiding the overhead of the
        #       BSseq() constructor)
        # (1)
        has_been_smoothed <- vapply(x, hasBeenSmoothed, logical(1L))
        if (any(!has_been_smoothed)) {
            x[has_been_smoothed] <- lapply(x[has_been_smoothed], function(xx) {
                # NOTE: Retain all other assays except coef and se.coef
                #       (including any non-standard ones added by the user)
                an <- setdiff(assayNames(xx), c("coef", "se.coef"))
                assays(xx) <- assays(xx)[an]
                xx@trans <- function() NULL
                xx@parameters <- list()
                xx
            })
        }
        # (2) and (3)
        is_matrix <- unlist(lapply(x, function(xx) {
            vapply(assays(xx), is.matrix, logical(1L))
        }), use.names = FALSE)
        is_DelayedMatrix <- unlist(lapply(x, function(xx) {
            vapply(assays(xx), function(a) is(a, "DelayedMatrix"), logical(1L))
        }), use.names = FALSE)
        if ((all(is_matrix) && hdf5) ||
            (!all(is_matrix) && !all(is_DelayedMatrix))) {
            x <- lapply(x, function(xx) {
                assays(xx, withDimnames = TRUE) <- mapply(function(a, an) {
                    # NOTE: Only convert those assays that need converting
                    if (is.matrix(a) && !is(a, "DelayedMatrix")) {
                        dn <- dimnames(a)
                        a <- .safeHDF5Array(a, "BSseq.", an)
                        dimnames(a) <- dn
                    }
                    a
                }, a = assays(xx), an = assayNames(xx))
                xx
            })
        }
        # Adapted from SummarizedExperiment:::.cbind.SummarizedExperiment
        assays <- do.call(cbind, lapply(x, slot, "assays"))
        metadata <- do.call(c, lapply(x, metadata))
        pData <- as(Reduce(combine, lapply(x, function(xx) {
            as.data.frame(pData(xx))
        })), "DataFrame")
        BSseq <- SummarizedExperiment(assays = as(assays, "SimpleList"),
                                      rowRanges = gr,
                                      colData = pData,
                                      metadata = metadata)
        BSseq <- as(BSseq, "BSseq")
    } else {
        # TODO: This will create a very long intermediate GRanges object
        #       and takes a while to run when the objects are long and/or there
        #       are many objects; is there a way to achieve the same outcome
        #       more efficiently (sort is the slowest part)
        gr <- sort(reduce(do.call(c, unname(lapply(x, granges))),
                          min.gapwidth = 0L))
        idx_list <- lapply(x, function(xx) queryHits(findOverlaps(xx, gr)))
        nrow <- length(gr)
        ncol <- length(sampleNames)
        M_list <- lapply(x, function(xx) getBSseq(xx, "M"))
        M <- .combineListMatrixLike(matrix_list = M_list,
                                    idx_list = idx_list,
                                    nrow = nrow,
                                    ncol = ncol,
                                    fill = 0L,
                                    hdf5 = hdf5)
        colnames(M) <- sampleNames
        Cov_list <- lapply(x, function(xx) getBSseq(xx, "Cov"))
        Cov <- .combineListMatrixLike(matrix_list = Cov_list,
                                      idx_list = idx_list,
                                      nrow = nrow,
                                      ncol = ncol,
                                      fill = 0L,
                                      hdf5 = hdf5)
        colnames(Cov) <- sampleNames
        # NOTE: all(sameGr) is FALSE so combined object cannot have any
        #       smoothing information (coef, se.coef, trans, parameters) copied
        #       over from the objects that are combined
        # NOTE: We construct the BSseq object without explicitly calling
        #       BSseq() because it does a lot of argument checking that is not
        #       necessary in this case and otherwise adds considerable
        #       overhead (especially when the objects are large and/or the
        #       assays are HDF5-backed)
        assays <- SimpleList(M = M, Cov = Cov)
        pData <- as(Reduce(combine, lapply(x, function(xx) {
            as.data.frame(pData(xx))
        })), "DataFrame")
        metadata <- do.call(c, lapply(x, metadata))
        BSseq <- SummarizedExperiment(assays = assays, rowRanges = gr,
                                      colData = pData, metadata = metadata)
        BSseq <- as(BSseq, "BSseq")
    }
    BSseq
}

# NOTE: No easy way to include a `hdf5` argument to combine() due to the
#       generic's definition in BiocGenerics. Instead, the returned value is a
#       HDF5-backed DelayedArray if any of the assays in x, y, or ... are
#       themselves DelayedArray objects (even if not necessarily HDF5-backed)
setMethod("combine", signature(x = "BSseq", y = "BSseq"), function(x, y, ...) {
    args <- list(x, y, ...)
    is_DelayedArray <- vapply(args, function(xx) {
        any(vapply(assays(xx), function(assay) {
            is(assay, "DelayedArray")
        }, logical(1L)))
    }, logical(1L))
    if (any(is_DelayedArray)) {
        hdf5 <- TRUE
    } else {
        hdf5 <- FALSE
    }
    combineList(x, y, ..., hdf5 = hdf5)
})
