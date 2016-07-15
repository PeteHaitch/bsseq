# TODO: Do really need both combine() and combineList(); combineList() is
#       demonstrably faster when > 2 BSseq objects and same speed when
#       2 BSseq objects. Maintaining both introduces redundancy and overhead.
#       If combine is required (after all, it is a BiocGeneric), then can we
#       define it in terms of combineList()?

# NOTE: No easy way to include a `hdf5` argument to combine() due to the
#       generic's definition in BiocGenerics. Instead, the returned value is a
#       HDF5-backed DelayedMatrix if any of the assays in x, y, or ... are
#       themselves DelayedMatrix objects.
# NOTE: combine()-ing HDF5-backed BSseq objects can result in assay elements
#       that are DelayedArray objects with multiple HDF5-backed seed, which
#       can be slow to process in subsequent computation. It may be a good idea
#       to realise the assays to disk as new object(s) in the HDF5 file
# TODO: Benchmark the above claim (which is currently anecdotal)
setMethod("combine", signature(x = "BSseq", y = "BSseq"), function(x, y, ...) {
    ## All of this assumes that we are place x and y "next" to each other,
    ##  ie. we are not combining the same set of samples sequenced at different times
    if (class(x) != class(y)) {
        stop(paste("objects must be the same class, but are ",
                   class(x), ", ", class(y), sep = ""))
    }
    # TODO: Should this be checking @parameters
    if (hasBeenSmoothed(x) && hasBeenSmoothed(y) && !all.equal(x@trans, y@trans)) {
        stop("'x' and 'y' need to be smoothed on the same scale")
    }
    pData <- combine(as(pData(x), "data.frame"), as(pData(y), "data.frame"))
    if (length(intersect(sampleNames(x), sampleNames(y))) != 0L) {
        stop("sampleNames of 'x' and 'y' must not overlap")
    }
    if (identical(granges(x), granges(y))) {
        # TODO: Could use the same cbind() trick as is used by combineList()
        #       but need to first settle on whether it should be allowed for a
        #       smoothed to be combined with a non-smoothed BSseq object
        # TODO: This assumes that x and y are both matrix-backed or
        #       DelayedMatrix-backed and will break if not true
        gr <- granges(x)
        M_x <- getBSseq(x, "M")
        M_y <- getBSseq(y, "M")
        M <- cbind(M_x, M_y)
        Cov_x <- getBSseq(x, "Cov")
        Cov_y <- getBSseq(y, "Cov")
        Cov <- cbind(Cov_x, Cov_y)
        if (!hasBeenSmoothed(x) || !hasBeenSmoothed(y)) {
            coef <- NULL
            se.coef <- NULL
            trans <- NULL
        } else {
            coef_x <- getBSseq(x, "coef")
            coef_y <- getBSseq(y, "coef")
            coef <- cbind(coef_x, coef_y)
            se.coef_x <- getBSseq(x, "se.coef")
            se.coef_y <- getBSseq(y, "se.coef")
            # NOTE: It's possible that only one object has non-NULL se.coef
            #       even if both objects have been smoothed. In this case,
            #       we create a matrix filled with NA so that the combined
            #       se.coef has the appropriate dimension
            # NOTE: These are always realised in memory and all elements
            #       set to NA_real_ (rather than 0)
            if (is.null(se.coef_x)) {
                se.coef_x <- matrix(NA_real_)
                se.coef_x <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x))
            }
            if (is.null(se.coef_y)) {
                se.coef_y <- matrix(NA_real_)
                se.coef_y <- matrix(NA_real_, nrow = nrow(y), ncol = ncol(y))
            }
            se.coef <- cbind(se.coef_x, se.coef_y)
            trans <- getBSseq(x, "trans")
        }
    } else {
        gr <- reduce(c(granges(x), granges(y)), min.gapwidth = 0L)
        ov_x <- findOverlaps(gr, granges(x))
        ov_y <- findOverlaps(gr, granges(y))
        sampleNames <- c(sampleNames(x), sampleNames(y))
        M_x <- getBSseq(x, "M")[subjectHits(ov_x), , drop = FALSE]
        M_y <- getBSseq(y, "M")[subjectHits(ov_y), , drop = FALSE]
        Cov_x <- getBSseq(x, "Cov")[subjectHits(ov_x), , drop = FALSE]
        Cov_y <- getBSseq(y, "Cov")[subjectHits(ov_y), , drop = FALSE]
        nrow <- length(gr)
        ncol <- length(sampleNames)
        # M <- .combineMatrixLike(x = M_x, y = M_y, idx_x = queryHits(ov_x),
        #                         idx_y = queryHits(ov_y), nrow = nrow,
        #                         ncol = ncol, fill = 0L)
        M  <- .combineListMatrixLike(matrix_list = list(M_x, M_y),
                                     idx_list = list(queryHits(ov_x),
                                                     queryHits(ov_y)),
                                     nrow = nrow, ncol = ncol, fill = 0L)
        # Cov <- .combineMatrixLike(x = Cov_x, y = Cov_y, idx_x = queryHits(ov_x),
        #                           idx_y = queryHits(ov_y), nrow = nrow,
        #                           ncol = ncol, fill = 0L)
        Cov  <- .combineListMatrixLike(matrix_list = list(Cov_x, Cov_y),
                                     idx_list = list(queryHits(ov_x),
                                                     queryHits(ov_y)),
                                     nrow = nrow, ncol = ncol, fill = 0L)
        colnames(M) <- sampleNames
        colnames(Cov) <- sampleNames
        if (!hasBeenSmoothed(x) || !hasBeenSmoothed(y)) {
            # TODO: This means if one object has been smoothed and the other
            #       hasn't then the smoothing is removed with no warning or
            #       message; discuss with Kasper (note that combineList()
            #       will only combine coef, se.coef, trans if in all objects
            #       have been smoothed and have the same GRanges; this seems
            #       reasonable to require/enforce)
            coef <- NULL
            se.coef <- NULL
            trans <- NULL
        } else {
            coef_x <- getBSseq(x, "coef")[subjectHits(ov_x), , drop = FALSE]
            coef_y <- getBSseq(y, "coef")[subjectHits(ov_y), , drop = FALSE]
            # TODO: fill = NA_real_ rather than 0L because coef is numeric
            #       and missing; discuss with Kasper but I think this is the
            #       correct choice
            # coef <- .combineMatrixLike(x = coef_x, y = coef_y,
            #                            idx_x = queryHits(ov_x),
            #                            idx_y = queryHits(ov_y),
            #                            nrow = nrow, ncol = ncol,
            #                            fill = NA_real_)
            coef  <- .combineListMatrixLike(matrix_list = list(coef_x, coef_y),
                                         idx_list = list(queryHits(ov_x),
                                                         queryHits(ov_y)),
                                         nrow = nrow, ncol = ncol,
                                         fill = NA_real_)
            colnames(coef) <- sampleNames
            if (is.null(getBSseq(x, "se.coef")) &&
                is.null(getBSseq(x, "se.coef"))) {
                se.coef <- NULL
            } else {
                se.coef_x <- getBSseq(x, "se.coef")[subjectHits(ov_x), , drop = FALSE]
                se.coef_y <- getBSseq(y, "se.coef")[subjectHits(ov_y), , drop = FALSE]
                # NOTE: These are always realised in memory and all elements
                #       set to NA_real_ (rather than 0)
                if (is.null(se.coef_x)) {
                    se.coef_x <- matrix(NA_real_, nrow = queryHits(ov_x),
                                        ncol = ncol(x))
                }
                if (is.null(se.coef_y)) {
                    se.coef_y <- matrix(NA_real_, nrow = queryHits(ov_y),
                                        ncol = ncol(y))
                }
                # TODO: fill = NA_real_ rather than 0L because se.coef is
                #       numeric and missing; discuss with Kasper but I think
                #       this is the correct choice
                # se.coef <- .combineMatrixLike(x = coef_x, y = coef_y,
                #                               idx_x = queryHits(ov_x),
                #                               idx_y = queryHits(ov_y),
                #                               nrow = nr, ncol = nc,
                #                               fill = NA_real_)
                se.coef  <- .combineListMatrixLike(matrix_list = list(se.coef_x,
                                                                      se.coef_y),
                                             idx_list = list(queryHits(ov_x),
                                                             queryHits(ov_y)),
                                             nrow = nrow, ncol = ncol,
                                             fill = NA_real_)
                colnames(se.coef) <- sampleNames
            }
            trans <- getBSseq(x, "trans")
        }
    }
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          pData = pData, trans = trans, rmZeroCov = FALSE)
})

# TODO: Check columns of output are correctly ordered
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
    stopifnot(all(sapply(x, class) == "BSseq"))
    gr <- getBSseq(x[[1]], "gr")
    trans <- getBSseq(x[[1]], "trans")
    sameTrans <- sapply(x[-1], function(xx) {
        identical(trans, getBSseq(xx, "trans"))
    })
    if (!all(sameTrans)) {
        stop("all elements of '...' in combineList needs to have the same trans")
    }
    sameGr <- sapply(x[-1], function(xx) {
        identical(gr, getBSseq(xx, "gr"))
    })
    if (length(sampleNames) > 1 &&
        length(Reduce(intersect, lapply(x, sampleNames))) != 0L) {
        stop("sampleNames of inputs must not overlap")
    }
    if (all(sameGr)) {
        # NOTE: Can shortcut and cbind the BSseq objects since they all have
        #       the same GRanges. But first need to ensure the assays have
        #       compatible classes and also act on the hdf5 argument.
        is_matrix <- unlist(lapply(x, function(xx) {
            vapply(assays(xx), is.matrix, logical(1L))
        }), use.names = FALSE)
        is_DelayedMatrix <- unlist(lapply(x, function(xx) {
            vapply(assays(xx), function(a) is(a, "DelayedMatrix"), logical(1L))
        }), use.names = FALSE)
        if ((all(is_matrix) && hdf5) ||
            (!all(is_matrix) && !all(is_DelayedMatrix))) {
            x <- lapply(x, function(xx) {
                assays(xx) <- endoapply(assays(xx), function(a) {
                    # NOTE: Only convert those assays that need converting
                    if (is.matrix(a) && !is(a, "DelayedMatrix")) {
                        a <- HDF5Array(a)
                    }
                    a
                })
                xx
            })
        }
        return(do.call(cbind, x))
    } else {
        gr <- sort(reduce(do.call(c, unname(lapply(x, granges))),
                          min.gapwidth = 0L))
        idx_list <- lapply(x, function(xx) queryHits(findOverlaps(xx, gr)))
        sampleNames <- do.call(c, unname(lapply(x, sampleNames)))
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
    }
    if (any(!sapply(x, hasBeenSmoothed)) || !(all(sameGr))) {
        coef <- NULL
        se.coef <- NULL
        trans <- NULL
    } else {
        coef <- do.call(cbind, lapply(x, getBSseq, "coef"))
        se.coef <- do.call(cbind, lapply(x, getBSseq, "se.coef"))
        if (hdf5) {
            coef <- HDF5Array(coef)
            if (!is.null(se.coef)) {
                se.coef <- HDF5Array(coef)
            }
        }
    }
    pData <- as(Reduce(combine, lapply(x, function(xx) {
        as.data.frame(pData(xx))
    })), "DataFrame")
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          pData = pData, trans = trans, rmZeroCov = FALSE)
}
