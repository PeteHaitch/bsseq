# TODO: Update all files that create .h5 files to use
#       HDF5Array(writeHDF5Dataset(x, file = .newBSseqHDF5Filename(), name =
#                 "callingFunction") syntax
# TODO: Do really need both combine() and combineList(); combineList() is
#       demonstrably faster when > 2 BSseq objects and same speed when
#       2 BSseq objects. Maintaining both introduces redundancy and overhead.
#       If combine is required (after all, it is a BiocGeneric), then can we
#       define it in terms of combineList()?

# Internal helper used by .combineMatrixLike and .combineListMatrixList
# While designed for DelayedMatrix objects, it will also work for matrix
# objects (although it is not the most efficient way to do so in this case
# since it (A) cbinds() vectors to form the final matrix-like object, and
# (B) returns a HDF5-backed DelayedMatrix object rather than a matrix object)
# NOTE: i is used to ensure names are unique when lapply()-ing this function
#       on a list of DelayedMatrix objects (see .combineListMatrixLike for an
#       example)
.combineDelayedMatrix <- function(m, idx, nrow, fill, file, i = 1) {
    # TODO: Could consider using mclapply(), but (A) is it worth it? and
    #       (B) would need to be careful in case .combineMatrixList() is
    #       itsself called from a function using mclapply(). Would also need
    #       to ensure that each object was written to a new .h5 (`file`)
    M <- lapply(seq_len(ncol(m)), function(j) {
        mm <- matrix(fill, nrow = nrow, ncol = 1L)
        mm[idx, 1L] <- as.array(m[, j, drop = FALSE])
        HDF5Array(writeHDF5Dataset(mm, file = file,
                                   name = paste0("combineDelayedMatrix.", i,
                                                 ".", j)))
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
            hdf5_file <- .newBSseqHDF5Filename()
            z <- HDF5Array(writeHDF5Dataset(z, file = hdf5_file,
                                            name = "combineListMatrixLike"))
        }
    } else if (any(is_DelayedMatrix) || hdf5) {
        # TODO: Could consider using mcmapply(), but (A) is it worth it? and
        #       (B) would need to be careful in case .combineMatrixList() is
        #       itself called from a function using mclapply(). Would also
        #       need to ensure that each object was written to a new .h5
        #       (`file`)
        hdf5_file <- .newBSseqHDF5Filename()
        z <- do.call(cbind, lapply(seq_along(matrix_list), function(i) {
            .combineDelayedMatrix(matrix_list[[i]], idx_list[[i]], nrow,
                                  fill, hdf5_file, i)
        }))
    } else {
        stop("Cannot combine list of objects with classes: ",
             paste0(vapply(matrix_list, class, character(length(matrix_list))),
                    collapse = ", "))
    }
    z
}

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
        # Adapted from SummarizedExperiment:::.cbind.SummarizedExperiment
        assays <- do.call(cbind, lapply(x, slot, "assays"))
        metadata <- do.call(c, lapply(x, metadata))
        pData <- as(Reduce(combine, lapply(x, function(xx) {
            as.data.frame(pData(xx))
        })), "DataFrame")
        BSseq <- SummarizedExperiment(assays = assays, rowRanges = gr,
                                      colData = pData, metadata = metadata)
        BSseq <- as(BSseq, "BSseq")
    } else {
        # TODO: This will create a very long intermediate GRanges object
        #       and takes a while to run when the objects are long and/or there
        #       are many objects; is there a way to achieve the same outcome
        #       more efficiently (see commented code for a candidate)
        gr <- sort(reduce(do.call(c, unname(lapply(x, granges))),
                          min.gapwidth = 0L))
        # gr <- sort(Reduce(function(x, y) reduce(c(x, y), min.gapwidth = 0L),
        #                   unname(lapply(x[-1], granges)), init = granges(x[[1]])))
        # TODO: Might be more efficient to compute elements of idx_list on the
        #       fly (e.g., inside call to .combineListMatrixList()) since we
        #       never simultaneously need all elements of idx_list
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
