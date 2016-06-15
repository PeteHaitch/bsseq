# NOTE: If any of the objects have an assay element that is a DelayedArray,
#       then the corresponding combined assay element is a HDF5Array.
#       Furthermore, if combine() is called with `...`, then each recursive
#       call will generate a new HDF5 file; this may be slow. An alternative
#       implementation would be to return a DelayedArray object (which may
#       comprise some combination of matrix, HDF5Matrix, DelayedMatrix objects)
#       but we have no way of automatically converting this to a HDF5Array at
#       the top-level function call.
# NOTE: No easy way to include a `hdf5` argument to combine() due to the
#       generic's definition in BiocGenerics
setMethod("combine", signature(x = "BSseq", y = "BSseq"), function(x, y, ...) {
    ## All of this assumes that we are place x and y "next" to each other,
    ##  ie. we are not combining the same set of samples sequenced at different times
    if (class(x) != class(y))
        stop(paste("objects must be the same class, but are ",
                   class(x), ", ", class(y), sep=""))
    if(hasBeenSmoothed(x) && hasBeenSmoothed(y) && !all.equal(x@trans, y@trans))
        stop("'x' and 'y' need to be smoothed on the same scale")
    pData <- combine(as(pData(x), "data.frame"), as(pData(y), "data.frame"))
    if(identical(granges(x), granges(y))) {
        gr <- granges(x)
        M_x <- getBSseq(x, "M")
        M_y <- getBSseq(y, "M")
        # NOTE: Have to realise DelayedArray objects because cannot do
        #       sub-assignment of DelayedArray into the matrix objects, `M` and
        #       `Cov`
        if (is(M_x, "DelayedArray") && is.matrix(M_y)) {
            M_y <- HDF5Array(M_y)
        }
        if (is.matrix(M_x) && is(M_y, "DelayedArray")) {
            M_x <- HDF5Array(M_x)
        }
        M <- cbind(M_x, M_y)
        Cov_x <- getBSseq(x, "Cov")
        Cov_y <- getBSseq(y, "Cov")
        if (is(Cov_x, "DelayedArray") && is.matrix(Cov_y)) {
            Cov_y <- HDF5Array(Cov_y)
        }
        if (is.matrix(Cov_x) && is(Cov_y, "DelayedArray")) {
            Cov_x <- HDF5Array(Cov_x)
        }
        Cov <- cbind(Cov_x, Cov_y)
        if(!hasBeenSmoothed(x) || !hasBeenSmoothed(y)) {
            coef <- NULL
            se.coef <- NULL
            trans <- NULL
        } else {
            # NOTE: Have to realise DelayedArray objects because cannot do
            #       sub-assignment of DelayedArray into the matrix object,
            #       `coef`
            coef_x <- getBSseq(x, "coef")
            coef_y <- getBSseq(y, "coef")
            if (is(coef_x, "DelayedArray") && is.matrix(coef_y)) {
                coef_y <- HDF5Array(coef_y)
            }
            if (is.matrix(coef_x) && is(coef_y, "DelayedArray")) {
                coef_x <- HDF5Array(coef_x)
            }
            coef <- cbind(coef_x, coef_y)

            se.coef_x <- getBSseq(x, "se.coef")
            se.coef_y <- getBSseq(y, "se.coef")
            # NOTE: Have to realise DelayedArray objects because cannot do
            #       sub-assignment of DelayedArray into the matrix object,
            #       `se.coef`
            if (is(se.coef_x, "DelayedArray") && is.matrix(se.coef_y)) {
                se.coef_y <- HDF5Array(se.coef_y)
            }
            if (is.matrix(se.coef_x) && is(se.coef_y, "DelayedArray")) {
                se.coef_x <- HDF5Array(se.coef_x)
            }
            se.coef <- cbind(se.coef_x, se.coef_y)
            trans <- getBSseq(x, "trans")
        }
    } else {
        gr <- reduce(c(granges(x), granges(y)), min.gapwidth = 0L)
        mm.x <- as.matrix(findOverlaps(gr, granges(x)))
        mm.y <- as.matrix(findOverlaps(gr, granges(y)))
        sampleNames <- c(sampleNames(x), sampleNames(y))
        ## FIXME: there is no check that the two sampleNames are disjoint.
        M <- Cov <- matrix(0, nrow = length(gr), ncol = length(sampleNames))
        M_hdf5 <- is(getBSseq(x, "M"), "DelayedArray") ||
            is(getBSseq(y, "M"), "DelayedArray")
        Cov_hdf5 <- is(getBSseq(x, "Cov"), "DelayedArray") ||
            is(getBSseq(y, "Cov"), "DelayedArray")
        colnames(M) <- colnames(Cov) <- sampleNames
        # NOTE: Have to realise DelayedArray objects because cannot do
        #       sub-assignment of DelayedArray into the matrix objects, `M` and
        #       `Cov`
        M_x <- getBSseq(x, "M")[mm.x[,2],]
        if (is(M_x, "DelayedArray")) {
            M_x <- as.array(M_x)
        }
        M[mm.x[,1], 1:ncol(x)] <- M_x
        M_y <- getBSseq(y, "M")[mm.y[,2],]
        if (is(M_y, "DelayedArray")) {
            M_y <- as.array(M_y)
        }
        M[mm.y[,1], ncol(x) + 1:ncol(y)] <- M_y
        if (M_hdf5) {
            M <- HDF5Array(M)
        }
        Cov_x <- getBSseq(x, "Cov")[mm.x[,2],]
        if (is(Cov_x, "DelayedArray")) {
            Cov_x <- as.array(Cov_x)
        }
        Cov[mm.x[,1], 1:ncol(x)] <- Cov_x
        Cov_y <- getBSseq(y, "Cov")[mm.y[,2],]
        if (is(Cov_y, "DelayedArray")) {
            Cov_y <- as.array(Cov_y)
        }
        Cov[mm.y[,1], ncol(x) + 1:ncol(y)] <- Cov_y
        if (Cov_hdf5) {
            Cov <- HDF5Array(Cov)
        }
        if(!hasBeenSmoothed(x) || !hasBeenSmoothed(y)) {
            coef <- NULL
            se.coef <- NULL
            trans <- NULL
        } else {
            trans <- x@trans
            coef <- matrix(0, nrow = length(gr), ncol = length(sampleNames))
            coef_hdf5 <- is(getBSseq(x, "coef"), "DelayedArray") ||
                is(getBSseq(y, "coef"), "DelayedArray")
            colnames(coef) <- rownames(pData)
            # NOTE: Have to realise DelayedArray objects because cannot do
            #       sub-assignment of DelayedArray into the matrix object,
            #       `coef`
            if(hasBeenSmoothed(x)) {
                coef_x <- getBSseq(x, "coef")[mm.x[,2],]
                if (is(coef_x, "DelayedArray")) {
                    coef_x <- as.array(M_x)
                }
                coef[mm.x[,1], 1:ncol(x)] <- coef_x
            }
            if(hasBeenSmoothed(y)) {
                coef_y <- getBSseq(y, "coef")[mm.y[,2],]
                if (is(coef_y, "DelayedArray")) {
                    coef_y <- as.array(M_y)
                }
                coef[mm.y[,1], ncol(x) + 1:ncol(x)] <- coef_y
            }
            if (coef_hdf5) {
                coef <- HDF5Array(coef)
            }
            if(is.null(getBSseq(x, "se.coef")) &&
               is.null(getBSseq(x, "se.coef"))) {
                se.coef <- NULL
            }
            else {
                se.coef <- matrix(0, nrow = length(gr),
                                  ncol = length(sampleNames))
                se.coef_hdf5 <- is(getBSseq(x, "se.coef"), "DelayedArray") ||
                    is(getBSseq(y, "se.coef"), "DelayedArray")
                colnames(se.coef) <- sampleNames(pData)
                # NOTE: Have to realise DelayedArray objects because cannot do
                #       sub-assignment of DelayedArray into the matrix object,
                #       `se.coef`
                if(!is.null(getBSseq(x, "se.coef"))) {
                    se.coef_x <- getBSseq(x, "se.coef")[mm.x[,2],]
                    if (is(se.coef_x, "DelayedArray")) {
                        se.coef_x <- as.array(M_x)
                    }
                    se.coef[mm.x[,1], 1:ncol(x)] <- se.coef_x
                }
                if(!is.null(getBSseq(y, "se.coef"))) {
                    se.coef_y <- getBSseq(y, "se.coef")[mm.y[,2],]
                    if (is(se.coef_y, "DelayedArray")) {
                        se.coef_y <- as.array(M_y)
                    }
                    se.coef[mm.y[,1], ncol(x) + 1:ncol(x)] <- se.coef_y
                }
                if (se.coef_hdf5) {
                    se.coef <- HDF5Array(se.coef)
                }
            }
        }
    }
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          pData = pData, trans = trans, rmZeroCov = FALSE)
})

# UP TO HERE
# NOTE: `hdf5` determines whether any DelayedArray objects created by
#       combineList() are written to disk as **new** HDF5 files. This is
#       regardless of whether any of the input objects were HDF5-backed.
# TODO: If using HDF5Array, should we write M and Cov to disk as a new
#       HDF5 file (i.e. call HDF5Array(M)) or keep them as a
#       DelayedMatrix object?
combineList <- function(x, ..., hdf5 = FALSE) {
    if(class(x) == "BSseq")
        x <- list(x, ...)
    stopifnot(all(sapply(x, class) == "BSseq"))
    gr <- getBSseq(x[[1]], "gr")
    trans <- getBSseq(x[[1]], "trans")
    sameTrans <- sapply(x[-1], function(xx) {
        identical(trans, getBSseq(xx, "trans"))
    })
    M_hdf5 <- any(vapply(x, function(xx) {
        is(getBSseq(xx, "M"), "DelayedArray")
    }, logical(1L)))
    Cov_hdf5 <- any(vapply(x, function(xx) {
        is(getBSseq(xx, "Cov"), "DelayedArray")
    }, logical(1L)))
    if(!all(sameTrans))
        stop("all elements of '...' in combineList needs to have the same trans")
    sameGr <- sapply(x[-1], function(xx) {
        identical(gr, getBSseq(xx, "gr"))
    })
    if(all(sameGr)) {
        # NOTE: Can't cbind matrix and DelayedMatrix so need to coerce matrix
        #       to a DelayedArray
        M <- do.call(cbind, lapply(x, function(xx) {
            M <- getBSseq(xx, "M")
            if (M_hdf5 && is.matrix(M)) {
                M <- DelayedArray(M)
            }
            M
        }))
        if (hdf5) {
            M <- HDF5Array(M)
        }
        Cov <- do.call(cbind, lapply(x, function(xx) {
            Cov <- getBSseq(xx, "Cov")
            if (Cov_hdf5 && is.matrix(Cov)) {
                Cov <- DelayedArray(Cov)
            }
            Cov
        }))
        if (hdf5) {
            M <- HDF5Array(M)
        }
    } else {
        gr <- sort(reduce(do.call(c, unname(lapply(x, granges))),
                          min.gapwidth = 0L))
        sampleNames <- do.call(c, unname(lapply(x, sampleNames)))
        M <- matrix(0, ncol = length(sampleNames), nrow = length(gr))
        colnames(M) <- sampleNames
        Cov <- M
        idxes <- split(seq(along = sampleNames),
                       rep(seq(along = x), times = sapply(x, ncol)))
        for(ii in seq(along = idxes)) {
            ov <- findOverlaps(gr, granges(x[[ii]]))
            idx <- idxes[[ii]]
            # NOTE: Have to realise DelayedArray objects because cannot do
            #       sub-assignment of DelayedArray into the matrix objects,
            #       `M` and `Cov`
            M_xx <- getBSseq(x[[ii]], "M")[subjectHits(ov),]
            if (is(M_xx, "DelayedArray")) {
                M_xx <- as.array(M_xx)
            }
            M[queryHits(ov), idx] <- M_xx
            Cov_xx <- getBSseq(x[[ii]], "Cov")[subjectHits(ov),]
            if (is(Cov_xx, "DelayedArray")) {
                Cov_xx <- as.array(Cov_xx)
            }
            Cov[queryHits(ov), idx] <- Cov_xx
        }
    }
    if (hdf5) {
        M <- HDF5Array(M)
        Cov <- HDF5Array(Cov)
    }
    if(any(!sapply(x, hasBeenSmoothed)) || !(all(sameGr))) {
        coef <- NULL
        se.coef <- NULL
        trans <- NULL
    } else {
        coef <- do.call(cbind, lapply(x, function(xx) getBSseq(xx, "coef")))
        se.coef <- do.call(cbind, lapply(x, function(xx) {
            getBSseq(xx, "se.coef")}
        ))
        if (hdf5) {
            coef <- HDF5Array(coef)
            se.coef <- HDF5Array(coef)
        }
    }
    pData <- as(Reduce(combine, lapply(x, function(xx) {
        as.data.frame(pData(xx))
    })), "DataFrame")
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          pData = pData, trans = trans, rmZeroCov = FALSE)
}
