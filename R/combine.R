# NOTE: If any of the objects have DelayedArray assay elements, then the
#       corresponding combined assay element is a HDF5Array
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
            # TODO: Check this conditional
            trans <- x@trans
            coef <- matrix(0, nrow = length(gr), ncol = length(sampleNames))
            coef_hdf5 <- is(getBSseq(x, "coef"), "DelayedArray") ||
                is(getBSseq(y, "coef"), "DelayedArray")
            colnames(coef) <- rownames(pData)
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
            if(is.null(getBSseq(x, "se.coef")) && is.null(getBSseq(x, "se.coef")))
                se.coef <- NULL
            else {
                se.coef <- matrix(0, nrow = length(gr),
                                  ncol = length(sampleNames))
                se.coef_hdf5 <- is(getBSseq(x, "se.coef"), "DelayedArray") ||
                    is(getBSseq(y, "se.coef"), "DelayedArray")
                colnames(se.coef) <- sampleNames(pData)
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
combineList <- function(x, ...) {
    if(class(x) == "BSseq")
        x <- list(x, ...)
    stopifnot(all(sapply(x, class) == "BSseq"))
    gr <- getBSseq(x[[1]], "gr")
    trans <- getBSseq(x[[1]], "trans")
    sameTrans <- sapply(x[-1], function(xx) {
        identical(trans, getBSseq(xx, "trans"))
    })
    if(!all(sameTrans))
        stop("all elements of '...' in combineList needs to have the same trans")
    sameGr <- sapply(x[-1], function(xx) {
        identical(gr, getBSseq(xx, "gr"))
    })
    if(all(sameGr)) {
        # TODO: Errors on list of DelayedMatrix objects (coming from
        #       read.bismark). Remove unname() once fixed.
        M <- do.call(cbind, unname(lapply(x, function(xx) getBSseq(xx, "M"))))
        # TODO: Errors on list of DelayedMatrix objects (coming from
        #       read.bismark). Remove unname9) once fixed
        Cov <- do.call(cbind, unname(lapply(x, function(xx) getBSseq(xx, "Cov"))))
    } else {
        gr <- sort(reduce(do.call(c, unname(lapply(x, granges))), min.gapwidth = 0L))
        sampleNames <- do.call(c, unname(lapply(x, sampleNames)))
        M <- matrix(0, ncol = length(sampleNames), nrow = length(gr))
        colnames(M) <- sampleNames
        Cov <- M
        idxes <- split(seq(along = sampleNames), rep(seq(along = x), times = sapply(x, ncol)))
        for(ii in seq(along = idxes)) {
            ov <- findOverlaps(gr, granges(x[[ii]]))
            idx <- idxes[[ii]]
            M[queryHits(ov), idx] <- getBSseq(x[[ii]], "M")[subjectHits(ov),]
            Cov[queryHits(ov), idx] <- getBSseq(x[[ii]], "Cov")[subjectHits(ov),]
        }
    }
    if(any(!sapply(x, hasBeenSmoothed)) || !(all(sameGr))) {
        coef <- NULL
        se.coef <- NULL
        trans <- NULL
    } else {
        coef <- do.call(cbind, lapply(x, function(xx) getBSseq(xx, "coef")))
        se.coef <- do.call(cbind, lapply(x, function(xx) getBSseq(xx, "se.coef")))
    }
    pData <- as(Reduce(combine, lapply(x, function(xx) as.data.frame(pData(xx)))), "DataFrame")
    BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = se.coef,
          pData = pData, trans = trans, rmZeroCov = FALSE)
}

