collapseBSseq <- function(BSseq, columns) {
    ## columns is a vector of new names, names(columns) is sampleNames or empty
    stopifnot(is.character(columns))
    if(is.null(names(columns)) && length(columns) != ncol(BSseq))
        stop("if `columns' does not have names, it needs to be of the same length as `BSseq` has columns (samples)")
    if(!is.null(names(columns)) && !all(names(columns) %in% sampleNames(BSseq)))
        stop("if `columns` has names, they need to be sampleNames(BSseq)")
    if(is.null(names(columns)))
        columns.idx <- 1:ncol(BSseq)
    else
        columns.idx <- match(names(columns), sampleNames(BSseq))
    sp <- split(columns.idx, columns)
    # TODO: .collapseDelayedMatrix() always return numeric; it may be
    #       worth coercing M and Cov to integer DelayedMatrix objects,
    #       which would halve storage requirements and impose some more
    #       structure on the BSseq class (via new validity method
    #       checks)
    # NOTE: Tries to be smart about how collapsed DelayedMatrix should
    #       be realized
    M <- getBSseq(BSseq, "M")
    if (.isHDF5ArrayBacked(M)) {
        M_BACKEND <- "HDF5Array"
    } else {
        M_BACKEND <- NULL
    }
    M <- .collapseDelayedMatrix(x = M,
                                sp = sp,
                                MARGIN = 1,
                                BACKEND = M_BACKEND)
    Cov <- getBSseq(BSseq, "Cov")
    if (.isHDF5ArrayBacked(Cov)) {
        Cov_BACKEND <- "HDF5Array"
    } else {
        Cov_BACKEND <- NULL
    }
    Cov <- .collapseDelayedMatrix(x = Cov,
                                  sp = sp,
                                  MARGIN = 1,
                                  BACKEND = Cov_BACKEND)
    BSseq(gr = getBSseq(BSseq, "gr"), M = M, Cov = Cov, sampleNames = names(sp))
}

chrSelectBSseq <- function(BSseq, seqnames = NULL, order = FALSE) {
    seqlevels(BSseq, pruning.mode = "coarse") <- seqnames
    if(order)
        BSseq <- orderBSseq(BSseq, seqOrder = seqnames)
    BSseq
}


orderBSseq <- function(BSseq, seqOrder = NULL) {
    if(!is.null(seqOrder))
        seqlevels(BSseq, pruning.mode = "coarse") <- seqOrder
    BSseq[order(granges(BSseq))]
}

# TODO: Add BACKEND argument?
getMeth <- function(BSseq, regions = NULL, type = c("smooth", "raw"),
                    what = c("perBase", "perRegion"), confint = FALSE,
                    alpha = 0.95, mc.cores = 1) {
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    if (type == "smooth" & !hasBeenSmoothed(BSseq)) {
        stop("'type = smooth' requires the object to have been smoothed.")
    }
    what <- match.arg(what)

    # NOTE: Computing meth is a delayed operation, so it is essentially free to
    #       do this regardless of the size of the BSseq object
    trans <- getBSseq(BSseq, type = "trans")
    if (type == "smooth") {
        meth <- trans(getBSseq(BSseq, type = "coef"))
    } else if (type == "raw") {
        meth <- getBSseq(BSseq, type = "M") / getBSseq(BSseq, type = "Cov")
    }

    if (is.null(regions)) {
        if (!confint) {
            return(meth)
        }
        z <- abs(qnorm((1 - alpha) / 2, mean = 0, sd = 1))
        if (type == "smooth") {
            upper <- trans(getBSseq(BSseq, type = "coef") +
                               z * getBSseq(BSseq, type = "se.coef"))
            lower <- trans(getBSseq(BSseq, type = "coef") -
                               z * getBSseq(BSseq, type = "se.coef"))
            return(list(meth = meth, lower = lower, upper = upper))
        } else if (type == "raw")
            # Function to compute confidence interval for a binomial
            # probability (p)
            p.conf <- function(p, n, alpha) {
                z <- abs(qnorm((1 - alpha) / 2, mean = 0, sd = 1))
                upper <- (p + z ^ 2 / (2 * n) + z *
                              sqrt((p * (1 - p) + z ^ 2 / (4 * n)) / n)) /
                    (1 + z ^ 2 / n)
                lower <- (p + z ^ 2 / (2 * n) - z *
                              sqrt((p * (1 - p) + z ^ 2 / (4 * n)) / n)) /
                    (1 + z ^ 2 / n)
                return(list(meth = p, lower = lower, upper = upper))
            }
        return(p.conf(p = meth,
                      n = getBSseq(BSseq, type = "Cov"),
                      alpha = alpha))
    }

    # At this point, regions have been specified
    if (class(regions) == "data.frame") {
        regions <- data.frame2GRanges(regions)
    }
    stopifnot(is(regions, "GenomicRanges"))
    if (confint) {
        stop("'confint = TRUE' is not supported by 'getMeth' when regions ",
             "is given")
    }
    list_of_hits <- as(findOverlaps(regions, BSseq), "List")
    if (what == "perBase") {
        out <- mclapply(list_of_hits, function(i) {
            meth[i, ]
        }, mc.cores = mc.cores)
        return(out)
    } else if (what == "perRegion") {
        out <- mclapply(list_of_hits, function(i) {
            colMeans(meth[i, ], na.rm = TRUE)
        }, mc.cores = mc.cores)

        out <- do.call(rbind, out)
        colnames(out) <- sampleNames(BSseq)
        .DelayedMatrix(out)
    }
}

# TODO: Add BACKEND argument?
getCoverage <- function(BSseq, regions = NULL, type = c("Cov", "M"),
                        what = c("perBase", "perRegionAverage",
                                 "perRegionTotal"),
                        mc.cores = 1) {
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    what <- match.arg(what)

    out <- getBSseq(BSseq, type = type)
    if (is.null(regions)) {
        if (what == "perBase") {
            return(out)
        }
        if (what == "perRegionTotal") {
            return(colSums(out))
        }
        if (what == "perRegionAverage") {
            return(colMeans(out))
        }
    }

    # At this point, regions have been specified
    if (class(regions) == "data.frame") {
        regions <- data.frame2GRanges(regions)
    }
    stopifnot(is(regions, "GenomicRanges"))

    list_of_hits <- as(findOverlaps(regions, BSseq), "List")
    if (what == "perBase") {
        out <- mclapply(list_of_hits, function(i) {
            out[i, ]
        }, mc.cores = mc.cores)
        return(out)
    } else if (what == "perRegionAverage") {
        out <- mclapply(list_of_hits, function(i) {
            colMeans(out[i, ])
        }, mc.cores = mc.cores)
    } else if (what == "perRegionTotal") {
        out <- mclapply(list_of_hits, function(i) {
            colSums(out[i, ])
        }, mc.cores = mc.cores)
    }
    out <- do.call(rbind, out)
    colnames(out) <- sampleNames(BSseq)
    .DelayedMatrix(out)
}
