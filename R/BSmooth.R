makeClusters <- function(hasGRanges, maxGap = 10 ^ 8) {
    chrOrder <- as.character(runValue(seqnames(hasGRanges)))
    if (anyDuplicated(chrOrder)) {
        stop("argument 'hasGRanges' is not properly order")
    }
    grBase <- granges(hasGRanges)
    clusters <- reduce(resize(grBase, width = 2 * maxGap + 1, fix = "center"))
    start(clusters) <- pmax(rep(1, length(clusters)), start(clusters))
    clusters.sp <- split(clusters, seqnames(clusters))
    # Are the clusters ordered within the chromosome? This is probably
    # guaranteed
    stopifnot(all(vapply(clusters.sp, function(cluster.gr) {
        if (length(cluster.gr) <= 1) {
            return(TRUE)
        } else {
            all(start(cluster.gr)[-length(cluster.gr)] < end(cluster.gr)[-1])
        }
    }, logical(1L))))
    clusters <- Reduce(c, clusters.sp[chrOrder])
    stopifnot(all(chrOrder == runValue(seqnames(clusters))))
    ov <- findOverlaps(grBase, clusters)
    clusterIdx <- split(as.matrix(ov)[,1], as.matrix(ov)[,2])
    names(clusterIdx) <- NULL
    clusterIdx
}

# TODO: Should hdf5 be the default? Should it guess/recommend the user switch
#       to hdf5 if they have large/many files?
# NOTE: hdf5 = TRUE only affects assays created by BSmooth(), i.e. the `coef`
#       and `se.coef` assays, and not existing assays such as `M` and `Cov`
BSmooth <- function(BSseq, ns = 70, h = 1000, maxGap = 10^8,
                    parallelBy = c("sample", "chromosome"),
                    mc.preschedule = FALSE, mc.cores = 1, keep.se = FALSE,
                    verbose = TRUE, hdf5 = FALSE) {
    # NOTE: .smooth() realises M and Cov as array objects on a per-sample or
    #       per-chromosome basis (depending on value of `parallelBy` in call to
    #       BSmooth())
    .smooth <- function(idxes, sampleIdx) {
        ## Assuming that idxes is a set of indexes into the BSseq object
        ## sampleIdx is a single character
        this_sample_chr <- c(sampleNames(BSseq)[sampleIdx],
                             as.character(seqnames(BSseq)[idxes[1]]))
        if (verbose >= 2) {
            cat(sprintf("[BSmooth]   smoothing start: sample:%s, chr:%s, nLoci:%s\n",
                        this_sample_chr[1], this_sample_chr[2], length(idxes)))
        }
        Cov <- getCoverage(BSseq, type = "Cov")[idxes, sampleIdx, drop = FALSE]
        M <- getCoverage(BSseq, type = "M")[idxes, sampleIdx, drop = FALSE]
        pos <- start(BSseq)[idxes]
        stopifnot(all(diff(pos) > 0))
        wh <- which(Cov != 0)
        nn <- ns / length(wh)
        if (length(wh) <= ns) {
            if (keep.se) {
                se.coef <- rep(NA_real_, length(Cov))
            } else {
                se.coef <- NULL
            }
            return(list(coef = rep(NA_real_, length(Cov)),
                        se.coef = se.coef,
                        trans = NULL, h = h, nn = nn))
        }
        # NOTE: At this point, it is much, much faster to realise M and Cov in
        #       memory than it is to use delayed operations. Also, at this point,
        #       M and Cov only contain data from a single sample/seqlevel, so are
        #       generally quite a bit smaller than the full object
        if (is(M, "DelayedArray") ) {
            M <- as.array(M)
        }
        if (is(Cov, "DelayedArray")) {
            Cov <- as.array(Cov)
        }
        sdata <- data.frame(pos = pos[wh],
                            M = pmin(pmax(M[wh], 0.01), Cov[wh] - 0.01),
                            Cov = Cov[wh])
        fit <- locfit(M ~ lp(pos, nn = nn, h = h), data = sdata,
                      weights = Cov, family = "binomial", maxk = 10000)
        pp <- preplot(fit, where = "data", band = "local",
                      newdata = data.frame(pos = pos))
        if (keep.se) {
            se.coef <- pp$se.fit
        } else {
            se.coef <- NULL
        }
        if (verbose >= 2) {
            cat(sprintf("[BSmooth]   smoothing end: sample:%s, chr:%s, nLoci:%s, nCoveredLoci:%s\n",
                        this_sample_chr[1], this_sample_chr[2], length(idxes),
                        nrow(sdata)))
        }
        # UP TO HERE: .smooth() could inherit `hdf5` and write `coef`` and
        #             `se.coef` to disk to reduce memory footprint
        return(list(coef = pp$fit, se.coef = se.coef, trans = pp$trans, h = h,
                    nn = nn))
    }
    stopifnot(class(BSseq) == "BSseq")
    parallelBy <- match.arg(parallelBy)
    if (verbose) {
        cat("[BSmooth] preprocessing ... ")
    }
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(BSseq, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f sec\n", stime))
    }
    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames

    ptime.outer1 <- proc.time()
    switch(parallelBy, "sample" = {
        if (verbose) {
            cat(sprintf("[BSmooth] smoothing by 'sample' (mc.cores = %d, mc.preschedule = %s)\n",
                                mc.cores, mc.preschedule))
        }
        out <- mclapply(seq_along(sampleNames), function(sIdx) {
            ptime1 <- proc.time()
            tmp <- lapply(clusterIdx, function(jj) {
                try(.smooth(idxes = jj, sampleIdx = sIdx))
            })
            coef <- do.call(c, lapply(tmp, function(xx) xx$coef))
            if (!is.null(coef)) {
                # UP TO HERE: Testing 2016-07-25 re writing .h5 files once vs. often
                if (hdf5) {
                    coef <- .safeHDF5Array(coef, "BSseq.", "coef")
                }
            }
            se.coef <- do.call(c, lapply(tmp, function(xx) xx$se.coef))
            if (!is.null(se.coef)) {
                # UP TO HERE: Testing 2016-07-25 re writing .h5 files once vs. often
                if (hdf5) {
                    se.coef <- .safeHDF5Array(coef, "BSseq.", "se.coef")
                }
            }
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                cat(sprintf("[BSmooth] sample %s (out of %d), done in %.1f sec\n",
                            sampleNames[sIdx], length(sampleNames), stime))
            }
            return(list(coef = coef, se.coef = se.coef))
        }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
        if (any(vapply(out, is, logical(1L), class2 = "try-error"))) {
            stop("BSmooth encountered smoothing errors")
        }
        coef <- do.call(cbind, lapply(out, function(xx) xx$coef))
        se.coef <- do.call(cbind, lapply(out, function(xx) xx$se.coef))
    }, "chromosome" = {
        if (verbose) {
            cat(sprintf("[BSmooth] smoothing by 'chromosome' (mc.cores = %d, mc.preschedule = %s)\n",
                                mc.cores, mc.preschedule))
        }
        out <- mclapply(seq_along(clusterIdx), function(ii) {
            ptime1 <- proc.time()
            tmp <- lapply(seq(along = sampleNames), function(sIdx) {
                .smooth(idxes = clusterIdx[[ii]], sampleIdx = sIdx)
            })
            coef <- do.call(cbind, lapply(tmp, function(xx) xx$coef))
            if (!is.null(coef)) {
                # UP TO HERE: Testing 2016-07-25 re writing .h5 files once vs. often
                if (hdf5) {
                    coef <- .safeHDF5Array(coef, "BSseq.", "coef")
                }
            }
            se.coef <- do.call(cbind, lapply(tmp, function(xx) xx$se.coef))
            if (!is.null(se.coef)) {
                # UP TO HERE: Testing 2016-07-25 re writing .h5 files once vs. often
                if (hdf5) {
                    se.coef <- .safeHDF5Array(coef, "BSseq.", "se.coef")
                }
            }
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if (verbose) {
                cat(sprintf("[BSmooth] chr idx %d (out of %d), done in %.1f sec\n",
                            ii, length(clusterIdx), stime))
            }
            return(list(coef = coef, se.coef = se.coef))
        }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
        if (any(sapply(out, is, class2 = "try-error"))) {
            stop("BSmooth encountered smoothing errors")
        }
        coef <- do.call(rbind, lapply(out, function(xx) xx$coef))
        se.coef <- do.call(rbind, lapply(out, function(xx) xx$se.coef))
    })
    ptime.outer2 <- proc.time()
    stime.outer <- (ptime.outer2 - ptime.outer1)[3]
    if (verbose) {
        cat(sprintf("[BSmooth] smoothing done in %.1f sec\n", stime.outer))
    }
    rownames(coef) <- NULL
    colnames(coef) <- sampleNames(BSseq)
    if (!is.null(se.coef)) {
        rownames(se.coef) <- NULL
        colnames(se.coef) <- sampleNames(BSseq)
    }

    if (!is.null(coef)) {
        # UP TO HERE: Testing 2016-07-25 re writing .h5 files once vs. often
        # if (hdf5) {
        #     coef <- .safeHDF5Array(coef, "BSseq.", "coef")
        # }
        assay(BSseq, "coef") <- coef
    }
    if (!is.null(se.coef)) {
        # UP TO HERE: Testing 2016-07-25 re writing .h5 files once vs. often
        # if (hdf5) {
        #     se.coef <- .safeHDF5Array(se.coef, "BSseq.", "se.coef")
        # }
        assay(BSseq, "se.coef") <- se.coef
    }
    # NOTE: The commented out version of mytrans() is equivalent to the new one.
    #       However, the new version is approximately 2x as fast and uses less
    #       memory since there is no need to create the index vectors `ix` and
    #       `ix2`. Moreover, there is now a plogis,DelayedArray-method in
    #       HDF5Array
    # mytrans <- function(x) {
    #     y <- x
    #     # TODO: If can't use plogis(), then use x[x < 0] and x[x > 0] instead
    #     #         of creating indices with which(); this requires
    #     #         1D-subsetting of HDF5Array objects (see email from Hervé on
    #     #         2016-06-19).
    #     ix <- which(x < 0)
    #     ix2 <- which(x >= 0)
    #     y[ix] <- exp(x[ix])/(1 + exp(x[ix]))
    #     y[ix2] <- 1/(1 + exp(-x[ix2]))
    #     y
    # }
    # TODO: This will require that all existing BSseq objects are updated using
    #       the new updateObject()
    mytrans <- plogis
    environment(mytrans) <- baseenv()
    BSseq@trans <- mytrans
    parameters <- list(smoothText =
                           sprintf("BSmooth (ns = %d, h = %d, maxGap = %d)",
                                   ns, h, maxGap),
                       ns = ns, h = h, maxGap = maxGap)
    BSseq@parameters <- parameters
    BSseq
}
