# NOTE: y should be transposed from its normal orientation, i.e., samples as
#       rows and loci as columns
.bsseq.lm.fit <- function(y, design, contrasts) {
    fit <- .lmFit(y, design)
    contrasts.fit(fit, contrasts)
}

# NOTE: Internal helper used by BSmooth.fstat()
.BSmooth.fstat <- function(tAllPs, parameters, design, contrasts,
                           verbose = TRUE, hdf5 = FALSE) {
    ptime1 <- proc.time()
    # TODO: Experiment with calling .bsseq.lm.fit() on chunks of tAllPs,
    #       realised serially. Should give identical results with reduced
    #       memory footprint
    # NOTE: Currently have to realise all of tAllPs in memory. tAllPs can be
    #       large because it has nrow = number of loci and ncol = number of
    #       samples. It could be useful, and might be necessary, to do this
    #       in such a way that the linear model is iteratively fit to chunks of
    #       tAllPs that are successively realised in memory
    if (is(tAllPs, "DelayedMatrix")) {
        tAllPs <- as.array(tAllPs)
    }
    fitC <- .bsseq.lm.fit(tAllPs, design, contrasts)
    ## Need
    ##   fitC$coefficients, fitC$stdev.unscaled, fitC$sigma, fitC$cov.coefficients
    ## actuall just need
    ##   tstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    ##   rawSds <- fitC$sigma
    ##   cor.coefficients <- cov2cor(fitC$cov.coefficients)
    if (hdf5) {
        rawSds <- .safeHDF5Array(as.matrix(fitC$sigma), "BSseq.", "rawSds")
    } else {
        rawSds <- as.matrix(fitC$sigma)
    }
    cor.coefficients <- cov2cor(fitC$cov.coefficients)
    if (hdf5) {
        rawTstats <- .safeHDF5Array(fitC$coefficients / fitC$stdev.unscaled
                                    / fitC$sigma, "BSseq.", "rawTstats")
    } else {
        rawTstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    }
    names(dimnames(rawTstats)) <- NULL
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f sec\n", stime))
    }

    # NOTE: rawSds and rawTstats are both vectors with length = nrow(BSseq),
    #       cor.coefficients is a n x n matrix, where n = ncol(contrasts);
    #       the first two *may* benefit from being stored as HDF5-backed
    #       column matrices, but the savings are minor since this scales with
    #       nrow(BSseq) rather than ncol(BSseq)
    list(rawTstats = rawTstats,
         rawSds = rawSds,
         cor.coefficients = cor.coefficients)
}

# TODO: The stats matrix *may* benefit from being stored as HDF5-backed matrix,
#       but the saving is minor since this scales with nrow(BSseq) rather than
#       ncol(BSseq); discuss with Kasper.
BSmooth.fstat <- function(BSseq, design, contrasts, verbose = TRUE,
                          hdf5 = FALSE) {
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))

    ## if(any(rowSums(getCoverage(BSseq)[, unlist(groups)]) == 0))
    ##     warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")

    if (verbose) {
        cat("[BSmooth.fstat] fitting linear models ... ")
    }
    tAllPs <- t(getMeth(BSseq, type = "smooth", what = "perBase",
                        confint = FALSE))
    stats <- .BSmooth.fstat(tAllPs = tAllPs,
                            design = design,
                            contrasts = contrasts,
                            verbose = verbose,
                            hdf5 = hdf5)
    parameters <- c(getBSseq(BSseq, "parameters"),
                    list(design = design, contrasts = contrasts))

    BSseqStat(gr = granges(BSseq), stats = stats, parameters = parameters)
}

# NOTE: Internal helper used by smoothSds()
.smoothSds <- function(clusterIdx, rawSds, k, qSd, mc.cores = 1, hdf5 = FALSE) {
    smoothSds <- do.call("c",
                         mclapply(clusterIdx, function(idx) {
                             # NOTE: Need to realise rawSds as an array because
                             #       .smoothSd() works with in-memory data
                             rawSds <- as.array(rawSds[idx, ])
                             .smoothSd(rawSds, k = k, qSd = qSd)
                         }, mc.cores = mc.cores))
    smoothSds <- matrix(smoothSds, ncol = 1)
    if (hdf5) {
        smoothSds <- .safeHDF5Array(smoothSds, "BSseq.", "smoothSds")
    }
    smoothSds
}

# NOTE: hdf5 = TRUE only affects output created by smoothSds(), i.e. the
#       column vector `smoothSds`
smoothSds <- function(BSseqStat, k = 101, qSd = 0.75, maxGap = 10 ^ 8,
                      mc.cores = 1, verbose = TRUE, hdf5 = FALSE) {
    if (is.null(maxGap)) {
        maxGap <- BSseqStat@parameters[["maxGap"]]
    }
    if (is.null(maxGap)) {
        stop("need to set argument 'maxGap'")
    }
    if (verbose) {
        cat("[smoothSds] preprocessing ... ")
    }
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(granges(BSseqStat), maxGap = maxGap)
    rawSds <- getStats(BSseqStat, what = "rawSds")
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        cat(sprintf("done in %.1f sec\n", stime))
    }
    smoothSds <- .smoothSds(clusterIdx = clusterIdx, rawSds = rawSds, k = k,
                            qSd = qSd, mc.cores = mc.cores)
    if ("smoothSds" %in% names(getStats(BSseqStat))) {
        BSseqStat@stats[["smoothSds"]] <- smoothSds
    } else {
        BSseqStat@stats <- c(getStats(BSseqStat), list(smoothSds = smoothSds))
    }
    BSseqStat
}

# NOTE: Required to quieten R CMD check
globalVariables("tstat")

# NOTE: Internal helper used by computeStat()
.computeStat <- function(rawTstats, rawSds, smoothSds, coef, cor.coefficients,
                         hdf5 = FALSE) {
    # TODO: Do I really need to explicitly subset rawSds and smoothSds
    tstats <- rawTstats[, coef, drop = FALSE] * rawSds[, 1L] / smoothSds[, 1L]
    # TODO: Need to realise tstats in case components are DelayedArray objects
    if (length(coef) > 1) {
        # TODO: Need to check this branch
        cor.coefficients <- cor.coefficients[coef, coef]
        stat <- matrix(as.numeric(classifyTestsF(tstats, cor.coefficients,
                                                 fstat.only = TRUE)), ncol = 1)
    } else {
        stat <- tstats
    }
    if (hdf5) {
        stat <- .safeHDF5Array(stat, "BSseq.", stat)
    }
    stat
}

# NOTE: hdf5 = TRUE only affects output created by smoothSds(), i.e. the
#       column vector `smoothSds`
computeStat <- function(BSseqStat, coef = NULL, hdf5 = FALSE) {
    stopifnot(is(BSseqStat, "BSseqStat"))

    rawTstats <- getStats(BSseqStat, what = "rawTstats")
    rawSds <- getStats(BSseqStat, what = "rawSds")
    smoothSds <- getStats(BSseqStat, what = "smoothSds")
    if (is.null(coef)) {
        coef <- seq_len(ncol(rawTstats))
    }
    if (length(coef) > 1) {
        cor.coefficients <- getStats(BSseqStat, what = "cor.coefficients")
        stat.type <- "fstat"
    } else {
        cor.coefficients <- NULL
        stat.type <- "tstat"
    }
    stat <- .computeStat(rawTstats = rawTstats,
                         rawSds = rawSds,
                         smoothSds = smoothSds,
                         coef = coef,
                         cor.coefficients = cor.coefficients,
                         hdf5 = hdf5)

    if ("stat" %in% names(getStats(BSseqStat))) {
        BSseqStat@stats[["stat"]] <- stat
        BSseqStat@stats[["stat.type"]] <- stat.type
    } else {
        BSseqStat@stats <- c(getStats(BSseqStat),
                             list(stat = stat, stat.type = stat.type))
    }
    BSseqStat
}

# TODO (longterm): Support HDF5-backed BSseqStat objects in localCorrectStat().
#                  Not a priority since localCorrectStat() is not yet used
#                  elsewhere in bsseq (which is why it is commented out)
# localCorrectStat <- function(BSseqStat, threshold = c(-15,15), mc.cores = 1, verbose = TRUE) {
#     compute.correction <- function(idx) {
#         xx <- start(BSseqTstat)[idx]
#         yy <- tstat[idx] ## FIXME
#         if(!is.null(threshold)) {
#             stopifnot(is.numeric(threshold) && length(threshold) == 2)
#             stopifnot(threshold[1] < 0 && threshold[2] > 0)
#             yy[yy < threshold[1]] <- threshold[1]
#             yy[yy > threshold[2]] <- threshold[2]
#         }
#         suppressWarnings({
#             drange <- diff(range(xx, na.rm = TRUE))
#         })
#         if(drange <= 25000)
#             return(yy)
#         tstat.function <- approxfun(xx, yy)
#         xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
#         yy.reg <- tstat.function(xx.reg)
#         fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0),
#                       family = "huber", maxk = 50000)
#         correction <- predict(fit, newdata = data.frame(xx.reg = xx))
#         yy - correction
#     }
#     maxGap <- BSseqStat$parameters$maxGap
#     if(verbose) cat("[BSmooth.tstat] preprocessing ... ")
#     ptime1 <- proc.time()
#     clusterIdx <- makeClusters(BSseqStat$gr, maxGap = maxGap)
#     ptime2 <- proc.time()
#     stime <- (ptime2 - ptime1)[3]
#     if(verbose) cat(sprintf("done in %.1f sec\n", stime))
#     stat <- BSseqStat$stat
#     stat.corrected <- do.call(c, mclapply(clusterIdx, compute.correction,
#                                           mc.cores = mc.cores))
#     BSseqTstat@stats <- cbind(getStats(BSseqTstat),
#                               "tstat.corrected" = stat.corrected)
#     BSseqTstat@parameters$local.local <- TRUE
#     BSseqTstat
# }

# NOTE: This uses a lot of 'dot' functions in place of their 'undotted'
#       counterparts. In general, the undotted functions do more processing
#       of their inputs and are consequently slower. In particular, in the
#       permutation setting, they often needlessly re-construct the same
#       objects each permutation. By constructing these objects only once
#       and passing them to the lower-level dot functions, we save a lot of
#       unnecessary computation.
.fstat.dmr.pipeline <- function(gr, tAllPs, parameters, design, contrasts,
                                clusterIdx, coef, cutoff, maxGap, k = 101,
                                qSd = 0.75, return_bstat = TRUE,
                                verbose = TRUE, hdf5 = FALSE) {
    stats <- .BSmooth.fstat(tAllPs = tAllPs,
                            design = design,
                            contrasts = contrasts,
                            hdf5 = hdf5)
    smoothSds <- .smoothSds(clusterIdx = clusterIdx,
                            rawSds = stats[["rawSds"]],
                            k = k,
                            qSd = qSd,
                            mc.cores = 1,
                            hdf5 = hdf5)
    stat <- .computeStat(rawTstats = stats[["rawTstats"]],
                         rawSds = stats[["rawSds"]],
                         smoothSds = smoothSds,
                         coef = coef,
                         cor.coefficients = stats[["cor.coefficients"]],
                         hdf5 = hdf5)
    dmrs <- .dmrFinder(dmrStat = stat,
                       cutoff = cutoff,
                       seqnames_as_char = as.character(seqnames(gr)),
                       positions = start(gr),
                       maxGap = maxGap,
                       verbose = max(as.integer(verbose) - 1L, 0L))
    if (!is.null(dmrs)) {
        ov <- findOverlaps(gr, data.frame2GRanges(dmrs))
        dmrs_stats <- .getRegionStats_BSseqStat(ov = ov, stat = stat)
        dmrs <- cbind(dmrs, dmrs_stats)
        dmrs <- dmrs[order(abs(dmrs[["areaStat"]]), decreasing = TRUE), ]
    }
    if (return_bstat) {
        parameters <- c(parameters,
                        list(design = design, contrasts = contrasts))
        bstat <- BSseqStat(gr = gr, stats = stats, parameters = parameters)
    } else {
        bstat <- NULL
    }
    list(bstat = bstat, dmrs = dmrs)
}

# TODO: How should the hdf5 argument work? Probably don't want to write all
#       outputs to disk
fstat.pipeline <- function(BSseq, design, contrasts, cutoff, fac, nperm = 1000,
                           coef = NULL, maxGap.sd = 10 ^ 8, maxGap.dmr = 300,
                           type = "dmrs", mc.cores = 1, verbose = TRUE,
                           hdf5 = FALSE) {
    type <- match.arg(type, c("dmrs", "blocks"))
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))

    permutationMatrix <- permuteAll(nperm, design)
    if (nrow(permutationMatrix) < nperm) {
        warning(paste0("Only ", nrow(permutationMatrix), " unique ",
                       "permutations exist (requested ", nperm),
                " permutations)")
        nperm <- nrow(permutationMatrix)
    }

    # NOTE: Certain objects are used in identifying the candidate DMRs and can
    #       be reused without changed when identifying DMRs in the permuted
    #       data. Constructing these once saves unnecessary computation
    gr <- rowRanges(BSseq)
    # TODO: Should tAllPs be realised at this point; i.e. can the forked
    #       processed share the same tAllPs object?
    tAllPs <- t(getMeth(BSseq, type = "smooth", what = "perBase",
                        confint = FALSE))
    parameters <- getBSseq(BSseq, "parameters")
    clusterIdx <- makeClusters(gr = gr, maxGap = maxGap.sd)
    if (is.null(coef)) {
        coef <- seq_len(ncol(design) - 1L)
    }

    # Compute bstat and identify candidate DMRs
    bstat_and_dmrs <- .fstat.dmr.pipeline(gr = gr,
                                          tAllPs = tAllPs,
                                          parameters = parameters,
                                          design = design,
                                          contrasts = contrasts,
                                          clusterIdx = clusterIdx,
                                          coef = coef,
                                          cutoff = cutoff,
                                          maxGap = maxGap.dmr,
                                          return_bstat = TRUE,
                                          verbose = verbose,
                                          hdf5 = hdf5)
    bstat <- bstat_and_dmrs[["bstat"]]
    dmrs <- bstat_and_dmrs[["dmrs"]]

    if (is.null(dmrs)) {
        stop("No DMRs identified. Consider reducing the 'cutoff' from (",
             paste0(cutoff, collapse = ", "), ")")
    }

    # Identify DMRs in permutated data
    # NOTE: More efficient to permute design matrix using idxMatrix[ii, ]
    #       than to permute the raw data with tAllPs[idxMatrix[ii, ]]
    nullDist <- .getNullDistribution_BSmooth.fstat(
        permutationMatrix = permutationMatrix,
        gr = gr,
        tAllPs = tAllPs,
        design = design,
        contrasts = contrasts,
        clusterIdx = clusterIdx,
        coef = coef,
        cutoff = cutoff,
        maxGap = maxGap.dmr,
        mc.cores = mc.cores)

    # Compute FWER for candidate DMRs
    fwer <- getFWER.fstat(null = c(list(dmrs), nullDist), type = type)
    dmrs$fwer <- fwer

    # Compute average methylation level in each group (`fac`) for each
    # candidate DMR
    meth <- getMeth(BSseq, regions = dmrs, what = "perRegion",
                    mc.cores = mc.cores)
    meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
    dmrs <- cbind(dmrs, meth)

    # Compute the maximum difference in group-wise mean methylation levels for
    # each candidate DMR
    dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)

    # Return BSseqStat object computed using unpermuted data, DMRs with FWER,
    # the permutation matrix, and the DMRs in the permuted data
    list(bstat = bstat, dmrs = dmrs, permutationMatrix = permutationMatrix,
         nullDist = nullDist)
}

fstat.comparisons.pipeline <- function(BSseq, design, contrasts, cutoff, fac,
                                       maxGap.sd = 10 ^ 8, maxGap.dmr = 300,
                                       verbose = TRUE) {
    bstat <- BSmooth.fstat(BSseq = BSseq,
                           design = design,
                           contrasts = contrasts,
                           verbose = verbose)
    bstat <- smoothSds(bstat, maxGap = maxGap.sd, verbose = verbose)
    # NOTE: Want to keep the bstat object corresponding to the original fstat
    #       and not that from the last t-tsts in the following lapply()
    bstat_f <- bstat
    if (verbose) {
        cat(paste0("[fstat.comparisons.pipeline] making ", ncol(contrasts),
                   " comparisons ... \n"))
    }
    l <- lapply(seq_len(ncol(contrasts)), function(coef) {
        if (verbose) {
            cat(paste0("[fstat.comparisons.pipeline] ",
                       colnames(contrasts)[coef], "\n"))
        }
        bstat <- computeStat(bstat, coef = coef)
        dmrs <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr,
                          verbose = verbose)
        if (!is.null(dmrs)) {
            meth <- getMeth(BSseq, dmrs, what = "perRegion")
            meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
            dmrs <- cbind(dmrs, meth)
            dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)
        }
        dmrs
    })
    names(l) <- colnames(contrasts)
    list(bstat = bstat, dmrs = l)
}

