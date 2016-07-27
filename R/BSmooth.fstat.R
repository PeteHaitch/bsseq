# TODO: The stats matrix *may* benefit from being stored as HDF5-backed matrix,
#       but the saving is minor since this scales with nrow(BSseq) rather than
#       ncol(BSseq); discuss with Kasper.
BSmooth.fstat <- function(BSseq, design, contrasts, verbose = TRUE,
                          hdf5 = FALSE) {
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))

    ## if(any(rowSums(getCoverage(BSseq)[, unlist(groups)]) == 0))
    ##     warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")

    if(verbose) cat("[BSmooth.fstat] fitting linear models ... ")
    ptime1 <- proc.time()
    allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                     confint = FALSE)
    if (is(allPs, "DelayedMatrix")) {
        # NOTE: Need to realise `allPs` as an array since we use
        #       limma::lmFit(allPs, design); actually, lmFit will do this
        #       internally via a call to as.matrix, but we explicitly do this
        #       ourselves to protect ourselves against changes to the internals
        #       of lmFit().
        allPs <- as.array(allPs)
    }

    fit <- lmFit(allPs, design)
    fitC <- contrasts.fit(fit, contrasts)
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
        rawTstats <-
            .safeHDF5Array(fitC$coefficients / fitC$stdev.unscaled / fitC$sigma,
                           "BSseq.", "rawTstats")
    } else {
        rawTstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    }
    names(dimnames(rawTstats)) <- NULL
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))

    parameters <- c(BSseq@parameters,
                    list(design = design, contrasts = contrasts))
    # NOTE: rawSds and rawTstats are both vectors with length = nrow(BSseq),
    #       cor.coefficients is a n x n matrix, where n = ncol(contrasts);
    #       the first two *may* benefit from being stored as HDF5-backed
    #       column matrices, but the savings are minor since this scales with
    #       nrow(BSseq) rather than ncol(BSseq)
    stats <- list(rawSds = rawSds,
                  cor.coefficients = cor.coefficients,
                  rawTstats = rawTstats)
    out <- BSseqStat(gr = rowRanges(BSseq),
                     stats = stats, parameters = parameters)
    out
}

# NOTE: hdf5 = TRUE only affects output created by smoothSds(), i.e. the
#       column vector `smoothSds`
smoothSds <- function(BSseqStat, k = 101, qSd = 0.75, mc.cores = 1,
                      maxGap = 10^8, verbose = TRUE, hdf5 = FALSE) {
    if(is.null(maxGap))
        maxGap <- BSseqStat@parameters[["maxGap"]]
    if(is.null(maxGap))
        stop("need to set argument 'maxGap'")
    if(verbose) cat("[smoothSds] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(granges(BSseqStat), maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    smoothSds <- do.call("c",
                         mclapply(clusterIdx, function(idx) {
                             # NOTE: Need to realise rawSds as an array because
                             #       .smoothSds() works with in-memory data
                             rawSds <- as.array(getStats(BSseqStat, what = "rawSds")[idx, ])
                             .smoothSd(rawSds, k = k, qSd = qSd)
                         }, mc.cores = mc.cores))
    if (hdf5) {
        smoothSds <- .safeHDF5Array(as.matrix(smoothSds), "BSseq.", "smoothSds")
    } else {
        smoothSds <- as.matrix(smoothSds)
    }
    if("smoothSds" %in% names(getStats(BSseqStat)))
        BSseqStat@stats[["smoothSds"]] <- smoothSds
    else
        BSseqStat@stats <- c(getStats(BSseqStat), list(smoothSds = smoothSds))
    BSseqStat
}

# NOTE: Required to quieten R CMD check
globalVariables("tstat")

# NOTE: hdf5 = TRUE only affects output created by smoothSds(), i.e. the
#       column vector `smoothSds`
computeStat <- function(BSseqStat, coef = NULL, hdf5 = FALSE) {
    stopifnot(is(BSseqStat, "BSseqStat"))
    if(is.null(coef)) {
        coef <- 1:ncol(getStats(BSseqStat, what = "rawTstats"))
    }
    tstats <- getStats(BSseqStat, what = "rawTstats")[, coef, drop = FALSE]
    tstats <- tstats * getStats(BSseqStat, what = "rawSds") /
        getStats(BSseqStat, what = "smoothSds")
    if(length(coef) > 1) {
        # TODO: Need to check this branch
        cor.coefficients <- getStats(BSseqStat, what = "cor.coefficients")[coef,coef]
        stat <- as.numeric(classifyTestsF(tstats, cor.coefficients,
                                          fstat.only = TRUE))
        stat.type <- "fstat"
    } else {
        stat <- tstats
        stat.type <- "tstat"
    }
    if (hdf5) {
        stat <- .safeHDF5Array(stat, "BSseq.", stat)
    }
    if("stat" %in% names(getStats(BSseqStat))) {
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

fstat.pipeline <- function(BSseq, design, contrasts, cutoff, fac, nperm = 1000,
                           coef = NULL, maxGap.sd = 10 ^ 8, maxGap.dmr = 300,
                           mc.cores = 1, hdf5 = FALSE) {
    bstat <- BSmooth.fstat(BSseq = BSseq, design = design,
                           contrasts = contrasts, hdf5 = hdf5)
    bstat <- smoothSds(bstat, hdf5 = hdf5)
    bstat <- computeStat(bstat, coef = coef, hdf5 = hdf5)
    dmrs <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr)
    if (is.null(dmrs)) {
        stop("No DMRs identified. Consider reducing the 'cutoff' from (",
             paste0(cutoff, collapse = ", "), ")")
    }
    idxMatrix <- permuteAll(nperm, design)
    nullDist <- getNullDistribution_BSmooth.fstat(BSseq = BSseq,
                                                  idxMatrix = idxMatrix,
                                                  design = design,
                                                  contrasts = contrasts,
                                                  coef = coef,
                                                  cutoff = cutoff,
                                                  maxGap.sd = maxGap.sd,
                                                  maxGap.dmr = maxGap.dmr,
                                                  mc.cores = mc.cores)
    fwer <- getFWER.fstat(null = c(list(dmrs), nullDist), type = "dmrs")
    dmrs$fwer <- fwer
    meth <- getMeth(BSseq, regions = dmrs, what = "perRegion")
    meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
    dmrs <- cbind(dmrs, meth)
    dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)
    list(bstat = bstat, dmrs = dmrs, idxMatrix = idxMatrix, nullDist = nullDist)
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

