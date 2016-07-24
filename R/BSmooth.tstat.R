# TODO: Should hdf5 be the default? The `stats` slot dominates the object's
#       size (nrow = number of loci, ncol = 5 or 6). Should it guess/recommend
#       the user switch to hdf5 nrow(BSseq) is large? If `hdf5 = TRUE` then
#       each column of `stats` is in its own .h5 file; this reduces the
#       memory usage of BSmooth.tstat() at the expense of some additional
#       running time.
# TODO: The real memory issue is that `allPs` (nrow = length(BSseq),
#       ncol = number of samples) is also realised in memory during the
#       function call; see TODO below.
# NOTE: hdf5 = TRUE only affects output created by BSmooth.tstat(), i.e. the
#       matrix-like object in the `stats` slot
# TODO: The stats matrix *may* benefit from being stored as HDF5-backed matrix,
#       but the saving is minor since this scales with nrow(BSseq) rather than
#       ncol(BSseq); discuss with Kasper.
BSmooth.tstat <- function(BSseq, group1, group2,
                          estimate.var = c("same", "paired", "group2"),
                          local.correct = TRUE, maxGap = NULL, qSd = 0.75,
                          k = 101, mc.cores = 1, verbose = TRUE, hdf5 = FALSE) {
    compute.correction <- function(idx, tstat, qSd = 0.75) {
        xx <- start(BSseq)[idx]
        yy <- tstat[idx]
        suppressWarnings({
            drange <- diff(range(xx, na.rm = TRUE))
        })
        if (drange <= 25000)
            return(yy)
        tstat.function <- approxfun(xx, yy)
        xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
        yy.reg <- tstat.function(xx.reg)
        fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0),
                      family = "huber", maxk = 50000)
        correction <- predict(fit, newdata = data.frame(xx.reg = xx))
        yy - correction
    }
    estimate.var <- match.arg(estimate.var)
    stopifnot(is(BSseq, "BSseq."))
    stopifnot(hasBeenSmoothed(BSseq))
    if (is.character(group1)) {
        stopifnot(all(group1 %in% sampleNames(BSseq)))
        group1 <- match(group1, sampleNames(BSseq))
    }
    if (is.numeric(group1)) {
        stopifnot(min(group1) >= 1 & max(group1) <= ncol(BSseq))
    } else stop("problems with argument 'group1'")
    if (is.character(group2)) {
        stopifnot(all(group2 %in% sampleNames(BSseq)))
        group2 <- match(group2, sampleNames(BSseq))
    }
    if (is.numeric(group2)) {
        stopifnot(min(group2) >= 1 & max(group2) <= ncol(BSseq))
    } else {
        stop("problems with argument 'group2'")
    }
    stopifnot(length(intersect(group1, group2)) == 0)
    stopifnot(length(group1) > 0)
    stopifnot(length(group2) > 0)
    stopifnot(length(group1) + length(group2) >= 3)
    if (estimate.var == "paired")
        stopifnot(length(group1) == length(group2))

    # TODO: rowSums(DelayedMatrix) can be very slow if @index slot has a large
    #       amount of i-subsetting. This is currently unavoidable (because of
    #       how the HDF5 file(s) are read from disk) but should at least be
    #       documented.
    if (any(rowSums(getCoverage(BSseq)[, c(group1, group2)]) == 0)) {
        warning("Computing t-statistics at locations where there is no data; ",
                "consider subsetting the 'BSseq' object first")
    }

    if (is.null(maxGap)) {
        maxGap <- BSseq@parameters$maxGap
    }
    if (is.null(maxGap)) {
        stop("need to set argument 'maxGap'")
    }
    if (verbose) cat("[BSmooth.tstat] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(BSseq, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) cat(sprintf("done in %.1f sec\n", stime))

    if (verbose) cat("[BSmooth.tstat] computing stats within groups ... ")
    ptime1 <- proc.time()
    allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                     confint = FALSE)
    if (is(allPs, "DelayedMatrix")) {
        # NOTE: Currently, need to realise allPs as an array since we use
        #       matrixStats::rowSds(allPs) or matrixStats::rowVars(allPs);
        # TODO: Need to make these matrixStats functions 'block processing
        #       friendly' to avoid realising the entire thing in memory in one
        #       go. Remember, dim(allPs) = dim(BSseq)!
        allPs <- as.array(allPs)
    }
    cn <- c("rawSds", "tstat.sd", "group2.means", "group1.means", "tstat")
    if (local.correct) {
        cn <- c(cn, "tstat.corrected")
    }
    if (hdf5) {
        group1.means <- .safeHDF5Array(
            as.matrix(rowMeans(allPs[, group1, drop = FALSE], na.rm = TRUE)),
            "BSseq.", "group1.means")
        group2.means <- .safeHDF5Array(
            as.matrix(rowMeans(allPs[, group2, drop = FALSE], na.rm = TRUE)),
            "BSseq.", "groups2.means")
    } else {
        # NOTE: Preallocate stats to avoid copies that would otherwise
        #       occur when cbind()-ing vectors to form stats
        stats <- matrix(0, nrow = length(BSseq), ncol = length(cn),
                        dimnames = list(NULL, cn))
        stats[, "group1.means"] <- rowMeans(allPs[, group1, drop = FALSE],
                                            na.rm = TRUE)
        stats[, "group2.means"] <- rowMeans(allPs[, group2, drop = FALSE],
                                        na.rm = TRUE)
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) cat(sprintf("done in %.1f sec\n", stime))

    if (verbose) cat("[BSmooth.tstat] computing stats across groups ... ")
    ptime1 <- proc.time()
    # NOTE: rawSds is only realised on disk after running .smoothSd()
    #       because .smoothSd() requires rawSds be realised in memory
    switch(estimate.var,
           "group2" = {
               if (hdf5) {
                   rawSds <- rowSds(allPs, cols = group2, na.rm = TRUE)

               } else {
                   stats[, "rawSds"] <- rowSds(allPs, cols = group2,
                                               na.rm = TRUE)
               }
               scale <- sqrt(1/length(group1) + 1/length(group2))
           },
           "same" = {
               if (hdf5) {
                   rawSds <- sqrt(((length(group1) - 1) *
                                       rowVars(allPs, cols = group1) +
                                       (length(group2) - 1) *
                                       rowVars(allPs, cols = group2)) /
                                      (length(group1) + length(group2) - 2))
               } else {
                   stats[, "rawSds"] <-
                       sqrt(((length(group1) - 1) *
                                 rowVars(allPs, cols = group1) +
                                 (length(group2) - 1) *
                                 rowVars(allPs, cols = group2)) /
                                (length(group1) + length(group2) - 2))
               }
               scale <- sqrt(1/length(group1) + 1/length(group2))
           },
           "paired" = {
               if (hdf5) {
                   rawSds <- rowSds(allPs[, group1, drop = FALSE] -
                                        allPs[, group2, drop = FALSE])
               } else {
                   stats[, "rawSds"] <-
                       rowSds(allPs[, group1, drop = FALSE] -
                                  allPs[, group2, drop = FALSE])
               }
               scale <- sqrt(1/length(group1))
           })
    if (hdf5) {
        tstat.sd <- .safeHDF5Array(
            as.matrix(do.call(c, mclapply(clusterIdx, function(idx) {
            scale * .smoothSd(rawSds[idx], k = k, qSd = qSd)
        }, mc.cores = mc.cores))),
        "BSseq.", "tstat.sd")
        rawSds <- .safeHDF5Array(as.matrix(rawSds), "BSseq.", "rawSds")
        tstat <- .safeHDF5Array((group1.means - group2.means) / tstat.sd,
                               "BSseq.", "tstat")
    } else {
        stats[, "tstat.sd"] <- do.call(c, mclapply(clusterIdx, function(idx) {
            scale * .smoothSd(stats[idx, "rawSds"], k = k, qSd = qSd)
        }, mc.cores = mc.cores))
        stats[, "tstat"] <- (stats[, "group1.means"] -
                                 stats[, "group2.means"]) / stats[, "tstat.sd"]
    }

    # TODO: This line doesn't do anything; I think it's meant to be
    #       tstat[is.na(tstat)[tstat.sd == 0]] <- TRUE or perhaps
    #       tstat[is.na(tstat) & tstat.sd == 0] <- TRUE
    # is.na(tstat)[tstat.sd == 0] <- TRUE
    if (local.correct) {
        if (hdf5) {
            tstat.corrected <- .safeHDF5Array(
                as.matrix(do.call(c, mclapply(clusterIdx,
                                              compute.correction,
                                              tstat = as.array(tstat),
                                              qSd = qSd,
                                              mc.cores = mc.cores))),
                "BSseq.", "tstat.corrected")
        } else {
            stats[, "tstat.corrected"] <-
                do.call(c, mclapply(clusterIdx,
                                    compute.correction,
                                    tstat = stats[, "tstat"],
                                    qSd = qSd,
                                    mc.cores = mc.cores))
        }
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) cat(sprintf("done in %.1f sec\n", stime))

    if (hdf5) {
        stats <- cbind(rawSds, tstat.sd, group2.means, group1.means, tstat)
        if (local.correct) {
            stats <- cbind(stats, tstat.corrected)
        }
        colnames(stats) <- cn
    }

    parameters <- c(BSseq@parameters,
                    list(tstatText = sprintf("BSmooth.tstat (local.correct = %s, maxGap = %d)",
                                             local.correct, maxGap),
                         group1 = group1, group2 = group2, k = k, qSd = qSd,
                         local.correct = local.correct, maxGap = maxGap))
    out <- BSseqTstat(gr = granges(BSseq),
                      stats = stats,
                      parameters = parameters)
    out
}
