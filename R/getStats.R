# TODO: getStats() return a data.frame if regions is specified but a
#       matrix/list if regions is NULL (this is because dmrFinder() calls
#       getStats() and wants the result as a data.frame); should the return
#       class be made consistent (and just implement the coercion in
#       dmrFinder())?
getStats <- function(bstat, regions = NULL, ...) {
    stopifnot(is(bstat, "BSseqTstat") || is(bstat, "BSseqStat"))
    if (is(bstat, "BSseqTstat")) {
        return(getStats_BSseqTstat(bstat, regions = regions, ...))
    }
    getStats_BSseqStat(bstat, regions = regions, ...)
}

.getRegionStats <- function(stat) {
    areaStat <- sum(stat)
    maxStat <- max(stat)
    c(areaStat, maxStat)
}

getStats_BSseqStat <- function(BSseqStat, regions = NULL, what = NULL) {
    stopifnot(is(BSseqStat, "BSseqStat"))
    if (!is.null(what)) {
        stopifnot(what %in% names(BSseqStat@stats))
        return(BSseqStat@stats[[what]])
    }
    if (is.null(regions)) {
        return(BSseqStat@stats)
    }
    ## Now we have regions and no what
    if (class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    ov <- findOverlaps(BSseqStat, regions)
    # NOTE: Rather than split(ov) [as in the original implementation of
    #       getStats_BSseqStat()], we split() the required element of the
    #       @stats slot. This is slightly more efficient if @stat is a matrix
    #       and **much** more efficient if @stat is a DelayedArray.
    #       split,DelayedArray-method returns a *List and thus realises the
    #       data in memory. And, of course, split.array() returns a list, which
    #       is already realised in memory
    stat <- BSseqStat@stats$stat[queryHits(ov), ]
    regionStats <- matrix(NA, ncol = 2, nrow = length(regions),
                          dimnames = list(NULL, c("areaStat", "maxStat")))
    tmp <- lapply(split(stat, subjectHits(ov)), .getRegionStats)
    regionStats[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    regionStats
}

.getRegionStats_ttest <- function(g1m, g2m, tsd) {
    group1.mean <- mean(g1m)
    group2.mean <- mean(g2m)
    meanDiff <- mean(g1m - g2m)
    tstat.sd <- mean(tsd)
    c(meanDiff, group1.mean, group2.mean, tstat.sd)
}

# TODO: Why doesn't this have a 'what' argument? Causes (minor) problems in
#       plotting.R since it requires slightly different handling to F-stat
getStats_BSseqTstat <- function(BSseqTstat, regions = NULL, stat = "tstat.corrected") {
    stopifnot(is(BSseqTstat, "BSseqTstat"))
    stopifnot(stat %in% colnames(BSseqTstat@stats))
    if (is.null(regions)) {
        return(BSseqTstat@stats)
    }
    if (class(regions) == "data.frame") {
        regions <- data.frame2GRanges(regions)
    }
    stopifnot(length(stat) == 1)
    stopifnot(is(regions, "GenomicRanges"))
    # NOTE: Rather than split() ov [as in the original implementation of
    #       getStats_BSseqTstat()], we split() the required columns of the
    #       @stats slot. This is slightly more efficient if @stat is a matrix
    #       and **much** more efficient if @stat is a DelayedArray.
    #       split,DelayedArray-method returns a *List and thus realises the
    #       data in memory. And, of course, split.array() returns a list, which
    #       is already realised in memory
    ov <- findOverlaps(BSseqTstat, regions)
    mat <- BSseqTstat@stats[queryHits(ov), ]
    regionStats <- matrix(NA, ncol = 2, nrow = length(regions),
                          dimnames = list(NULL, c("areaStat", "maxStat")))
    tmp <- lapply(split(mat[, stat], subjectHits(ov)), .getRegionStats)
    regionStats[as.integer(names(tmp)), ] <- do.call(rbind, tmp)
    out <- as.data.frame(regionStats)
    if (!stat %in% c("tstat.corrected", "tstat")) {
        return(out)
    }
    stats_ttest <- matrix(NA, ncol = 4, nrow = length(regions))
    colnames(stats_ttest) <- c("meanDiff", "group1.mean", "group2.mean", "tstat.sd")
    tmp <- mapply(.getRegionStats_ttest,
                  g1m = split(mat[, "group1.means"], subjectHits(ov)),
                  g2m = split(mat[, "group2.means"], subjectHits(ov)),
                  tsd = split(mat[, "tstat.sd"], subjectHits(ov)),
                  SIMPLIFY = FALSE)
    stats_ttest[as.integer(names(tmp)), ] <- do.call(rbind, tmp)
    out <- cbind(out, as.data.frame(stats_ttest))
    out
}
