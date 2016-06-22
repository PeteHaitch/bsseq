# DONE
getStats <- function(bstat, regions = NULL, ...) {
    stopifnot(is(bstat, "BSseqTstat") || is(bstat, "BSseqStat"))
    if(is(bstat, "BSseqTstat"))
        return(getStats_BSseqTstat(bstat, regions = regions, ...))
    getStats_BSseqStat(bstat, regions = regions, ...)
}

# TODO
getStats_BSseqStat <- function(BSseqStat, regions = NULL, what = NULL) {
    stopifnot(is(BSseqStat, "BSseqStat"))
    if(!is.null(what)) {
        stopifnot(what %in% names(BSseqStat@stats))
        return(BSseqStat@stats[[what]])
    }
    if(is.null(regions))
        return(BSseqStat@stats)
    ## Now we have regions and no what
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    ov <- findOverlaps(BSseqStat, regions)
    ov.sp <- split(queryHits(ov), subjectHits(ov))
    ## We need areaStat and maxStat
    ## Could need averageMeth in each group?
    ## Need to have a specific design for that
    ## Could at least get contrast coefficient
    getRegionStats <- function(idx) {
        areaStat <- sum(BSseqStat@stats$stat[idx], na.rm = TRUE)
        maxStat <- max(BSseqStat@stats$stat[idx], na.rm = TRUE)
        c(areaStat, maxStat)
    }
    regionStats <- matrix(NA, ncol = 2, nrow = length(regions))
    colnames(regionStats) <- c("areaStat", "maxStat")
    tmp <- lapply(ov.sp, getRegionStats)
    regionStats[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    regionStats
}

# DONE
getStats_BSseqTstat <- function(BSseqTstat, regions = NULL, stat = "tstat.corrected") {
    stopifnot(is(BSseqTstat, "BSseqTstat"))
    if(is.null(regions))
        return(BSseqTstat@stats)
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    stopifnot(stat %in% colnames(BSseqTstat@stats))
    stopifnot(length(stat) == 1)
    stopifnot(is(regions, "GenomicRanges"))
    # NOTE: Rather than split() ov, we split() the required columns of the
    #       @stat slot. This is slightly more efficient if @stat is a matrix
    #       and **much** more efficient if @stat is a DelayedArray
    ov <- findOverlaps(BSseqTstat, regions)
    mat <- BSseqTstat@stats[queryHits(ov), ]
    mat_stat <- mat[, stat]
    # NOTE: split,DelayedArray-method returns a *List and thus realises the
    #       data in memory. And of course, split.array() returns a list, which
    #       is already realised in memory
    mat_stat_split <- split(mat_stat, subjectHits(ov))
    getRegionStats <- function(stat) {
        areaStat <- sum(stat)
        maxStat <- max(stat)
        c(areaStat, maxStat)
    }
    stats <- matrix(NA, ncol = 2, nrow = length(regions))
    colnames(stats) <- c("areaStat", "maxStat")
    tmp <- lapply(mat_stat_split, getRegionStats)
    stats[as.integer(names(tmp)),] <- do.call(rbind, tmp)
    out <- as.data.frame(stats)
    if(! stat %in% c("tstat.corrected", "tstat"))
        return(out)
    mat_group1.means <- BSseqTstat@stats[queryHits(ov), "group1.means"]
    mat_group1.means_split <- split(mat_group1.means, subjectHits(ov))
    mat_group2.means <- BSseqTstat@stats[queryHits(ov), "group2.means"]
    mat_group2.means_split <- split(mat_group2.means, subjectHits(ov))
    mat_tstat.sd <- BSseqTstat@stats[queryHits(ov), "tstat.sd"]
    mat_tstat.sd_split <- split(mat_tstat.sd, subjectHits(ov))
    getRegionStats_ttest <- function(g1m, g2m, tsd) {
        group1.mean <- mean(g1m)
        group2.mean <- mean(g2m)
        meanDiff <- mean(g1m - g2m)
        tstat.sd <- mean(tsd)
        c(meanDiff, group1.mean, group2.mean, tstat.sd)
    }
    stats_ttest <- matrix(NA, ncol = 4, nrow = length(regions))
    colnames(stats_ttest) <- c("meanDiff", "group1.mean", "group2.mean", "tstat.sd")
    tmp <- mapply(getRegionStats_ttest,
                  g1m = split(mat_group1.means, subjectHits(ov)),
                  g2m = split(mat_group2.means, subjectHits(ov)),
                  tsd = split(mat_tstat.sd, subjectHits(ov)),
                  SIMPLIFY = FALSE)
    stats_ttest[as.integer(names(tmp)), ] <- do.call(rbind, tmp)
    out <- cbind(out, as.data.frame(stats_ttest))
    out
}
