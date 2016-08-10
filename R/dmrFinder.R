# NOTE: Internal helper used by dmrFinder()
.dmrFinder <- function(dmrStat, cutoff, seqnames_as_char, positions, maxGap,
                       verbose = FALSE) {
    direction_logical <- dmrStat >= cutoff[2]
    direction <- as.integer(direction_logical)
    idx <- dmrStat <= cutoff[1]
    direction[idx] <- -1L
    direction[is.na(direction)] <- 0L
    regions <- regionFinder3(x = direction,
                             chr = seqnames_as_char,
                             positions = positions,
                             maxGap = maxGap,
                             verbose = verbose)
    if (is.null(regions$down) && is.null(regions$up)) {
        return(NULL)
    }
    if (verbose) {
        cat("[dmrFinder] creating dmr data.frame\n")
    }
    regions <- do.call(rbind, regions)
    rownames(regions) <- NULL
    regions[["width"]] <- regions[["end"]] - regions[["start"]] + 1L
    regions[["invdensity"]] <- regions[["width"]] / regions[["n"]]
    regions
}

dmrFinder <- function(bstat, cutoff = NULL, qcutoff = c(0.025, 0.975),
                      maxGap = 300, stat = "tstat.corrected", verbose = TRUE) {
    if (is(bstat, "BSseqTstat")) {
        if (!stat %in% colnames(bstat@stats)) {
            stop("'stat' needs to be a column of 'bstat@stats'")
        }
        dmrStat <- bstat@stats[, stat]
    }
    if (is(bstat, "BSseqStat")) {
        dmrStat <- getStats(bstat, what = "stat")
    }
    if (is(dmrStat, "DelayedArray")) {
        # NOTE: Realise in memory since dmrStat is relatively small, being a
        #       column vector with nrow = length(bstat). Subsequent operations
        #       are much simplified by having dmrStat as a regular vector
        dmrStat <- as.vector(dmrStat)
    }
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if (is.null(cutoff)) {
        cutoff <- quantile(dmrStat, qcutoff)
    }
    if (length(cutoff) == 1) {
        cutoff <- c(-cutoff, cutoff)
    }
    seqnames_as_char <- as.character(seqnames(bstat))
    positions <- start(bstat)
    regions <- .dmrFinder(dmrStat = dmrStat,
                          cutoff = cutoff,
                          seqnames_as_char = seqnames_as_char,
                          positions = positions,
                          maxGap = maxGap,
                          verbose = subverbose)
    if (is(bstat, "BSseqTstat")) {
        stats <- getStats(bstat, regions, stat = stat)
        regions <- cbind(regions, stats)
        regions$direction <- ifelse(regions$meanDiff > 0, "hyper", "hypo")
    } else if (is(bstat, "BSseqStat")) {
        stats <- getStats(bstat, regions)
        regions <- cbind(regions, stats)
    }
    regions[order(abs(regions$areaStat), decreasing = TRUE), ]
}

clusterMaker <- function(chr, pos, order.it = TRUE, maxGap = 300){
    nonaIndex <- which(!is.na(chr) & !is.na(pos))
    Indexes <- split(nonaIndex, chr[nonaIndex])
    clusterIDs <- rep(NA, length(chr))
    LAST <- 0
    for (i in seq(along = Indexes)) {
        Index <- Indexes[[i]]
        x <- pos[Index]
        if (order.it) {
            Index <- Index[order(x)]
            x <- pos[Index]
        }
        y <- as.numeric(diff(x) > maxGap)
        z <- cumsum(c(1, y))
        clusterIDs[Index] <- z + LAST
        LAST <- max(z) + LAST
    }
    clusterIDs
}

regionFinder3 <- function(x, chr, positions, keep, maxGap = 300,
                          verbose = TRUE) {
    getSegments3 <- function(x, cid, verbose = TRUE) {
        if (verbose) {
            cat("[regionFinder3] segmenting\n")
        }
        segments <- cumsum(c(1L, diff(x) != 0)) + cumsum(c(1L, diff(cid) != 0))
        names(segments) <- NULL
        if (verbose) {
            cat("[regionFinder3] splitting\n")
        }
        out <- list(up = split(which(x > 0), segments[x > 0]),
                    down = split(which(x < 0), segments[x < 0]),
                    nochange = split(which(x == 0), segments[x == 0]))
        names(out[[1]]) <- NULL
        names(out[[2]]) <- NULL
        names(out[[3]]) <- NULL
        out
    }
    cid0 <- clusterMaker(chr = chr, pos = positions, maxGap = maxGap)
    segments <- getSegments3(x = x, cid = cid0, verbose = verbose)
    out <- vector("list", 2)
    names(out) <- c("up", "down")
    for (ii in 1:2) {
        idx <- segments[[ii]]
        if (length(idx) > 0) {
            out[[ii]] <-
                data.frame(chr = sapply(idx, function(jj) chr[jj[1]]),
                           start = sapply(idx, function(jj) min(positions[jj])),
                           end = sapply(idx, function(jj) max(positions[jj])),
                           idxStart = sapply(idx, function(jj) min(jj)),
                           idxEnd = sapply(idx, function(jj) max(jj)),
                           cluster = sapply(idx, function(jj) cid0[jj[1]]),
                           n = sapply(idx, length),
                           stringsAsFactors = FALSE)
        }
    }
    out
}

