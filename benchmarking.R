# BENCHMARKING
# -------------------------------------------------------------------------
# NOTE: Faster to extract stats and then split rather than split and
#       then extract, especially when BSseqTstat@stats is a DelayedArray
stats_hdf5 <- BSseqTstat@stats
stats_array <- as.array(stats_hdf5)

extract_then_split_array <- function() {
    s <- stats_array[queryHits(ov), stat]
    base::tapply(s, subjectHits(ov), function(x) {
        c(sum(x), max(x))
    }, simplify = FALSE)
}
extract_then_split_hdf5 <- function() {
    s <- stats_hdf5[queryHits(ov), stat]
    tapply(s, subjectHits(ov), function(x) {
        c(sum(x), max(x))
    }, simplify = FALSE)
}
extract_then_split_array_old <- function() {
    lapply(split(stats_array[queryHits(ov), stat], subjectHits(ov)),
           function(x) {
               areaStat <- sum(x)
               maxStat <- max(x)
               c(areaStat, maxStat)
           })
}
extract_then_split_hdf5_realised <- function() {
    s <- as.array(stats_hdf5[queryHits(ov), stat])
    tapply(s, subjectHits(ov), function(x) {
        c(sum(x), max(x))
    }, simplify = FALSE)
}

extract_then_split_hdf5_old <- function() {
    lapply(split(stats_hdf5[queryHits(ov), stat], subjectHits(ov)),
           function(x) {
               areaStat <- sum(x)
               maxStat <- max(x)
               c(areaStat, maxStat)
           })
}
all.equal(extract_then_split_array(), as.list(extract_then_split_hdf5()))
all.equal(extract_then_split_array(), extract_then_split_hdf5_realised())
all.equal(extract_then_split_array(), extract_then_split_array_old())
all.equal(extract_then_split_array(), as.list(extract_then_split_hdf5_old()))


split_then_extract_array <- function() {
    lapply(split(queryHits(ov), subjectHits(ov)), function(idx, stats) {
        mat <- stats[idx,, drop=FALSE]
        areaStat <- sum(mat[, stat])
        maxStat <- max(mat[, stat])
        c(areaStat, maxStat)
    }, stats = stats_array)
}
split_then_extract_hdf5 <- function() {
    lapply(split(queryHits(ov), subjectHits(ov)), function(idx, stats) {
        mat <- stats[idx,, drop=FALSE]
        areaStat <- sum(mat[, stat])
        maxStat <- max(mat[, stat])
        c(areaStat, maxStat)
    }, stats = stats_hdf5)
}
all.equal(extract_then_split_array(), split_then_extract_array())
all.equal(extract_then_split_array(), as.list(split_then_extract_hdf5()))

# NOTE: Not testing split_then_extract_hdf5() because it's so slow
microbenchmark::microbenchmark(extract_then_split_array(),
                               extract_then_split_array_old(),
                               extract_then_split_hdf5(),
                               extract_then_split_hdf5_old(),
                               extract_then_split_hdf5_realised(),
                               split_then_extract_array(),
                               times = 10)
# -------------------------------------------------------------------------
