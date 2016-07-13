# NOTE: All elements of the list in the stats slot are now matrix-like (i.e.
#       2-dimensional) instead of some being vector-like (i.e. 1-dimensional)
# TODO: Perhaps the stats slot should be a matrixOrDelayedMatrix? The main
#       complication is that the number of columns now depends on the number of
#       coefficients whereas for the BSseqTstat the number of columns is fixed.
#       I don't know that this buys us anything except consistency with
#       BSseqTstat (and the change will require some code re-writing).
setClass("BSseqStat", contains = "hasGRanges",
         representation(stats = "list",
                        parameters = "list")
         )

setValidity("BSseqStat", function(object) {
    msg <- NULL
    if(is.null(names(object@stats)) || any(names(object@stats) == "") ||
       anyDuplicated(names(object@stats)))
        msg <- validMsg(msg, "the 'stats' list needs to be named with unique names.")
    for(name in c("rawTstats", "rawSds", "smoothSds", "stat")) {
        if(name %in% names(object@stats) &&
           (!is.matrix(object@stats[[name]])) &&
           !is(object@stats[[name]], "DelayedArray"))
            msg <- validMsg(msg, sprintf("component '%s' of slot 'stats' have to be a matrix or a DelayedArray", name))
        if(name %in% names(object@stats) && nrow(object@stats[[name]]) != length(object@gr))
            msg <- validMsg(msg, sprintf("component '%s' of slot 'stats' have to have the same number of rows as slot 'gr' is long", name))
    }
    if(is.null(msg)) TRUE else msg
})

setMethod("show", signature(object = "BSseqStat"),
          function(object) {
              cat("An object of type 'BSseqStat' with\n")
              cat(" ", length(object), "methylation loci\n")
              cat("based on smoothed data:\n")
              cat(" ", object@parameters$smoothText, "\n")
              cat("with parameters\n")
              cat(" ", object@parameters$tstatText, "\n")
          })

setMethod("[", "BSseqStat", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    statnames <- names(x@stats)
    names(statnames) <- statnames
    x@stats <- lapply(statnames, function(nam) {
        if(nam %in% c("rawSds", "rawTstats", "modelCoefficients", "smoothSds", "stat")) {
            stopifnot(is.matrix(x@stats[[nam]]) ||
                          is(x@stats[[nam]], "DelayedArray"))
            return(x@stats[[nam]][i, , drop = FALSE])
        }
        x@stats[[nam]]
    })
    x
})

BSseqStat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    out <- new("BSseqStat")
    out@gr <- gr
    out@stats <- stats
    out@parameters <- parameters
    out
}


## summary.BSseqStat <- function(object, ...) {
##     quant <- quantile(getStats(object)[, "tstat.corrected"],
##                       prob = c(0.0001, 0.001, 0.01, 0.5, 0.99, 0.999, 0.9999))
##     quant <- t(t(quant))
##     colnames(quant) <- "quantiles"
##     out <- list(quantiles = quant)
##     class(out) <- "summary.BSseqStat"
##     out
## }

## print.summary.BSseqStat <- function(x, ...) {
##     print(as.matrix(x$quantiles))
## }

## plot.BSseqStat <- function(x, y, ...) {
##     tstat <- getStats(x)[, "tstat"]
##     plot(density(tstat), xlim = c(-10,10), col = "blue", main = "")
##     if("tstat.corrected" %in% colnames(getStats(x))) {
##         tstat.cor <- getStats(x)[, "tstat.corrected"]
##         lines(density(tstat.cor), col = "black")
##         legend("topleft", legend = c("uncorrected", "corrected"), lty = c(1,1),
##                col = c("blue", "black"))
##     } else {
##         legend("topleft", legend = c("uncorrected"), lty = 1,
##                col = c("blue"))
##     }
## }

