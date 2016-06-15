# TODO: Is using a class union the best way to achieve this?
# TODO: Is this the place to put it?
setClassUnion("matrixOrDelayedMatrix", c("matrix", "DelayedMatrix"))

setClass("BSseqTstat", contains = "hasGRanges",
         representation(stats = "matrixOrDelayedMatrix",
                        parameters = "list")
         )
setValidity("BSseqTstat", function(object) {
    msg <- NULL
    if(length(object@gr) != nrow(object@stats))
        msg <- c(msg, "length of 'gr' is different from the number of rows of 'stats'")
    if(is.null(msg)) TRUE else msg
})

setMethod("show", signature(object = "BSseqTstat"),
          function(object) {
              cat("An object of type 'BSseqTstat' with\n")
              cat(" ", length(object), "methylation loci\n")
              cat("based on smoothed data:\n")
              cat(" ", object@parameters$smoothText, "\n")
              cat("with parameters\n")
              cat(" ", object@parameters$tstatText, "\n")
          })

setMethod("[", "BSseqTstat", function(x, i, ...) {
    if(missing(i))
        stop("need [i] for subsetting")
    if(missing(i))
        return(x)
    x@gr <- x@gr[i]
    x@stats <- x@stats[i,, drop = FALSE]
    x
})

BSseqTstat <- function(gr = NULL, stats = NULL, parameters = NULL) {
    out <- new("BSseqTstat")
    out@gr <- gr
    out@stats <- stats
    out@parameters <- parameters
    out
}

summary.BSseqTstat <- function(object, ...) {
    # TODO: quantile,DelayedArray-method to avoid explicit realisation
    tstat_corrected <- getStats(object)[, "tstat.corrected"]
    if (is(tstat_corrected, "DelayedArray")) {
        tstat_corrected <- as.array(tstat_corrected)
    }
    quant <- quantile(tstat_corrected,
                      prob = c(0.0001, 0.001, 0.01, 0.5, 0.99, 0.999, 0.9999))
    quant <- t(t(quant))
    colnames(quant) <- "quantiles"
    out <- list(quantiles = quant)
    class(out) <- "summary.BSseqTstat"
    out
}

print.summary.BSseqTstat <- function(x, ...) {
    print(as.matrix(x$quantiles))
}

plot.BSseqTstat <- function(x, y, ...) {
    # TODO: density,DelayedArray-method to avoid explicit realisation
    tstat <- getStats(x)[, "tstat"]
    if (is(tstat, "DelayedArray")) {
        tstat <- as.array(tstat)
    }
    plot(density(tstat), xlim = c(-10,10), col = "blue", main = "")
    if("tstat.corrected" %in% colnames(getStats(x))) {
        tstat_corrected <- getStats(x)[, "tstat.corrected"]
        if (is(tstat_corrected, "DelayedArray")) {
            tstat_corrected <- as.array(tstat_corrected)
        }
        lines(density(tstat_corrected), col = "black")
        legend("topleft", legend = c("uncorrected", "corrected"), lty = c(1,1),
               col = c("blue", "black"))
    } else {
        legend("topleft", legend = c("uncorrected"), lty = 1,
               col = c("blue"))
    }
}

