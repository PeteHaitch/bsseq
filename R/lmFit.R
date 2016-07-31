# TODO: Need to resolve potential licensing issue: limma is GPL2 and bsseq is
#       Artistic-2.0, but the code in this file is adapted from R/lmfit.R in
#       limma. If we can't use this code, then we can just fall back to using
#       limma::lmFit() with the attendant loss in speed and increase in memory
#       usage.

#-------------------------------------------------------------------------------
# Adapted from limma::lmFit for special case(s) used in bsseq
# NOTE: y should be transposed, i.e., y should have samples as rows and loci as
#       columns
# Limitations:
#   - `y` must be a matrix
#   - The return value only contains those elements needed in downstream bsseq
#     functionality
#       - `genes` is always set to NULL
#       - `Amean` is not computed or returned
#   - Many options to limma::lmFit() are not supported
#       - `ndups` is not supported
#       - `spacing` is not supported
#       - `block` is not supported
#       - `correlation` is not supported
#       - `weights` is not supported
#       - `method` is not supported

# Benchmarking
# > design
# (Intercept) BSseq$Sexmale
# 1           1             1
# 2           1             1
# 3           1             0
# 4           1             1
# 5           1             0
# 6           1             0
# attr(,"assign")
# [1] 0 1
# attr(,"contrasts")
# attr(,"contrasts")$`BSseq$Sex`
# [1] "contr.treatment"
# > dim(allPs)
# [1] 23059530        6
# > tAllPs <- t(allPs)
# > microbenchmark(lmFit(allPs, design), .lmFit(tAllPs, design), times = 10)
# Unit: seconds
# expr      min       lq     mean   median       uq
# lmFit(allPs, design) 17.35394 19.22649 21.33923 20.52108 24.05139
# .lmFit(tAllPs, design) 11.16108 16.98138 18.31080 19.25468 20.44507
# max neval
# 26.25427    10
# 23.27585    10
# > profvis(lmFit(allPs, design))
# Memory: 0/9412.3 MB'
# > profvis(.lmFit(tAllPs, design))
# Memory: 0/7125.2 MB'

# Summary: .lmFit() saves 7-15% of time and ~25% of memory usage over
#          limma::lmFit()

.lmFit <- function(y, design) {
    stopifnot(is.matrix(y))
    if (is.null(design)) {
        design <- matrix(1, ncol(y), 1)
    } else {
        design <- as.matrix(design)
        if (mode(design) != "numeric") {
            stop("design must be a numeric matrix")
        }
        if (nrow(design) != nrow(y)) {
            stop("row dimension of design doesn't match column dimension of ",
                 "data object")
        }
    }
    ne <- nonEstimable(design)
    if (!is.null(ne)) {
        cat("Coefficients not estimable:", paste(ne, collapse = " "), "\n")
    }

    fit <- .lm.series(M = y, design = design)

    if (NCOL(fit$coef) > 1) {
        n <- rowSums(is.na(fit$coef))
        n <- sum(n > 0 & n < NCOL(fit$coef))
        if (n > 0) {
            warning("Partial NA coefficients for ", n, " probe(s)",
                    call. = FALSE)
        }
    }

    fit$genes <- NULL
    fit$method <- "ls"
    fit$design <- design
    new("MArrayLM", fit)
}

#-------------------------------------------------------------------------------
# Adapted from limma::lm.series for special case(s) used in bsseq
# NOTE: M should be transposed from it's orientation in limma::lm.series(),
#       i.e., M should have samples as rows and loci as columns
# Limitations:
#   - All limitations inherited from .lmFit()
#   - The return value only contains those elements needed in downstream bsseq
#     functionality
#       - `df.residual` is returned as a vector with length 1 rather than with
#         length(df.residual) == ncol(tAllPs). `df.residual` is constant
#         across all loci, so no information is lost. Furthermore, bsseq does
#         not currently use `df.residual` in any downstream computations.
#       - `pivot` is not returned

# Benchmarking
# > microbenchmark(.lm.series(M = tAllPs, design = design), lm.series(M = allPs, design = design), times = 20)
# Unit: seconds
# expr      min       lq     mean
# .lm.series(M = tAllPs, design = design) 10.91780 14.63206 18.99563
# lm.series(M = allPs, design = design) 12.40072 14.78253 19.63373
# median       uq      max neval
# 16.63868 19.39843 38.41019    20
# 18.81319 25.80417 31.37283    20
# > profvis(lm.series(M = allPs, design = design))
# Memory: 0 / 8620.2 MB
# > profvis(.lm.series(M = tAllPs, design = design))
# Memory 0 / 6509.4 MB

# Summary: We can save ~12% running time and ~25% of memory usage by using
#          .lm.series() insetad of limma::lm.series()

.lm.series <- function(M, design = NULL) {
    narrays <- nrow(M)

    # NOTE: Redundant if called from .lmFit()
    if (is.null(design)) {
        design <- matrix(1, narrays, 1)
    } else {
        design <- as.matrix(design)
    }

    nbeta <- ncol(design)
    coef.names <- colnames(design)
    if (is.null(coef.names)) {
        coef.names <- paste("x", 1:nbeta, sep = "")
    }

    ngenes <- ncol(M)
    # NOTE: Initisalise with NA_real_ rather than NA (which is logical) to
    #       ensure the storage mode of these matrices is numeric
    stdev.unscaled <- beta <- matrix(NA_real_, ngenes, nbeta,
                                     dimnames = list(colnames(M), coef.names))

    # Check whether QR-decomposition is constant for all genes
    # If so, fit all genes in one sweep
    NoProbeWts <- all(is.finite(M))
    if (NoProbeWts) {
        # NOTE: If wanting to be really crazy, could use .lm.fit(). This would
        #       require some post-processing of the .lm.fit() output and it is
        #       not obvious that it is worth the effort involved
        fit <- lm.fit(design, M)
        if (fit$df.residual > 0) {
            if (is.matrix(fit$effects)) {
                fit$sigma <- sqrt(colMeans(
                    fit$effects[(fit$rank + 1):narrays, , drop = FALSE]^2))
            } else {
                fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 1):narrays]^2))
            }
        } else {
            fit$sigma <- rep(NA, ngenes)
        }

        # NOTE: In bsseq we only need `sigma`, `cov.coefficients`,
        #       `coefficients`, `stdev.unscaled`
        fit$fitted.values <- fit$residuals <- fit$effects <- NULL
        fit$coefficients <- t(fit$coefficients)
        fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
        est <- fit$qr$pivot[1:fit$qr$rank]
        dimnames(fit$cov.coefficients) <- list(coef.names[est], coef.names[est])
        stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)),
                                        ngenes, fit$qr$rank, byrow = TRUE)
        fit$stdev.unscaled <- stdev.unscaled
        dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
        fit$pivot <- fit$qr$pivot
        return(fit)
    } else {
        stop("'M' contains non-finite values")
    }
}
