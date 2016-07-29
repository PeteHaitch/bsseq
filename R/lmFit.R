#-------------------------------------------------------------------------------
# Adapted from limma::lmFit for special case(s) used in bsseq
# NOTE: y should be transposed, i.e., samples as rows and loci as columns
# Limitations:
# - `y` must be a matrix
# - Many options to limma::lmFit() are not supported
#   - `ndups` is not supported
#   - `spacing` is not supported
#   - `block` is not supported
#   - `correlation` is not supported
#   - `weights` is not supported
#   - `method` is not supported

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
# > tallPs <- t(allPs)

# > microbenchmark(lmFit(allPs, design), .lmFit(tallPs, design), times = 10)
# Unit: seconds
# expr      min       lq     mean   median       uq
# lmFit(allPs, design) 14.00243 14.90291 20.45310 20.19567 22.85180
# .lmFit(tallPs, design) 11.53216 12.86627 17.97016 16.47999 17.81626
# max neval
# 30.41207    10
# 36.08204    10
# profvis(lmFit(allPs, design)) reports 'Memory: 0/7916.9 MB'
# profvis(.lmFit(allPs, design)) reports 'Memory 0/6509.4 MB'

# Summary: .lmFit() saves ~4 seconds (~25% of time) and ~1.4 GB (~18% of
#          memory) over lmFit()
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
            stop("row dimension of design doesn't match column dimension of data object")
        }
    }
    ne <- limma::nonEstimable(design)
    if (!is.null(ne)) {
        cat("Coefficients not estimable:", paste(ne, collapse = " "), "\n")
    }

    fit <- .lm.series(M = y,
                      design = design)

    if (NCOL(fit$coef) > 1) {
        n <- rowSums(is.na(fit$coef))
        n <- sum(n > 0 & n < NCOL(fit$coef))
        if (n > 0) {
            warning("Partial NA coefficients for ", n, " probe(s)",
                    call. = FALSE)
        }
    }

    # NOTE: probes are always set to NULL in .lmFit()
    fit$genes <- NULL
    fit$method <- "ls"
    fit$design <- design
    new("MArrayLM", fit)
}

#-------------------------------------------------------------------------------
# Adapted from limma::lm.series for special case(s) used in bsseq
# NOTE: M should be transposed from it's orientation in limma::lm.series(),
#       i.e., samples as rows and loci as columns
# Limitations:
# - All limitations inherited from .lmFit()

# Benchmarking
# > microbenchmark(.lm.series(M = tallPs, design = design), lm.series(M = allPs, design = design), times = 20)
# Unit: seconds
# expr       min       lq     mean
# .lm.series(M = tallPs, design = design)  9.850865 13.20839 17.84259
# lm.series(M = allPs, design = design) 11.403523 15.95992 19.87436
# median       uq      max neval
# 16.80557 20.86947 36.10210    20
# 19.32464 21.44028 34.83261    20

# Summary: We can save ~2-3 seconds by using .lm.series(). It's important to
#          note there was considerable run-to-run variation, but .lm.series()
#          always came out ahead of lm.series()
.lm.series <- function(M, design = NULL) {
    narrays <- nrow(M)

    # NOTE: Redundanct if call from .lmFit()
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
        # TODO: If wanting to be really ballsy, could use .lm.fit() but would
        #       require some post-processing of the .lm.fit() output
        fit <- lm.fit(design, M)
        if (fit$df.residual > 0) {
            if (is.matrix(fit$effects)) {
                fit$sigma <- sqrt(colMeans(
                    fit$effects[(fit$rank + 1):narrays, , drop = FALSE] ^ 2))
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
        # TODO: Do we need df.residual to be repeated in order to use with
        #       classifyTestsF()? If not, just record once
        fit$df.residual <- rep.int(fit$df.residual, ngenes)
        dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
        fit$pivot <- fit$qr$pivot
        return(fit)
    } else {
        stop("'M' contains non-finite values")
    }
}
