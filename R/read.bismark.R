# TODO: Should hdf5 be the default? Should it guess/recommend the user switch
#       to hdf5 if they have large/many files?
read.bismark <- function(files,
                         sampleNames,
                         rmZeroCov = FALSE,
                         strandCollapse = TRUE,
                         fileType = c("cov", "oldBedGraph", "cytosineReport"),
                         mc.cores = 1,
                         verbose = TRUE,
                         hdf5 = FALSE) {
    ## Argument checking
    if (anyDuplicated(files)) {
        stop("duplicate entries in 'files'")
    }
    if (length(sampleNames) != length(files) || anyDuplicated(sampleNames)) {
        stop("argument 'sampleNames' has to have the same length as argument 'files', without duplicate entries")
    }
    fileType <- match.arg(fileType)
    if (verbose) {
        message(paste0("Assuming file type is ", fileType))
    }
    if (strandCollapse && (fileType == "cov" || fileType == "oldBedGraph")) {
        stop("'strandCollapse' must be 'FALSE' if 'fileType' is '", fileType, "'")
    }
    ## SummarizedExperiment validator makes calls to identical() that will fail
    ## due to how sampleNames are propagated sometimes with and without names().
    ## To make life simpler, we simply strip the names()
    sampleNames <- unname(sampleNames)

    ## Process each file
    idxes <- seq_along(files)
    names(idxes) <- sampleNames
    allOut <- mclapply(idxes, function(ii) {
        if (verbose) {
            cat(sprintf("[read.bismark] Reading file '%s' ... ", files[ii]))
        }
        ptime1 <- proc.time()
        if (fileType == "cov" || fileType == "oldBedGraph") {
            out <- read.bismarkCovRaw(thisfile = files[ii],
                                      thisSampleName = sampleNames[ii],
                                      rmZeroCov = rmZeroCov,
                                      hdf5 = hdf5)
        } else if (fileType == "cytosineReport") {
            out <- read.bismarkCytosineReportRaw(thisfile = files[ii],
                                                 thisSampleName = sampleNames[ii],
                                                 rmZeroCov = rmZeroCov,
                                                 keepContext = FALSE,
                                                 hdf5 = hdf5)
        }
        if (strandCollapse) {
            out <- strandCollapse(out)
        }
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if (verbose) {
            cat(sprintf("done in %.1f secs\n", stime))
        }
        out
    }, mc.cores = mc.cores)

    if (verbose) {
        cat(sprintf("[read.bismark] Joining samples ... "))
    }
    ptime1 <- proc.time()
    if (length(allOut) > 1L) {
        allOut <- combineList(allOut, hdf5 = hdf5)
    } else {
        allOut <- allOut[[1L]]
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3L]
    if (verbose) {
        cat(sprintf("done in %.1f secs\n", stime))
    }
    allOut
}

read.bismarkCovRaw <- function(thisfile,
                               thisSampleName,
                               rmZeroCov,
                               hdf5 = FALSE) {

    ## data.table::fread() can't read directly from a gzipped file so, if
    ## necessary, gunzip the file to a temporary location.
    if (isGzipped(thisfile)) {
        thisfile <- gunzip(thisfile,
                           temporary = TRUE,
                           overwrite = TRUE,
                           remove = FALSE)
    }

    ## Read in the file
    ## NOTE: Surprisingly, faster to use read.table() than data.table::fread()
    ##       to check number of columns
    if (ncol(read.table(thisfile, header = FALSE, nrows = 1)) != 6L) {
        stop("File does not appear to be in 'cov' format (ncol != 6)")
    }
    out <- fread(thisfile, header = FALSE, select = c(1, 2, 5, 6))

    ## Create GRanges instance from 'out'
    gr <- GRanges(seqnames = out[[1L]],
                  ranges = IRanges(start = out[[2L]], width = 1L))

    ## Create BSseq instance from 'out'
    BSseq(gr = gr,
          M = as.matrix(out[[5L]]),
          Cov = as.matrix(out[[5L]] + out[[6L]]),
          sampleNames = thisSampleName,
          rmZeroCov = rmZeroCov,
          hdf5 = hdf5)
}

read.bismarkCytosineReportRaw <- function(thisfile,
                                          thisSampleName,
                                          rmZeroCov,
                                          keepContext = FALSE,
                                          hdf5 = FALSE) {

    ## NOTE: keepContext not yet implemented
    if (keepContext) {
        stop("'keepContext' argument not yet implemented")
    }

    ## data.table::fread() can't read directly from a gzipped file so, if
    ## necessary, gunzip the file to a temporary location.
    if (isGzipped(thisfile)) {
        thisfile <- gunzip(thisfile,
                           temporary = TRUE,
                           overwrite = TRUE,
                           remove = FALSE)
    }

    ## Read in the file
    ## NOTE: Surprisingly, faster to use read.table() than data.table::fread()
    ##       to check number of columns
    if (ncol(read.table(thisfile, header = FALSE, nrows = 1)) != 7L) {
        stop("File does not appear to be in 'cytosineReport' format (ncol != 7)")
    }
    out <- fread(thisfile, header = FALSE, select = c(1, 2, 3, 4, 5))

    ## Create GRanges instance from 'out'
    gr <- GRanges(seqnames = out[[1L]],
                  ranges = IRanges(start = out[[2L]], width = 1),
                  strand = out[[3L]])

    ## Create BSseq instance from 'out'
    BSseq(gr = gr,
          sampleNames = thisSampleName,
          M = as.matrix(out[[4L]]),
          Cov = as.matrix(out[[4L]] + out[[5L]]),
          rmZeroCov = rmZeroCov,
          hdf5 = hdf5)

}
