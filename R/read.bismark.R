# TODO: Should hdf5 be the default? Should it guess/recommend the user switch
#       to hdf5 if they have large/many files?
# TODO: Could add a seqinfo argument to read.bismark() that could be used
#       to ensure the output of read.bismark() was sorted. The default
#       value would infer the seqinfo from the bismark input files,
#       although this could not guarantee the final object returned by
#       read.bismark() was sorted
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
    # TODO: If `fileType = "cytosineReport"` and `rmZeroCov = TRUE` then
    #       may not want to remove zero coverage sites until the combined
    #       BSseq object is formed, otherwise lose the opportunity to do the
    #       'cbind() shortcut' in combineList()all
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

    setnames(out, paste0("V", c(1, 2, 5, 6)), c("seqnames", "pos", "M", "U"))
    # NOTE: Sort out in place using data.table wizardy. This is technically a
    #       breaking change (and should be documented) but one that I think is
    #       worth introducing. Sort order is based on sortSeqlevels() and is
    #       limited to individual files, i.e. there is no guarantee that the
    #       object returned by read.bismark() is sorted.
    # NOTE: Drop columns of 'out' once they are no longer required
    seqlevels <- out[, sortSeqlevels(unique(seqnames))]
    out[, "seqnames" := factor(seqnames, seqlevels)]
    setkey(out, seqnames, pos)
    seqnames <- out[, .N, by = seqnames]
    out[, "seqnames" := NULL]
    seqnames <- Rle(seqnames[, seqnames], seqnames[, N])
    strand <- strand(Rle("*", nrow(out)))
    M <- as.matrix(out[, M])
    out[, M := NULL]
    Cov <- M + out[, U]
    out[, U := NULL]
    ## Create GRanges instance from 'out'
    gr <- GRanges(seqnames = seqnames,
                  ranges = IRanges(out[, pos], width = 1L),
                  strand = strand)
    out[, pos := NULL]

    ## Create BSseq instance from 'out'
    BSseq(gr = gr,
          sampleNames = thisSampleName,
          M = M,
          Cov = Cov,
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
    setnames(out, paste0("V", 1:5), c("seqnames", "pos", "strand", "M", "U"))
    # NOTE: Sort out in place using data.table wizardy. This is technically a
    #       breaking change (and should be documented) but one that I think is
    #       worth introducing. Sort order is based on sortSeqlevels() and is
    #       limited to individual files, i.e. there is no guarantee that the
    #       object returned by read.bismark() is sorted.
    # NOTE: Drop columns of 'out' once they are no longer required
    seqlevels <- out[, sortSeqlevels(unique(seqnames))]
    out[, c("seqnames", "strand") :=
            list(factor(seqnames, seqlevels), factor(strand, levels(strand())))]
    setkey(out, seqnames, strand, pos)
    seqnames_strand <- out[, .N, by = list(seqnames, strand)]
    out[, c("seqnames", "strand") := NULL]
    seqnames <- Rle(seqnames_strand[, seqnames], seqnames_strand[, N])
    strand <- Rle(seqnames_strand[, strand], seqnames_strand[, N])
    M <- as.matrix(out[, M])
    out[, M := NULL]
    Cov <- M + out[, U]
    out[, U := NULL]
    ## Create GRanges instance from 'out'
    gr <- GRanges(seqnames = seqnames,
                  ranges = IRanges(out[, pos], width = 1L),
                  strand = strand)
    out[, pos := NULL]

    ## Create BSseq instance from 'out'
    BSseq(gr = gr,
          sampleNames = thisSampleName,
          M = M,
          Cov = Cov,
          rmZeroCov = rmZeroCov,
          hdf5 = hdf5)
}
# NOTE: Required to quieten R CMD check
globalVariables(c("N", "U"))
