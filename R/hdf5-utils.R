# TODO: Some of these probably belong in HDF5Array package

setHDF5DumpDir <- function(dir = tempdir()) {
    # TODO: A bunch of argument checks
    assign("dir", dir, envir = HDF5Array:::.HDF5_dump_settings_envir)
}

# NOTE: Currently errors if no directory has been set. Could make tempdir() the
#       default, but at this stage want to encourage/force the user to
#       think about where this should be to avoid surprises down the line
#       (e.g., session ending and losing .h5 files).
getHDF5DumpDir <- function() {
    z <- try(get("dir", envir = HDF5Array:::.HDF5_dump_settings_envir),
             silent = TRUE)
    if (inherits(z, "try-error")) {
        stop("Must set HDF5 dump directory before creating HDF5-backed ",
             "objects see '?setHDF5DumpDir'")
    }
    z
}

.newHDF5Filename <- function(pattern = "file", dir = getHDF5DumpDir()) {
    tempfile(pattern = pattern, tmpdir = dir, ".h5")
}

# A "parallel-safe" way to write an array or DelayedArray to disk as a '.h5'
# file. By that I mean unlike HDF5Array() which always writes to the same file,
# .safeHDF5Array() always writes to a new file so multiple calls in parallel to
# .safeHDF5Array() can't clobber each other by writing to the same file.
# TODO (longterm): Is there some fancy way to infer what package called
#      safeHDF5Array() and set that as `pattern`?
.safeHDF5Array <- function(x, pattern = "file", name = basename(tempfile(""))) {
    stopifnot(is.array(x) || is(x, "DelayedArray"))
    file <- .newHDF5Filename(pattern = pattern)
    HDF5Array(writeHDF5Dataset(x = x, file = file, name = name))
}
