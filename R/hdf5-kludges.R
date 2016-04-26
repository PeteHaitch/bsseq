# Various kludges to workaround missing functionality from HDF5Array package.
# Remember, these are *kludges* and are not meant to be the best way to
# implement the desired behaviour.

# NOTE: is.infinite is a Primitive
is.infinite <- function(x) {
    if (is(x, "HDF5Array") || is(x, "DelayedArray")) {
        x <- as.array(x)
    }
    base::is.infinite(x)
}
