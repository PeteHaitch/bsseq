\name{BSmooth.fstat}
\alias{BSmooth.fstat}
\title{
    Compute F-statistics based on smoothed whole-genome bisulfite
    sequencing data.
}
\description{
    Compute F-statistics based on smoothed whole-genome bisulfite
    sequencing data.
}
\usage{
BSmooth.fstat(BSseq, design, contrasts, verbose = TRUE)
}
\arguments{
    \item{BSseq}{An object of class \code{BSseq}.}
    \item{design}{The design matrix of the bisulfite-sequencing experiment,
        with rows corresponding to arrays and columns to coefficients to be
        estimated.}
    \item{contrasts}{Numeric matrix with rows corresponding to columns in
        \code{design} and columns containing contrasts. May be a vector if
        there is only one contrast.}
    \item{verbose}{Should the function be verbose?}
}
\details{
    \strong{TODO}
}
\value{
    An object of class \linkS4class{BSseqStat}.
}
\author{
    Kasper Daniel Hansen \email{khansen@jhsph.edu}
}
\seealso{
    \code{\link{BSmooth}} for the input object and
    \linkS4class{BSseq} for its class.
    \linkS4class{BSseqTstat} describes the return class. This
    function is likely to be followed by the use of
    \code{\link{dmrFinder}}.
}
\examples{
  \dontrun{
    ## TODO: Kind of a forced example; use data from a multifactorial experiment
    if(require(bsseqData)) {
        data(keepLoci.ex)
        data(BS.cancer.ex.fit)
        BS.cancer.ex.fit <- updateObject(BS.cancer.ex.fit)
        design <- model.matrix(~0 + BS.cancer.ex.fit$Type)
        colnames(design) <- gsub("BS\\\.cancer\\\.ex\\\.fit\\\$Type", "",
                                 colnames(design))
        contrasts <- makeContrasts(
            cancer_vs_normal = cancer - normal,
            levels = design
        )
        ## Remember to subset the BSseq object, see vignette for explanation
        BS.stat <- BSmooth.fstat(BS.cancer.ex.fit[keepLoci.ex,],
                                 design,
                                 contrasts)
        BS.stat

        # Using HDF5Array
        library(HDF5Array)
        # TODO: Subset before calling HDF5Array() otherwise the resulting
        #       DelayedMatrix carries around the i-subsetting in the `@index`
        #       slot, which must be realised when any delayed ops are
        #       performed; contacting Hervé to fix
        BS.cancer.ex.fit.hdf5 <- BS.cancer.ex.fit[keepLoci.ex, ]
        assays(BS.cancer.ex.fit.hdf5) <- endoapply(assays(BS.cancer.ex.fit.hdf5),
                                                   HDF5Array)
        BS.stat.hdf5 <- BSmooth.fstat(BS.cancer.ex.fit.hdf5, design, contrasts,
                                      hdf5 = TRUE)
        BS.stat.hdf5
    }
  }
}
\keyword{internal}
