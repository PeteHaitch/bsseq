# HDF5Array notes

Bits and pieces noticed during refactor of _bsseq_ to work with 
_HDF5Array_-based assays.

## Things to be aware of

`x` is a _HDF5Array_ object

- Many (all?) operations convert `x` to a _DelayedMatrix_ object
    - E.g., `[,HDF5Array,*-method`, `dimnames<-,HDF5Array,*-method`
    - Use `as.array()` to realise the result in memory
    - Use `HDF5Dataset()` to realise the result on disk
- The location of automatically created HDF5 datasets is set by 
`setHDF5DumpFile()`, which uses `tempfile()`
    - Will want to use fast, local disk for intermediate files but shared 
    (possibly slower) disk for more permanent files
- Try to wrap all operations in a _DelayedArray_, then realise in memory (using 
`as.array()`) or on disk (using `writeHDF5Dataset()`) at the last possible 
moment. There may be occasions where its better to realise a result in memory 
and then wrap it in a _DelayedArray_. E.g., if all seeds are themselves 
_DelayedArray_ objects and are in memory, then it might make sense to realise 
the result in memory and wrap the result in a _DelayedArray_ rather than to 
store all the seeds. It might be hard to know when to do this automatically, 
but I think it should be an easy option for a developer/user to invoke. Here's 
an illustration
- You can a actually have a matrix, `x`, where 
`length(x) > .Machine$integer.max`; the real limit is that 
`nrow(x) <= .Machine$integer.max && ncol(x) <= .Machine$integer.max` must be 
`TRUE`. Nonetheless, such a matrix is obviously going to be **huge** in memory.

```r
x <- matrix(1:1000000)
a <- DelayedArray(x)
b <- DelayedArray(x * 10)
d <- DelayedArray(x * 100)
e1 <- (a + b) / d
e2 <- as.array(e1)
e3 <- DelayedArray(e2)
pryr::object_size(e1)
pryr::object_size(e2)
pryr::object_size(e3)
```

## Questions

### Asked Hervé 2016-04-24

- [x] Should a _DelayedMatrix_ be realised as a _HDF5Array_ (or _array_) before including as an assay element in a _SummarizedExperiment_ object?

Best to realise as _HDF5Array_ upon construction of the BSseq object, i.e., as 
a final step of the `BSseq()` constructor.

From Hervé:

> You can stick a DelayedArray object directly as an assay in a
SummarizedExperiment object. Or you can realize it before by calling
HDF5Array() on it (equivalent to calling HDF5Array(HDF5Dataset())
on it). Both work. But because realization is expensive you would
normally do the former. The general advice is to delay realization
as much as possible.

> Note that, in theory, a DelayedArray object that carries a lot of 
delayed operations on it will be a little bit slower when it comes
to realizing it than a DelayedArray object with few delayed operations
on it. But this is only in theory. It could be that in practice you
won't observe any significant difference.

> When you call HDF5Array() on it, you realize all the delayed operations
to an HDF5 dataset on disk, and get back a DelayedArray object that
points to the new dataset and doesn't carry any delayed operation (i.e.
a pristine DelayedArray object). Note that an HDF5Array object is just
a pristine DelayedArray object whose seed is an HDF5Dataset object.
In other words, you're starting fresh with a new pristine DelayedArray
object. Because the object is pristine, it's in sync with the data on
disk. As soon as you start performing delayed operations on it, it goes
out of sync with the data on disk.

> This new pristine object is in theory a little bit smaller than before
realization (because the stack of delayed operation is empty), and will
also be slightly faster to realize next time, but that counts very
little in the balance compared to the high price you just paid to make
it pristine.

### Asked Hervé 2016-06-16

- [x] `pmax2()`, `pmin2()` exist but not `pmax()` and `pmin()`; why?
> There are several problems with `base::pmax()` and `base::pmin()` that make 
them hard to reuse here.

- [x] If combining (e.g., with `cbind()`/`rbind()` or as part of 
`bsseq::combine()`/ `bsseq::combineList()`) _DelayedArray_ objects with a 
_HDF5Dataset_ `@seed`, how to decide if the method should write a new 
HDF5 file or return a _DelayedArray_ object? E.g., is there a substantial 
overhead for parsing _DelayedArray_ objects comprising multiple HDF5 files in 
subsequent calls?
> The answer is the same as for your "should a DelayedArray be realized
as a HDF5Dataset before going in a SummarizedExperiment object"
question. The general advice is to always delay realization as much
as possible. That's because the price of realization is much much
bigger than the overhead of operating on a DelayedArray with multiple
seeds.

- [x] Relatedly, might there be need for a method to check if a 
_DelayedArray_ has a _HDF5Dataset_ for its `seed`? E.g., might want to 
`cbind()` a bunch of _DelayedArray_ objects and write to disk as new 
_HDF5Array_ iff any of the objects where themselves HDF5-backed. This has 
implications for the code where I do `is(x, "DelayedArray")` because I am 
currently assuming this is equivalent to "is `x` a HDF5-backed _DelayedArray_".
> I'm not sure we need this. But if you show me a concrete use case where
this might be useful I could change my mind.

This is similar to the issue I'm having when an operation makes a new 
_DelayedArray_ with multiple seed (see 'Things to be aware of' above).

- [x] If a package offers HDF5Array and array/matrix options, how to best 
specify and/or determine that function/method should return HDF5-backed result. 
Currently, using `hdf5` argument to relevant functions. Perhaps there should be 
a package option? Perhaps function/method should default to `hdf5 = TRUE` if 
input includes HDF5-backed data?
> A good question. For which I'm not sure I have a good answer. I've
plans to modify summarizeOverlaps() to let the user choose if they
want the returned SummarizedExperiment object to be in memory (the
default) or HDF5Array-backed. For this one I think I'll add an
extra argument to the function.

### TODO: Ask Hervé

- [ ] How to realise DelayedArray object in memory but wrap result in a 
DelayedArray? E.g.,

```r
x <- matrix(1:1000000)
a <- DelayedArray(x)
b <- DelayedArray(x * 10)
d <- DelayedArray(x * 100)
e1 <- (a + b) / d
e2 <- as.array(e1)
e3 <- DelayedArray(e2)
# Better of realising in memory and then wrapping in DelayedArray (size e3 < e1)
pryr::object_size(e1)
pryr::object_size(e2)
pryr::object_size(e3)
```

- [ ] Why does `rowMeans,DelayedMatrix` not use `base::rowMeans()` internally 
but rather uses '`base::rowSums() / ncol()`'?
- [ ] How to make function, e.g., `matrixStats::rowVars()`, 'block processing 
friendly'?
- [ ] What should the behaviour of `cbind()`, `rbind()`, `combine()`, 
`arbind()`, `acbind()` be if some of inputs are array objects and some are 
DelayedArray objects? E.g., when combining SummarizedExperiment objects where 
some have their assays in memory and others have their assays on disk.

## Missing methods

Need methods implements for _HDF5Array_ and probably _DelayedArray_ classes:

### Asked Hervé 2016-04-24

- [x] `is.infinite()`
- [x] `split()`
- [x] `arbind()`, `acbind()`
- [x] `which()`

Added as of `v1.1.8`.

## Bugs

### Asked Hervé 2016-04-24

- [x] `do.call(cbind, list(x, x))` works but `do.call(cbind, list(a = x, b = x))` errors, where `x` is a _HDF5Matrix_ or a _DelayedMatrix_

Fixed as of `v1.1.8`.

### Asked Hervé 2016-06-16

- [x] `rbind,DelayedMatrix-method` and `cbind,DelayedMatrix-method` fail if 
`...` includes `NULL` objects. Contrast with `rbind()` and `cbind()` if `...` 
includes _array_ and `NULL` objects
> Unfortunately I'm not sure I can do anything about this. The `NULL`
arguments break dispatch ...

So I need to remove `NULL` elements before binding, e.g., 
`do.call(rbind, (lst[sapply(lst, function(x) !is.null(x))]))`

### TODO: Ask Hervé

- [ ] `X1 <- HDF5Array(x)` ensures `dimnames(X1)` are identical to `dimnames(x)`, 
      however, `X2 <- HDF5Array(writeHDF5Dataset(x))` does not (`dimnames(X2)` 
      are `NULL`). Can this be fixed?
    

## Wishlist

### Asked Hervé 2016-04-24

- [x] An efficient helper function to create a _HDF5Array_ from a _data.frame_ with numeric columns (currently have to go via `as.matrix()`, which incurs a copy, I think)

Added as of `v1.1.8`.

### Asked Hervé 2016-06-16

- [x] Could `rowSums()` (resp. `colSums`) be optimised if `@index` has a
**large** amount of i-subsetting (resp. j-subsetting). E.g.:

```r
library(HDF5Array)
h <- HDF5Array(matrix(1:200000, ncol = 8))
i <- sample(nrow(h), nrow(h) / 2)
ii <- sort(i)
system.time(rowSums(h))
system.time(rowSums(h[i, ]))
system.time(rowSums(h)[i])
system.time(rowSums(h[ii, ]))
system.time(rowSums(h)[ii])
```
> This is caused by an inefficiency in rhdf5::h5read() when 'index' is 
specified ... Surprisingly, reading a subset of the dataset is slower than 
reading the full dataset!

HDF5ARray may get its own `h5read()`-like function, optimised for its specific 
use case(s)

- [x] Is 1D-style subsetting ever going to be supported? E.g. this is an op 
used a bit in __bsseq__:

```r
library(HDF5Array)
m <- matrix(1:10, ncol = 2)
M <- HDF5Array(m)
m[which(m > 3)]
# Errors
M[which(M > 3)]
```
> Sure. I'll add this.

- [x] `as.integer,DelayedArray-method` for when _DelayedArray_ contains logical
> We can't have as.integer() do this because as.integer() should always
return an integer vector (i.e. is.integer() must be TRUE on the
result). This is a very strong expectation. Breaking this rule is too
dangerous e.g. for code that calls as.integer() on an object before
passing it to a C function.

A really neat alternative is to do the following (a delayed op)

```r
x <- DelayedArray(matrix(sample(c(TRUE, FALSE), 1000, replace = TRUE)))
x + 0L
```

- [x] Is it possible to serialise a _HDF5Array_? `saveRDS()` will store 
original `@seed@file` which is not guaranteed to exist across sessions.

This proved much more complicated than I first thought. There are 
two different scenarios to consider:

1. "Going home for the day". I want to save my SE and load it again tomorrow
to continue work on the same system.
2. "Sharing with others". I want to share with collaborators/colleagues/BioC
community a file that contains everything they need to create the
SE on their machine.

For (1), it makes sense to set a non-temporary, absolute path for the .h5 files
using setHDF5DumpFile(). Then the serialised SE is small on disk and the .h5
files will still exist at the end of the R session, so tomorrow when we load the
SE everything "just works" [there's probably some complications to worry about
if we use a relative path in setHDF5DumpFile() and load the object in
a different
working dir]. In this case, there's no point in carrying additional
copies of the
(possibly large) .h5 file(s) in the serialised object written to disk.

For (2), perhaps its as simple as bundling up the saved SE object (.rda or .rds)
with the .h5 files in a tarball (or similar).

See email exchange for full details and brainstorming. This is going to require 
a fair bit of work to implement in full generality.

### TODO: Ask Hervé

- [ ] `DelayedArray,vector` so I don't have to keep doing 
`DelayedArray(as.matrix(x))`
- [ ] `setHDF5DumpDir()` as well as `setHDF5DumpFile()`

## Files checked

- [x] `BSmooth.fstat.R`
- [x] `BSmooth.R`
- [x] `BSmooth.tstat.R`
- [x] `BSseq_class.R`
- [x] `BSseq_utils.R`
- [x] `BSseqStat_class.R`
- [x] `BSseqTstat_class.R`
- [x] `combine.R`
- [x] `dmrFinder.R`
- [x] `fisher.R`
- [x] `getStats.R`
- [x] `gof_stats.R` 
- [x] `hasGRanges.R`
- [x] `permutations.R`
- [x] `plotting.R`
- [x] `read.bismark.R`
- [x] `read.bsmooth.R`
- [x] `utils.R`

__All files checked for compatability issues with HDF5Array-backed BSseq objects.__ 
Many patches were added and in some places the code much simplified. However, 
there remain a number of __TODO__s, some of these are code optimisations, 
others are design decisions I need to discuss with Kasper. I also really want 
to add more unit tests to BSseq.

## Misc. TODOs

- [ ] `is()` vs. `inherits()` (motivated by [https://github.com/Bioconductor-mirror/BiocParallel/commit/420aeff4a222415908a4fd9028d907c473f42043](https://github.com/Bioconductor-mirror/BiocParallel/commit/420aeff4a222415908a4fd9028d907c473f42043))
    - [ ] Relatedly, should I be checking for _HDF5[Array|Matrix]_ or 
         _DelayedArray[Array|Matrix]_?
- [ ] Default value of hdf5? Currently FALSE to preserve previous behaviour, 
      but an 'auto' option could be useful, e.g., if input uses HDF5 then so 
      should output or if output is going to be big then use HDF5 (e.g., 
      output of read.bismark())
- [ ] Should bsseq DEPEND on HDF5Array so that HDF5Array is automatically added 
      to the search path?
    - [ ] Relatedly, `library(HDF5Array)` or `require(HDF5Array)` or something 
          else when using HDF5Array in examples.
- [ ] FWER of DMRs sounds like it could be candidate for optimisation
- [ ] Parallel writes to a `.h5` file are not possible

```r
# This errors on a disk with slow access and seems to create a dodgy .h5 file 
# on a fast disk
mclapply(1:4, function(x) writeHDF5Dataset(DelayedArray(matrix(1:1000000)), file = getHDF5DumpFile(), name = paste0("kraken_", x)), mc.cores = 4)
```
    - [ ] Instead, can write each array to its own h5 file. Perhaps even 
          better, write each samples `M`, `Cov`, `coef`, and `se.coef` - will 
          need to be careful with writing `coef` and `se.coef`
            
```
# What about if we write to a new .h5 file for each array
mclapply(1:4, function(x) {
    writeHDF5Dataset(DelayedArray(matrix(1:1000000)), 
                     file = tempfile("matrix", tmpdir = ".", fileext = ".h5"),
                     name = "M")
}, mc.cores = 4)
```
- [ ] Use of `HDF5Array()` vs. `writeDataHDF5Dataset()`; `HDF5Array()` (via 
      `HDF5Dataset`) errors if `file` is a DelayedArray unless `name = NA` and 
      `type` is missing. I think this means I want to use 
      `writeDataHDF5Dataset()` whenever it needs a `file` (other than 
      `file = getHDF5DumpFile()` and/or when it needs a `name` (other than 
      `name = getHDF5DumpName()`), e.g., whenever writes may be occuring in 
      parallel.
      - [ ] Will need to wrap the result in a call to `HDF5Array()` whenever 
            it is to be used in downstream stuff.
