# HDF5Array notes

Bits and pieces noticed during refactor of _bsseq_ to work with _HDF5Array_-based assays.

## Things to be aware of

`x` is a _HDF5Array_ object

- Many operations convert `x` to a _DelayedMatrix_ object
    - E.g., `[,HDF5Array,*-method`, `dimnames<-,HDF5Array,*-method`
    - Use `as.array()` to realise the result in memory
    - Use `HDF5Dataset()` to realise the result on disk
- The location of automatically created HDF5 datasets is set by `setHDF5DumpFile()`, which uses `tempfile()`
    - Will want to use fast, local disk for intermediate files but shared (possibly slower) disk for more permanent files
- _HDF5Array_ do not support character arrays
    - ~`HDF5Array(matrix(letters[1:10], ncol = 2))` returns an error because it things it's being passed a filename.~ 
    - No longer true
    
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

### TODO: Ask Hervé

- [ ] `pmax2()`, `pmin2()` exist but not `pmax()` and `pmin()`; why?
- [ ] `rowSums()` is painfully slow on test data `bsseqData::BS.cancer.ex.fit` 
with ~600k rows
- [ ] If combining _DelayedArray_ objects (e.g., with `cbind()`/`rbind()` or as 
part of `bsseq::combine()`/ `bsseq::combineList()`), should the method write a 
new HDF5 file or return a _DelayedArray_ object? E.g., is there a substantial 
overhead for parsing multiple HDF5 files in subsequent calls?
- [ ] How to check if a _DelayedArray_ has a `HDF5Dataset` for its `seed`? E.g., 
might want to `cbind()` a bunch of _DelayedArray_ objects and write to disk as 
new _HDF5Array_ iff any of the objects where themselves HDF5-backed. This has 
implications for the code where I do `is(x, "DelayedArray")` because I am 
currently assuming this is equivalent to "is `x` a HDF5-backed _DelayedArray_".
- [ ] If package offers HDF5Array and array/matrix options, how to best specify 
and/or determine that function/method should return HDF5-backed result. 
Currently, using `hdf5` argument to relevant functions. Perhaps there should be 
a package option? Perhaps function/method should default to `hdf5 = TRUE` if 
input includes HDF5-backed data?

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

## Wishlist

### Asked Hervé 2016-04-24

- [x] An efficient helper function to create a _HDF5Array_ from a _data.frame_ with numeric columns (currently have to go via `as.matrix()`, which incurs a copy, I think)

Added as of `v1.1.8`.

### TODO: Ask Hervé

- [ ] Suggest that `rowSums()` (resp. `colSums`) be optimised if `@index` has a **large** amount of i-subsetting (resp. j-subsetting). E.g.:

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

- [ ] Is 1D-style subsetting ever going to be supported? E.g.:

```r
library(HDF5Array)
m <- matrix(1:10, ncol = 2)
M <- HDF5Array(m)
m[which(m > 3)]
M[which(M > 3)]
```
- [ ] How to add additional DelayedArray-based methods? E.g.,
    - [ ] `quantile()`
    - [ ] `density()`

- [ ] `as.integer,DelayedArray-method` for when _DelayedArray_ contains logical
- [ ] Is it possible to serialise a _HDF5Array_? `saveRDS()` will store 
original `@seed@file` which is not guaranteed to exist across sessions.

## Files checked

- [ ] `BSmooth.fstat.R`
- [x] `BSmooth.R`
- [x] `BSmooth.tstat.R`
- [x] `BSseq_class.R`
- [x] `BSseq_utils.R`
- [ ] `BSseqStat_class.R`
- [x] `BSseqTstat_class.R`
- [x] `combine.R`
- [ ] `dmrFinder.R` (need to first update `BSseqStat_class.R`)
- [x] `fisher.R`
- [ ] `getStats.R` (need to first update `BSseqStat_class.R`)
- [x] `gof_stats.R`
- [x] `hasGRanges.R`
- [ ] `permutations.R`
- [ ] `plotting.R`
- [x] `read.bismark.R`
- [ ] `read.bsmooth.R`
- [ ] `utils.R`
