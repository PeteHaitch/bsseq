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
    - `HDF5Array(matrix(letters[1:10], ncol = 2))` returns an error because it things it's being passed a filename
    
## Questions

### Asked Hervé 2016-04-24

- [x] Should a _DelayedMatrix_ be realised as a _HDF5Array_ (or _array_) before including as an assay element in a _SummarizedExperiment_ object?

## Missing methods

Need methods implements for _HDF5Array_ and probably _DelayedArray_ classes:

### Asked Hervé 2016-04-24

- [ ] `is.infinite()`
- [ ] `split()`
- [ ] `arbind()`, `acbind()`
- [ ] `which()`

## Bugs

### Asked Hervé 2016-04-24

- [ ] `do.call(cbind, list(x, x))` works but `do.call(cbind, list(a = x, b = x))` errors, where `x` is a _HDF5Matrix_ or a _DelayedMatrix_

## Wishlist

### Asked Hervé 2016-04-24

- An efficient helper function to create a _HDF5Array_ from a _data.frame_ with numeric columns (currently have to go via `as.matrix()`, which incurs a copy, I think)

## Files checked

- [ ] `BSmooth.fstat.R`
- [x] `BSmooth.R`
- [ ] `BSmooth.tstat.R`
- [x] `BSseq_class.R`
- [ ] `BSseq_utils.R`
- [ ] `BSseqStat_class.R`
- [ ] `BSseqTstat_class.R`
- [ ] `combine.R`
- [ ] `dmrFinder.R`
- [ ] `fisher.R`
- [ ] `getStats.R`
- [ ] `gof_stats.R`
- [ ] `hasGRanges.R`
- [ ] `permutations.R`
- [ ] `plotting.R`
- [x] `read.bismark.R`
- [ ] `read.bsmooth.R`
- [ ] `utils.R`
