#include <R.h>
#include <Rinternals.h>

SEXP makeBinEndpoints(SEXP ss, SEXP max_bins, SEXP min_loci, SEXP min_width) {

    // Allocate vectors
    // TODO: This overallocates s and e; alternatively, re-write using C++ with
    //       std::vector and its push_back() method
    int cmax_bins = asInteger(max_bins);
    SEXP s = PROTECT(allocVector(INTSXP, cmax_bins));
    SEXP e = PROTECT(allocVector(INTSXP, cmax_bins));
    SEXP val = PROTECT(allocVector(VECSXP, 2));

    // Clear contents of s and e
    memset(INTEGER(s), 0, cmax_bins * sizeof(int));
    memset(INTEGER(e), 0, cmax_bins * sizeof(int));

    // Get pointers
    int *pss, *ps, *pe;
    pss = INTEGER(ss);
    ps = INTEGER(s);
    pe = INTEGER(e);

    // Coerce from length one R vectors into C scalars
    int cmin_loci, cmin_width;
    cmin_loci = asInteger(min_loci);
    cmin_width = asInteger(min_width);

    // Loop over ss and construct appropriate s and e
    int i = 0;
    int j = 0;
    int k = 0;
    R_xlen_t n = xlength(ss);
    while (i < n) {
        k = cmin_loci;
        while ((i + k) < n) {
            if ((pss[i + k - 1] - pss[i]) >= cmin_width) {
                ps[j] = pss[i];
                pe[j] = pss[i + k - 1L];
                j = j + 1;
                i = i + k;
                break;
            } else {
                k = k + 1;
            }
        }
        i = i + 1;
    }

    // Construct list, val
    // TODO: Figure out how to return s[0:j] and e[0:j] (or whatever the last
    //       valid element is)
    // TODO: Figure out how to include names on val
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, e);

    // Don't forget to UNPROTECT() every PROTECT()
    UNPROTECT(3);

    return val;
}
