#include <Rcpp.h>
using namespace Rcpp;

// TODO: It should be possible to do some template magic to simplify these 4
//       definitions into 1 (or 2).
// TODO: Move to R's C API to avoid dependency on Rcpp?

//' Internal helper functions to check that the Cov and M matrix objects in a
//' BSseq object are valid.
//'
//' Equivalent to
//' any(M < 0) || any(M > Cov) || any(is.na(M)) || any(is.na(Cov)) ||
//' any(is.infinite(Cov)) but
//' - **Much** more memory efficient because there are no intermediate vectors
//' - Many times faster when M and Cov are valid (requiring 5 less iterations
//'   over M or Cov and no memory allocations) and many, many times faster if
//'   M or Cov are invalid for any reason (especially so if the invalid element
//'   occurs towards the beginning of M and Cov or is triggered by one of the
//'   later checks, e.g, Inf value in Cov)
//'
//' IMPORTANT: Do not pass DelayedMatrix objects to these functions, they are
//' only valid for matrix objects
//'
//' @param Cov The Cov matrix
//' @param M The M matrix
//'
//' @keywords internal
//'
//' @return TRUE if M and Cov are valid, FALSE otherwise
//'
//' @keywords internal
//'
// [[Rcpp::export(".validIntegerCovAndIntegerM")]]
bool validIntegerCovAndIntegerM(IntegerMatrix Cov, IntegerMatrix M) {
    // TODO: Check dimensions are compatible
    int length = Cov.length();
    for (int i = 0; i < length; i++) {
        // NOTE: No need to check is.infinite(Cov[i]) because by definition an
        //       IntegerMatrix cannot contain Inf
        if (M[i] < 0 || M[i] > Cov[i] || IntegerMatrix::is_na(M[i]) ||
            IntegerMatrix::is_na(Cov[i])) {
            return false;
        }
    }
    return true;
}

// [[Rcpp::export(".validIntegerCovAndNumericM")]]
bool validIntegerCovAndNumericM(IntegerMatrix Cov, NumericMatrix M) {
    // TODO: Check dimensions are compatible
    int length = Cov.length();
    for (int i = 0; i < length; i++) {
        // NOTE: No need to check is.infinite(Cov[i]) because by definition an
        //       IntegerMatrix cannot contain Inf
        if (M[i] < 0 || M[i] > Cov[i] || NumericMatrix::is_na(M[i]) ||
            IntegerMatrix::is_na(Cov[i])) {
            return false;
        }
    }
    return true;
}

// [[Rcpp::export(".validNumericCovAndIntegerM")]]
bool validNumericCovAndIntegerM(NumericMatrix Cov, IntegerMatrix M) {
    // TODO: Check dimensions are compatible
    int length = Cov.length();
    for (int i = 0; i < length; i++) {
        // TODO: See http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2014-January/007033.html
        // for why is_infinite() is kinda weird; see if can be fixed
        if (M[i] < 0 || M[i] > Cov[i] || IntegerMatrix::is_na(M[i]) ||
            NumericMatrix::is_na(Cov[i]) ||
            traits::is_infinite<REALSXP>(Cov[i])) {
            return false;
        }
    }
    return true;
}

// [[Rcpp::export(".validNumericCovAndNumericM")]]
bool validNumericCovAndNumericM(NumericMatrix Cov, NumericMatrix M) {
    // TODO: Check dimensions are compatible
    int length = Cov.length();
    for (int i = 0; i < length; i++) {
        // TODO: See http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2014-January/007033.html
        // for why is_infinite() is kinda weird; see if can be fixed
        if (M[i] < 0 || M[i] > Cov[i] || NumericMatrix::is_na(M[i]) ||
            NumericMatrix::is_na(Cov[i]) ||
            traits::is_infinite<REALSXP>(Cov[i])) {
            return false;
        }
    }
    return true;
}
