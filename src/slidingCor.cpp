#include <Rcpp.h>
using namespace Rcpp;

// Optimized version of
// https://github.com/cran/dvmisc/blob/master/src/sliding_cor_c.cpp
//
// Removes redundant sum operation and keeps track of sum_longvec_current
// as a sliding sum, so instead of re-summing every window, just subtracts
// the first entry and adds the next entry.
//
// Also checks if the standard deviation is close to 0, in which case
// an NA value is returned. This fixes a bug in the original function where
// it would attempt to divide 0/0 and floating point error would lead to
// potentially spurious correlations.
//
// [[Rcpp::export]]
NumericVector slidingCorCpp(const NumericVector shortvec,
                              const NumericVector longvec, double sd_shortvec) {

  // Get vector lengths and initialize output vector
  int length_longvec = longvec.size();
  int n = shortvec.size();
  int n_minus1 = n - 1;
  int out_length = length_longvec - n_minus1;
  NumericVector out(out_length);

  // Calculate sum of short vector divided by n nsquared
  double term2 = sum(shortvec) / n / n_minus1;

  // Loop through longvec. For each segment calculate the sum of the products
  // with shortvec and record the covariance
  NumericVector longvec_current(n);
  longvec_current = longvec[Range(0, n_minus1)];
  double sum_longvec_current = sum(longvec_current);
  double mean_longvec_current = sum_longvec_current / n;
  double sum_products = 0;
  double ss_longvec_current = 0;
  for (int b = 0; b < n; ++b) {
    double longvec_current_b = longvec[b];
    sum_products += longvec_current_b * shortvec[b];
    ss_longvec_current += pow(longvec_current_b - mean_longvec_current, 2);
  }
  double sd_longvec_current = sqrt(ss_longvec_current / n_minus1);
  if (sd_longvec_current < 1e-10) {
    out[0] = NA_REAL;
  } else {
    out[0] = (sum_products / n_minus1 - sum_longvec_current * term2) /
      sd_shortvec / sd_longvec_current;
  }
  for (int a = 1; a < out_length; ++a) {
    sum_longvec_current -= longvec[a - 1];
    sum_longvec_current += longvec[a + n_minus1];
    mean_longvec_current = sum_longvec_current / n;
    sum_products = 0;
    ss_longvec_current = 0;
    for (int b = 0; b < n; ++b) {
      double longvec_current_b = longvec[a + b];
      sum_products += longvec_current_b * shortvec[b];
      ss_longvec_current += pow(longvec_current_b - mean_longvec_current, 2);
    }
    sd_longvec_current = sqrt(ss_longvec_current / n_minus1);
    if (sd_longvec_current < 1e-10) {
      out[a] = NA_REAL;
    } else {
      out[a] = (sum_products / n_minus1 - sum_longvec_current * term2) /
        sd_shortvec / sd_longvec_current;
    }
  }
  return out;
}
