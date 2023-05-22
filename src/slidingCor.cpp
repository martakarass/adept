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
// an NA value is returned for both the cor and sd. NA is returned for the sd
// instead of 0 so later values are also NA and not +-Inf.
//
// Also returns the longvec standard deviations for reuse in
// slidingCorCpp().
//
// [[Rcpp::export]]
List slidingCorStoreSdCpp(const NumericVector shortvec,
                          const NumericVector longvec,
                          double sd_shortvec) {
  // Get vector lengths and initialize output vector
  int length_longvec = longvec.size();
  int n = shortvec.size();
  int n_minus1 = n - 1;
  int out_length = length_longvec - n_minus1;
  NumericVector out(out_length);
  NumericVector sds(out_length);

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
    sds[0] = NA_REAL; // could also set to 0 to get Infs downstream
  } else {
    out[0] = (sum_products / n_minus1 - sum_longvec_current * term2) /
      sd_shortvec / sd_longvec_current;
    sds[0] = sd_longvec_current;
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
      sds[a] = NA_REAL;
    } else {
      out[a] = (sum_products / n_minus1 - sum_longvec_current * term2) /
        sd_shortvec / sd_longvec_current;
      sds[a] = sd_longvec_current;
    }
  }
  return List::create(_["cor"] = out, _["sds"] = sds);
}


// Same as slidingCorStoreSdCpp() but with precomputed longvec sds.
//
// Requires sliding standard deviations of longvec, for example those output
// by slidingCorStoreSdCpp, to be passed in.
//
// Does not return sd_longvec_current again.
//
// [[Rcpp::export]]
NumericVector slidingCorCpp(const NumericVector shortvec,
                            const NumericVector longvec,
                            double sd_shortvec,
                            const NumericVector sd_longvec_current) {
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
  double sum_products = 0;
  for (int b = 0; b < n; ++b) {
    double longvec_current_b = longvec[b];
    sum_products += longvec_current_b * shortvec[b];
  }
  out[0] = (sum_products / n_minus1 - sum_longvec_current * term2) /
    sd_shortvec / sd_longvec_current[0];

  for (int a = 1; a < out_length; ++a) {
    sum_longvec_current -= longvec[a - 1];
    sum_longvec_current += longvec[a + n_minus1];
    sum_products = 0;
    for (int b = 0; b < n; ++b) {
      double longvec_current_b = longvec[a + b];
      sum_products += longvec_current_b * shortvec[b];
    }
    out[a] = (sum_products / n_minus1 - sum_longvec_current * term2) /
      sd_shortvec / sd_longvec_current[a];
  }
  return out;
}
