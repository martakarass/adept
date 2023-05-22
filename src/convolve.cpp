#include <Rcpp.h>
using namespace Rcpp;

// Taken from https://dirk.eddelbuettel.com/code/rcpp/Rcpp-attributes.pdf
// [[Rcpp::export]]
NumericVector convolveCpp(const NumericVector a, const NumericVector b) {
  int na = a.size(), nb = b.size();
  int nab = na + nb - 1;
  NumericVector xab(nab);
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++)
      xab[i + j] += a[i] * b[j];
  return xab;
}
