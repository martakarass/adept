#include <Rcpp.h>
using namespace Rcpp;

// Get pmax and index of max from a list of vectors.
//
// Adapted from https://stackoverflow.com/a/66020829
//
// Vectors should be of the same length. Otherwise, will only go up to the
// length of the first vector.
//
// Returns a list where $pmax gives the same output as base R's pmax(),
// and $idx gives the index of the vector from which the max value was found.
//
// [[Rcpp::export]]
List pmaxIdxCpp(List args) {
  NumericVector tmp = args[0];
  NumericVector out = clone(tmp);
  NumericVector out_idx (out.length(), 1);
  int n_arg = args.length();
  int n_vec = out.length();
  for (int i = 1; i < n_arg; ++i) {
    NumericVector pa = args[i];
    for (int j = 0; j < n_vec; ++j) {
      if (pa[j] > out[j]) {
        out[j] = pa[j];
        out_idx[j] = i + 1;
      }
    }
  }
  return List::create(_["pmax"]=out,
                      _["idx"]=out_idx);
}
