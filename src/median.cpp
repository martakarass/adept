#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double medianCpp(NumericVector x){
  return(median(x));
}
