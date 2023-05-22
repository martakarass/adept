#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double median_cpp(NumericVector x){
  return(median(x));
}
