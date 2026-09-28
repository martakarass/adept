#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <float.h>
#include <math.h>

static SEXP adept_sliding_cov(SEXP short_arg, SEXP long_arg) {
  const int n = LENGTH(short_arg);
  const int n_long = LENGTH(long_arg);
  const int n_out = n_long - n + 1;
  const double *short_vec = REAL(short_arg);
  const double *long_vec = REAL(long_arg);
  SEXP out = PROTECT(allocVector(REALSXP, n_out));
  double sum_short = 0.0;
  for (int j = 0; j < n; ++j) sum_short += short_vec[j];
  const double term2 = sum_short / n / (n - 1);

  for (int i = 0; i < n_out; ++i) {
    double sum_long = 0.0;
    double sum_products = 0.0;
    for (int j = 0; j < n; ++j) {
      const double value = long_vec[i + j];
      sum_long += value;
      sum_products += value * short_vec[j];
    }
    REAL(out)[i] = sum_products / (n - 1) - sum_long * term2;
    if ((i & 4095) == 4095) R_CheckUserInterrupt();
  }
  UNPROTECT(1);
  return out;
}

static SEXP adept_sliding_cor(SEXP short_arg, SEXP long_arg) {
  const int n = LENGTH(short_arg);
  const int n_long = LENGTH(long_arg);
  const int n_out = n_long - n + 1;
  const double *short_vec = REAL(short_arg);
  const double *long_vec = REAL(long_arg);
  SEXP out = PROTECT(allocVector(REALSXP, n_out));

  double sum_short = 0.0;
  for (int j = 0; j < n; ++j) sum_short += short_vec[j];
  const double mean_short = sum_short / n;
  double ss_short = 0.0;
  for (int j = 0; j < n; ++j) {
    const double d = short_vec[j] - mean_short;
    ss_short += d * d;
  }
  const double sd_short = sqrt(ss_short / (n - 1));
  const double term2 = sum_short / n / (n - 1);
  const double min_sd = sqrt(DBL_EPSILON);
  if (sd_short < min_sd) {
    for (int i = 0; i < n_out; ++i) REAL(out)[i] = NA_REAL;
    UNPROTECT(1);
    return out;
  }

  for (int i = 0; i < n_out; ++i) {
    double sum_long = 0.0;
    double sum_products = 0.0;
    for (int j = 0; j < n; ++j) {
      const double value = long_vec[i + j];
      sum_long += value;
      sum_products += value * short_vec[j];
    }
    const double mean_long = sum_long / n;
    double ss_long = 0.0;
    for (int j = 0; j < n; ++j) {
      const double d = long_vec[i + j] - mean_long;
      ss_long += d * d;
    }
    const double sd_long = sqrt(ss_long / (n - 1));
    if (sd_long < min_sd) {
      REAL(out)[i] = NA_REAL;
    } else {
      REAL(out)[i] = (sum_products / (n - 1) - sum_long * term2) /
        sd_short / sd_long;
    }
    if ((i & 4095) == 4095) R_CheckUserInterrupt();
  }
  UNPROTECT(1);
  return out;
}

static SEXP adept_moving_mean(SEXP x_arg, SEXP window_arg) {
  const int n = LENGTH(x_arg);
  const int window = asInteger(window_arg);
  const int n_out = n - window + 1;
  const double *x = REAL(x_arg);
  const double inverse = 1.0 / window;
  SEXP out = PROTECT(allocVector(REALSXP, n_out));
  double sum = 0.0;
  for (int i = 0; i < window; ++i) sum += x[i];
  REAL(out)[0] = sum * inverse;
  for (int i = window; i < n; ++i) {
    sum += x[i] - x[i - window];
    REAL(out)[i - window + 1] = sum * inverse;
    if ((i & 16383) == 16383) R_CheckUserInterrupt();
  }
  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef CallEntries[] = {
  {"adept_sliding_cov", (DL_FUNC) &adept_sliding_cov, 2},
  {"adept_sliding_cor", (DL_FUNC) &adept_sliding_cor, 2},
  {"adept_moving_mean", (DL_FUNC) &adept_moving_mean, 2},
  {NULL, NULL, 0}
};

void R_init_adept(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
