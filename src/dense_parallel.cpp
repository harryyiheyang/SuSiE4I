// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

//----------------------------------------
// 1. blockwise_crossprod: X^T X by row blocks
//----------------------------------------
// [[Rcpp::export]]
arma::mat blockwise_crossprod(const arma::mat& X, int n_threads = 1, int block_size = 10000) {
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat XtX(p, p, arma::fill::zeros);

  for (int i = 0; i < n; i += block_size) {
    int end = std::min(i + block_size, n);
    arma::mat Xi = X.rows(i, end - 1);
    XtX += Xi.t() * Xi;
  }

  return XtX;
}

//----------------------------------------
// 2. blockwise_crossprod2: X^T Z by row blocks
//----------------------------------------
// [[Rcpp::export]]
arma::mat blockwise_crossprod2(const arma::mat& X, const arma::mat& Z, int n_threads = 1, int block_size = 10000) {
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  if (X.n_rows != Z.n_rows) {
    Rcpp::stop("X and Z must have the same number of rows");
  }

  int n = X.n_rows;
  int p = X.n_cols;
  int q = Z.n_cols;
  arma::mat XtZ(p, q, arma::fill::zeros);

  for (int i = 0; i < n; i += block_size) {
    int end = std::min(i + block_size, n);
    arma::mat Xi = X.rows(i, end - 1);
    arma::mat Zi = Z.rows(i, end - 1);
    XtZ += Xi.t() * Zi;
  }

  return XtZ;
}

//----------------------------------------
// 3. large_scale: column-wise centering and scaling (in-place, but returned)
//----------------------------------------
// [[Rcpp::export]]
arma::mat large_scale(arma::mat X, bool center = true, bool scale = true, int n_threads = 1) {
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  int p = X.n_cols;

#pragma omp parallel for
  for (int j = 0; j < p; ++j) {
    arma::subview_col<double> col = X.col(j);  // modifies matrix X directly

    if (center) {
      double mu = arma::mean(col);
      col -= mu;
    }

    if (scale) {
      double sd = arma::stddev(col);
      if (sd > 0) {
        col /= sd;
      }
    }
  }

  return X;
}
