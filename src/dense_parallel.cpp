// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

static int column_block_width(int p, int n_threads) {
  int nb = std::max(1, 4 * n_threads);
  return std::max(32, (p + nb - 1) / nb);
}

//----------------------------------------
// 1. X^T X, parallel over column blocks of the output (upper triangle)
//----------------------------------------
// [[Rcpp::export]]
arma::mat blockwise_crossprod_cpp(const arma::mat& X, int n_threads = 4, int block_size = 10000) {
  const arma::uword n = X.n_rows;
  const int p = X.n_cols;
  arma::mat XtX(p, p, arma::fill::zeros);
  const int cb = column_block_width(p, n_threads);
  const int nb = (p + cb - 1) / cb;
  double* base = const_cast<double*>(X.memptr());

#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
  for (int b = 0; b < nb; ++b) {
    const int c0 = b * cb;
    const int c1 = std::min(c0 + cb, p);
    const arma::mat A(base, n, c1, false, true);
    const arma::mat B(base + static_cast<arma::uword>(c0) * n, n, c1 - c0, false, true);
    XtX.submat(0, c0, c1 - 1, c1 - 1) = A.t() * B;
  }

  return arma::symmatu(XtX);
}

//----------------------------------------
// 2. X^T Z, parallel over row blocks of the output (columns of X)
//----------------------------------------
// [[Rcpp::export]]
arma::mat blockwise_crossprod2_cpp(const arma::mat& X, const arma::mat& Z, int n_threads = 4, int block_size = 10000) {
  if (X.n_rows != Z.n_rows) {
    Rcpp::stop("X and Z must have the same number of rows");
  }

  const arma::uword n = X.n_rows;
  const int p = X.n_cols;
  arma::mat XtZ(p, Z.n_cols, arma::fill::zeros);
  const int cb = column_block_width(p, n_threads);
  const int nb = (p + cb - 1) / cb;
  double* base = const_cast<double*>(X.memptr());

#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
  for (int b = 0; b < nb; ++b) {
    const int c0 = b * cb;
    const int c1 = std::min(c0 + cb, p);
    const arma::mat A(base + static_cast<arma::uword>(c0) * n, n, c1 - c0, false, true);
    XtZ.rows(c0, c1 - 1) = A.t() * Z;
  }

  return XtZ;
}

//----------------------------------------
// 3. large_scale: column-wise centering and scaling (in-place, but returned)
//----------------------------------------
// [[Rcpp::export]]
arma::mat large_scale_cpp(arma::mat X, bool center = true, bool scale = true, int n_threads = 4) {
  int p = X.n_cols;

#pragma omp parallel for num_threads(n_threads)
  for (int j = 0; j < p; ++j) {
    arma::subview_col<double> col = X.col(j);

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
