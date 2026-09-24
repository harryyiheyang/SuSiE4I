// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

// Row-chunked X' diag(w) X and X' M. Only one block_size x p chunk is copied
// at a time; the weighted n x p matrix is never formed. X'WX is accumulated in
// place by dsyrk (threaded by an optimized BLAS). With R's reference BLAS,
// use_omp splits each chunk's update into OpenMP column blocks instead.
// [[Rcpp::export]]
Rcpp::List weighted_crossprod_cpp(const arma::mat& X, const arma::vec& w,
                                  const arma::mat& M, int block_size = 10000,
                                  int n_threads = 1, bool use_omp = false) {
  const arma::uword n = X.n_rows;
  const int p = X.n_cols;
  const arma::uword m = M.n_cols;
  if (w.n_elem != n) Rcpp::stop("length(w) must equal nrow(X).");
  if (m > 0 && M.n_rows != n) Rcpp::stop("nrow(M) must equal nrow(X).");
  if (arma::any(w < 0)) Rcpp::stop("w must be non-negative.");

  arma::mat XtWX(p, p, arma::fill::zeros);
  arma::mat XtM(p, m, arma::fill::zeros);
  if (n == 0 || p == 0) return Rcpp::List::create(Rcpp::Named("XtWX") = XtWX, Rcpp::Named("XtM") = XtM);

  const arma::uword bs = static_cast<arma::uword>(std::max(1, block_size));
  const int nt = std::max(1, n_threads);
  const int cb = std::max(32, (p + 4 * nt - 1) / (4 * nt));
  const int nb = (p + cb - 1) / cb;
  const double one = 1.0;
  arma::blas_int pp = p;

  for (arma::uword r0 = 0; r0 < n; r0 += bs) {
    const arma::uword r1 = std::min(r0 + bs, n);
    arma::mat Xi = X.rows(r0, r1 - 1);
    if (m > 0) XtM += Xi.t() * M.rows(r0, r1 - 1);
    Xi.each_col() %= arma::sqrt(w.subvec(r0, r1 - 1));
    arma::blas_int k = r1 - r0;

    if (use_omp && nt > 1) {
#pragma omp parallel for schedule(dynamic) num_threads(nt)
      for (int b = 0; b < nb; ++b) {
        const int c0 = b * cb;
        const int c1 = std::min(c0 + cb, p);
        arma::blas_int mm = c1, nn = c1 - c0;
        arma::blas::gemm<double>("T", "N", &mm, &nn, &k, &one, Xi.memptr(), &k,
                                 Xi.colptr(c0), &k, &one, XtWX.colptr(c0), &pp);
      }
    } else {
      arma::blas::syrk<double>("U", "T", &pp, &k, &one, Xi.memptr(), &k,
                               &one, XtWX.memptr(), &pp);
    }
  }

  return Rcpp::List::create(Rcpp::Named("XtWX") = arma::symmatu(XtWX),
                            Rcpp::Named("XtM") = XtM);
}
