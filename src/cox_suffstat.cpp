// cox_suffstat.cpp
// Breslow risk-set quantities for Cox fine-mapping sufficient statistics.
// XtX = X' diag(a) X - B' diag(dev) B and Xty = X' M are assembled on the R
// side from these pieces, without forming an n x n matrix or copying X.
//
// Breslow ties: all individuals sharing the same time share the risk-set sums
// S0 / S1 evaluated at the end of their tied block. B has one row per event
// time (not per event); dev holds the number of events at that time.

#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Returns:
//   a    (n)       per-individual cumulative weight Lambda0(t_i) * exp(eta_i)
//   M    (n)       martingale residuals status_i - a_i
//   B    (nb x p)  risk-set weighted means of X at each event time
//   dev  (nb)      number of events at each event time
//   d    int       number of events
// The per-column risk-set scans are independent and run in parallel with
// OpenMP; they call no BLAS.
// [[Rcpp::export]]
Rcpp::List cox_riskset(const arma::mat& X,
                       const arma::vec& eta,
                       const arma::vec& time,
                       const arma::ivec& status,
                       int n_threads = 1) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (eta.n_elem != n || time.n_elem != n || status.n_elem != n) {
    Rcpp::stop("eta, time, and status must have one value per row of X.");
  }
  const int nt = std::max(1, n_threads);

  // descending sort by time; risk set then grows monotonically as we advance
  arma::uvec ord = arma::sort_index(time, "descend");

  // eta shift cancels in every risk-set ratio and in a_i.
  const double eta_max = eta.max();
  arma::vec w = arma::exp(eta - eta_max);

  // ----- pass 1: tied blocks, S0 at each block end, events per block -----
  std::vector<arma::uword> blk_end;   // exclusive end (in sorted order) of event blocks
  std::vector<double> blk_S0, blk_time;
  std::vector<double> blk_dev;
  double S0 = 0.0;
  arma::uword i = 0;
  while (i < n) {
    arma::uword j = i;
    const double t_block = time(ord(i));
    arma::uword dev = 0;
    while (j < n && time(ord(j)) == t_block) {
      S0 += w(ord(j));
      if (status(ord(j)) == 1) ++dev;
      ++j;
    }
    if (dev > 0) {
      blk_end.push_back(j);
      blk_S0.push_back(S0);
      blk_time.push_back(t_block);
      blk_dev.push_back(static_cast<double>(dev));
    }
    i = j;
  }
  const arma::uword nb = blk_end.size();
  double d = 0.0;
  for (arma::uword b = 0; b < nb; ++b) d += blk_dev[b];

  // ----- pass 2: per-column risk-set means at event blocks -----
  arma::mat B(nb, p, arma::fill::none);
#pragma omp parallel for schedule(static) num_threads(nt)
  for (arma::uword c = 0; c < p; ++c) {
    const double* xc = X.colptr(c);
    double S1 = 0.0;
    arma::uword r = 0;
    for (arma::uword b = 0; b < nb; ++b) {
      for (; r < blk_end[b]; ++r) {
        const arma::uword idx = ord(r);
        S1 += w(idx) * xc[idx];
      }
      B(b, c) = S1 / blk_S0[b];
    }
  }

  // ----- a_i = exp(eta_i) * Lambda0(t_i), Breslow increment per event time -----
  arma::uvec ord_asc = arma::sort_index(time, "ascend");
  arma::vec a(n, arma::fill::zeros);
  double Lam = 0.0;
  long bptr = static_cast<long>(nb) - 1;
  for (arma::uword r = 0; r < n; ++r) {
    const arma::uword idx = ord_asc(r);
    const double ti = time(idx);
    while (bptr >= 0 && blk_time[bptr] <= ti) {
      Lam += blk_dev[bptr] / blk_S0[bptr];
      --bptr;
    }
    a(idx) = w(idx) * Lam;
  }
  arma::vec M = arma::conv_to<arma::vec>::from(status) - a;

  return Rcpp::List::create(
    Rcpp::Named("a")   = a,
    Rcpp::Named("M")   = M,
    Rcpp::Named("B")   = B,
    Rcpp::Named("dev") = arma::vec(blk_dev),
    Rcpp::Named("d")   = static_cast<int>(d)
  );
}
