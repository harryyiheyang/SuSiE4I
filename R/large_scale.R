#' large_scale: Center and scale columns of a matrix (fast version)
#'
#' This function wraps the C++ version of scale(), implemented with RcppArmadillo.
#' It centers and/or scales each column of a numeric matrix.
#'
#' @param X A numeric matrix.
#' @param center Logical. If TRUE, center each column by its mean.
#' @param scale Logical. If TRUE, scale each column by its standard deviation.
#' @param n_threads The number of threads. Defaults to 4.
#'
#' @return A matrix with the same dimensions as X, standardized accordingly.
#' @export
large_scale <- function(X, n_threads = 4L, center = TRUE, scale = TRUE) {
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (!is.numeric(n_threads) || length(n_threads) != 1L ||
      !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a positive numeric scalar.")
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("scale must be TRUE or FALSE.")
  }

  DN <- dimnames(X)
  X <- large_scale_cpp(X, center, scale, as.integer(n_threads))
  dimnames(X) <- DN
  X
}
