#' large_scale: Center and scale columns of a matrix (fast version)
#'
#' This function wraps the C++ version of scale(), implemented with RcppArmadillo.
#' It centers and/or scales each column of a numeric matrix.
#'
#' @param X A numeric matrix.
#' @param center Logical. If TRUE, center each column by its mean.
#' @param scale Logical. If TRUE, scale each column by its standard deviation.
#' @param n_threads The number of cores. Defaults to 1.
#'
#' @return A matrix with the same dimensions as X, standardized accordingly.
#' @export
large_scale <- function(X,  n_threads=1, center = TRUE, scale = TRUE) {
  if (!is.matrix(X)) stop("X must be a numeric matrix")
  NAM=colnames(X)
  X=.Call('_SuSiE4I_large_scale', X, center, scale,n_threads)
  colnames(X)=NAM
  return(X)
}
