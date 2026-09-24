#' Crossprod using the system BLAS or an OpenMP fallback
#'
#' Computes `crossprod(X)` or `crossprod(X, Z)`. With an optimized BLAS
#' (OpenBLAS, MKL, FlexiBLAS, ...), it calls `base::crossprod()`, whose
#' threading is controlled by the BLAS (e.g. `OPENBLAS_NUM_THREADS`). With R's
#' reference BLAS and `n_threads > 1`, the output is split into column blocks
#' computed in parallel with OpenMP.
#'
#' @param X Numeric matrix.
#' @param Z Optional second matrix. If provided, computes `crossprod(X, Z)`.
#' @param n_threads Number of OpenMP threads for the reference-BLAS path.
#'   Defaults to 4.
#' @param block_size Unused; retained for backward compatibility.
#' @export
blockwise_crossprod <- function(X, Z = NULL, n_threads = 4L, block_size = 10000L) {
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (!is.null(Z) && (!is.matrix(Z) || !is.numeric(Z))) {
    stop("Z must be a numeric matrix.")
  }
  if (!is.numeric(n_threads) || length(n_threads) != 1L ||
      !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a positive numeric scalar.")
  }
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size < 1) {
    stop("block_size must be a positive numeric scalar.")
  }
  n_threads <- as.integer(n_threads)
  block_size <- as.integer(block_size)

  if (n_threads <= 1L || !reference_blas()) {
    if (is.null(Z)) return(base::crossprod(X))
    return(base::crossprod(X, Z))
  }

  if (is.null(Z)) {
    return(blockwise_crossprod_cpp(X, n_threads, block_size))
  }
  blockwise_crossprod2_cpp(X, Z, n_threads, block_size)
}

reference_blas <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      blas <- tryCatch(extSoftVersion()[["BLAS"]], error = function(e) "")
      cached <<- !nzchar(blas) || grepl("Rblas", blas, fixed = TRUE)
    }
    cached
  }
})
