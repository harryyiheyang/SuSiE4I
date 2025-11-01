#' Blockwise crossproduct (X'X or X'Z)
#'
#' Computes X'X or X'Z in memory-efficient blocks.
#'
#' @param X Numeric matrix.
#' @param Z Optional second matrix. If provided, computes X'Z.
#' @param n_threads Number of threads.
#' @param block_size Block size for computation.
#' @export
blockwise_crossprod <- function(X, Z = NULL, n_threads = 1L, block_size = 10000L) {
if (is.null(Z)) {
return(.Call(`_SuSiE4I_blockwise_crossprod`, X, n_threads, block_size))
} else {
return(.Call(`_SuSiE4I_blockwise_crossprod2`, X, Z, n_threads, block_size))
}
}
