#' SuSiE4I: Iterative SuSiE fitting for main and interaction effects (dense matrix only)
#'
#' This function performs iterative estimation of additive (Z) and main (X) effects, as well as
#' interaction effects between selected X components and Z, using SuSiE on sufficient statistics.
#'
#' For large sample size, it uses a blockwise strategy to compute `crossprod(X)` efficiently
#' without requiring the full X'X in memory at once.
#'
#' The algorithm alternates between updating the effects from Z (via linear projection),
#' from X (via susie_suff_stat), and interaction terms (via susie_suff_stat after constructing selected interactions).
#' After convergence, a post-hoc linear model of y ~ Z is fit with offset(etaX + etaW) to estimate
#' the significance of Z effects adjusted for X and interactions.
#'
#' @param X An n × p matrix of covariates.
#' @param Z An n × q matrix of covariates. If Z = NULL, only GxG interaction will be consider.
#' @param y An n-vector of responses.
#' @param n_threads Integer. Number of threads used in parallel computations, such as blockwise matrix cross-product or matrix–vector multiplication. Default is 1 (no parallelism). Increasing this can speed up computation on multicore machines when \code{is.large = TRUE}.
#' @param L_main The number of SuSiE effects for main effect X (default = 10).
#' @param L_int The number of SuSiE effects for interactions (default = 5).
#' @param select_env A indicator of whether selecting the environmental factors (default = F).
#' @param L_env If select_env = T, the number of SuSiE effects for interactions (default = 10).
#' @param noint_env The indices of environmental factors that would not be included in GEI search (default = NULL). When select_env = T, this parameter will not be considered.
#' @param max_iter Maximum number of outer iterations (default = 10).
#' @param max_eps Convergence threshold on η (default = 1e-5).
#' @param min_iter Minimum number of outer iterations (default = 3).
#' @param susie_iter Maximum number of iterations for each SuSiE fit (default = 300).
#' @param coverage_main The coverage to defining a credible set in fine-mapping main effects (default to 0.95).
#' @param coverage_env If select_env = T,the coverage to defining a credible set in fine-mapping environmental effects (default to 0.95).
#' @param coverage_int The coverage to defining a credible set in fine-mapping interaction effects (default to 0.85).
#' @param crossprodX The cross-product of X (default = F, and SuSiE4X will calculate it).
#' @param returnModel A logical indicator of whether to return the final model design matrix (default = False).
#' @param verbose A logical indicator of whether to display iteration diagnostics (default = TRUE).
#' @param ... Additional arguments passed to \code{susie_suff_stat}, such as \code{prior_variance}, \code{null_weight}, \code{standardize}, etc. Note that these tuning parameters will only be used in the fine-mapping of marginal effects.
#'
#' @return A list containing results from the iterative SuSiE fitting procedure:
#' \describe{
#'   \item{iter}{Total number of outer iterations completed before convergence or reaching maximum.}
#'   \item{error}{Numeric vector of max change in fitted values across iterations.}
#'   \item{fitX}{SuSiE result for main effects.}
#'   \item{fitW}{SuSiE result for interaction terms.}
#'   \item{fitZ}{Final linear model of y ~ Z with offset.}
#'   \item{fitJoint}{Final linear model of y ~ Z + credible sets of direct and interaction effects.}
#' }
#'
#' @importFrom Matrix colMeans crossprod
#' @importFrom stats var lm coef
#' @importFrom susieR susie_suff_stat coef.susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply
#' @importFrom graphics text

#'
#' @export
SuSiE4I <- function(X, Z=NULL, y, crossprodX=NULL,
                    n_threads = 4, L_main = 10, L_int = 5,
                    select_env=F, L_env=10, noint_env=NULL,
                    coverage_main = 0.95, coverage_int = 0.85, coverage_env=0.95,
                    max_iter = 15, max_eps = 1e-5, min_iter = 3,
                    susie_iter = 300, verbose = TRUE, returnModel=F,...) {


X=as.matrix(X)
X=large_scale(X)

if (is.null(colnames(X))){
colnames(X) <- paste0("X", seq_len(ncol(X)))
}

if(is.null(Z)==F){
Z=as.matrix(Z)
Z=large_scale(Z)
nameZ=colnames(Z)
if (is.null(colnames(Z))){
nameZ <- paste0("Z", seq_len(ncol(Z)))
colnames(Z)=nameZ
}

if(select_env==F){
Run_GGE(X = X, Z = Z, y = y, crossprodX=crossprodX,
        Lmain = L_main, Lint = L_int,noint_env=noint_env,
        coverage.main = coverage_main, coverage.int = coverage_int,
        max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
        susie.iter = susie_iter, verbose = verbose,n_threads=n_threads,
        returnModel=returnModel,...)
}else{
Run_GGE_Select(X = X, Z = Z, y = y, crossprodX=crossprodX,
        Lmain = L_main, Lint = L_int, Lenv=L_env,
        coverage.main = coverage_main, coverage.int = coverage_int,coverage.env=coverage_env,
        max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
        susie.iter = susie_iter, verbose = verbose,n_threads=n_threads,
        returnModel=returnModel,...)
}
}else{
Run_GG(X = X, y = y, crossprodX=crossprodX,
      Lmain = L_main, Lint = L_int,
      coverage.main = coverage_main, coverage.int = coverage_int,
      max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
      susie.iter = susie_iter, verbose = verbose,n_threads=n_threads,
      returnModel=returnModel,...)
}

}
