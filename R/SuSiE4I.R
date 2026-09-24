#' SuSiE4I: Iterative SuSiE fitting for main and interaction effects
#'
#' Fits main and interaction effects for Gaussian, GLM/mgcv,
#' ordered-categorical, and Cox outcomes.
#' Gaussian outer-loop joint refits use an exact sufficient-statistic ridge
#' solve with the current penalty `phi / V`. After the final design is fixed,
#' the dispersion is iterated to its ridge fixed point and one final mgcv fit
#' is produced for the returned formal model and downstream LBF calculations.
#'
#' @param X An n by p numeric predictor matrix.
#' @param Z Optional n by q environmental-covariate matrix.
#' @param y Response vector, ordered factor, or `survival::Surv` object.
#' @param status Optional event indicator for Cox models when `y` is a
#'   survival-time vector.
#' @param family Outcome family. May be a supported family object or dispatch
#'   string. Cumulative-link models use `"clm_logit"` or `"clm_probit"`.
#' @param mgcv_model `NULL`, `"gam"`, or `"bam"` for mgcv refits. `NULL`
#'   selects `"gam"`. For Gaussian outcomes this controls only the final
#'   formal fit; intermediate joint refits use the exact ridge solver.
#' @param crossprodX Optional precomputed cross-product of `X` for Gaussian paths.
#' @param scale_data Whether to standardize `X` and `Z`.
#' @param n_threads Number of threads used for cross-products.
#' @param L_main Number of main-effect SuSiE components.
#' @param L_int Number of interaction SuSiE components.
#' @param select_env Whether to fine-map columns of `Z`. Supported for all
#'   outcome paths.
#' @param L_env Number of environmental SuSiE components.
#' @param noint_env Indices of `Z` columns excluded from interaction construction.
#' @param include_x_squared Whether to include squared main-effect summaries in
#'   the interaction design.
#' @param susie_para_main Named `susieR::susie_ss()` options for main effects.
#'   `prior_variance` is an absolute coefficient prior variance; the legacy
#'   `scaled_prior_variance` name is accepted with a warning and interpreted
#'   identically. Set `estimate_prior_variance = FALSE` to fix V or `TRUE` to
#'   initialize its estimation. Through `min_iter`, the prior-variance value and
#'   estimation switch are ignored in favor of fixed `scaled_prior_variance = 2`
#'   and `estimate_prior_variance = FALSE`; other `susie_para_main` settings
#'   remain active. User prior-variance settings take effect after `min_iter`.
#'   This setting applies only to the main stage.
#' @param susie_para_int SuSiE options for interaction effects, with the same
#'   prior-variance semantics as `susie_para_main`; it is independent of the
#'   main-stage setting.
#' @param susie_para_env SuSiE options for environmental effects, with the same
#'   prior-variance semantics as `susie_para_main`; it is independent of the
#'   main-stage setting.
#' @param max_iter Maximum number of outer iterations.
#' @param max_eps Outer-loop convergence threshold.
#' @param min_iter End of the fixed-V warm-up and minimum outer iterations
#'   before convergence can be declared.
#' @param L.init Number of `X` variables used in the low-dimensional fitting step.
#' @param x_noncs_var Minimum projected residual-variance fraction required to
#'   include `Main_noncs_res`.
#' @param w_noncs_var Minimum residual-variance fraction required to include
#'   `Int_noncs_res`.
#' @param noncs_max_abs_cor Maximum allowed absolute correlation between a
#'   non-CS term and applicable refit covariates.
#' @param suff_block_size Row-block size for sufficient-statistic cross-products.
#' @param verbose Whether to print iteration diagnostics.
#' @param returnModel Whether to return the final design matrix.
#'
#' @return A list containing the fitted SuSiE objects, final joint model,
#'   `main_discoveries`, `interaction_discoveries`, optional
#'   `environment_discoveries`, component summaries, diagnostics, and optional
#'   model data. Non-CS refit terms are nuisance variables rather than
#'   discoveries; refit p-values are descriptive rather than post-selection
#'   inference.
#'
#' @importFrom Matrix crossprod
#' @importFrom stats var lm glm coef binomial gaussian cor cov2cor reformulate sd
#' @importFrom mgcv gam bam nb tw betar scat
#' @importFrom susieR susie_ss coef.susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixEigen
#' @importFrom graphics plot text
#' @useDynLib SuSiE4I, .registration = TRUE
#'
#' @export
SuSiE4I <- function(X, Z = NULL, y, status = NULL, family = NULL,
                    mgcv_model = NULL,
                    crossprodX = NULL, scale_data = TRUE,
                    n_threads = 4,
                    L_main = 10, L_int = 5,
                    select_env = FALSE, L_env = 10, noint_env = NULL,
                    include_x_squared = FALSE,
                    susie_para_main = NULL,
                    susie_para_int = NULL,
                    susie_para_env = NULL,
                    max_iter = 10, max_eps = 1e-5, min_iter = 2,
                    L.init = 1,
                    x_noncs_var = 0.1,
                    w_noncs_var = 0.1,
                    noncs_max_abs_cor = 0.9,
                    suff_block_size = 10000L,
                    verbose = TRUE, returnModel = FALSE) {

if (missing(y) || is.null(y)) stop("y must be provided.")
susie_para_main <- .resolve_susie_para(susie_para_main, "susie_para_main")
susie_para_int <- .resolve_susie_para(susie_para_int, "susie_para_int")
susie_para_env <- .resolve_susie_para(susie_para_env, "susie_para_env")
suff_block_size <- validate_suff_block_size(suff_block_size)

is_surv_y <- inherits(y, "Surv")
if (is_surv_y) {
status <- as.integer(y[, 2])
y <- as.numeric(y[, 1])
}

X <- as.matrix(X)
if (!is.numeric(X)) stop("X must be numeric.")
if (ncol(X) == 0) stop("X has zero columns.")
n <- nrow(X)
if (length(y) != n) stop("Length(y) must equal nrow(X).")
if (scale_data) X <- large_scale(X)
if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

if (!is.null(Z)) {
Z <- as.matrix(Z)
if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
if (scale_data) Z <- large_scale(Z)
if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
bad_z <- grepl("^Main_", colnames(Z))
if (any(bad_z)) {
warning("Renaming Z column(s) starting with 'Main_' to 'MaIn_' to avoid collision with main-effect labels.")
colnames(Z)[bad_z] <- sub("^Main_", "MaIn_", colnames(Z)[bad_z])
}
}

is_binary_response <- function(v) {
vv <- unique(stats::na.omit(v))
length(vv) <= 2L && all(vv %in% c(0, 1))
}

family_string <- NULL
if (is.character(family)) {
if (length(family) != 1L) stop("family must be a scalar string or a GLM family object.")
family_string <- tolower(family)
}

if (is.null(family)) {
if (is_surv_y || !is.null(status)) {
family_string <- "cox"
} else if (is.ordered(y)) {
stop("For an ordered response, specify family as clm_logit, clm_probit, clm_cloglog, clm_loglog, clm_cauchit, or ocat.")
} else if (is_binary_response(y)) {
family <- binomial(link = "logit")
} else {
family <- gaussian()
}
}

if (!is.null(family_string) && family_string %in% c("cox", "coxph", "survival")) {
if (is.null(status)) stop("status is required for Cox models unless y is a Surv object.")
status <- as.integer(status)
if (length(status) != n) stop("Length(status) must equal nrow(X).")
if (is.null(Z)) {
return(Run_GG_Cox(
X = X, y = y, status = status,
Lmain = L_main, Lint = L_int,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
L.init = L.init,
include_x_squared = include_x_squared,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
if (select_env) {
return(Run_GGE_Select_Cox(
X = X, Z = Z, y = y, status = status,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int, Lenv = L_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
susie_para_env = susie_para_env,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
return(Run_GGE_Cox(
X = X, Z = Z, y = y, status = status,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int, noint_env = noint_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
L.init = L.init,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}

clm_links <- c("logit", "probit", "cloglog", "loglog", "cauchit")
clm_link <- NULL
if (!is.null(family_string) && startsWith(family_string, "clm_")) {
clm_link <- sub("^clm_", "", family_string)
if (!clm_link %in% clm_links) {
stop("Unsupported CLM family. Use clm_logit, clm_probit, clm_cloglog, clm_loglog, or clm_cauchit.")
}
}
if (identical(family_string, "clm")) {
stop("family = 'clm' is incomplete. Specify clm_logit, clm_probit, clm_cloglog, clm_loglog, or clm_cauchit.")
}
is_clm_string <- !is.null(clm_link)
is_ocat_string <- !is.null(family_string) &&
family_string %in% c("ocat", "ordinal", "ordered", "ordered.categorical",
                     "ordered_categorical")
if (is_clm_string || is_ocat_string || ocat_is_family(family)) {
ocat_family <- if (ocat_is_family(family)) family else NULL
ordinal_link <- if (is_clm_string) clm_link else "logit"
if (identical(ordinal_link, "logit")) {
if (is.null(ocat_family)) {
y_info <- ocat_prepare_response(y)
ocat_family <- mgcv::ocat(R = y_info$ncat)
}
if (is.null(Z)) {
return(Run_GG_OCAT(
X = X, y = y, family = ocat_family,
include_x_squared = include_x_squared,
mgcv_model = mgcv_model,
Lmain = L_main, Lint = L_int,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
L.init = L.init,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
if (select_env) {
return(Run_GGE_Select_OCAT(
X = X, Z = Z, y = y, family = ocat_family,
include_x_squared = include_x_squared,
mgcv_model = mgcv_model,
Lmain = L_main, Lint = L_int, Lenv = L_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
susie_para_env = susie_para_env,
L.init = L.init,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
return(Run_GGE_OCAT(
X = X, Z = Z, y = y, family = ocat_family,
include_x_squared = include_x_squared,
mgcv_model = mgcv_model,
Lmain = L_main, Lint = L_int, noint_env = noint_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
L.init = L.init,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}

if (is.null(Z)) {
return(Run_GG_CLM(
X = X, y = y, clm_link = ordinal_link,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
if (select_env) {
return(Run_GGE_Select_CLM(
X = X, Z = Z, y = y, clm_link = ordinal_link,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int, Lenv = L_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
susie_para_env = susie_para_env,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
return(Run_GGE_CLM(
X = X, Z = Z, y = y, clm_link = ordinal_link,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int, noint_env = noint_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}

if (!is.null(family_string)) {
if (family_string %in% c("gaussian", "normal", "linear")) {
family <- gaussian()
} else if (family_string %in% c("binomial", "logit", "logistic")) {
family <- binomial(link = "logit")
} else if (family_string %in% c("negbin", "nb", "negative.binomial", "negative_binomial")) {
family <- mgcv::nb(theta = NULL)
} else if (family_string %in% c("lognormal", "lnormal", "log-normal")) {
stop("family = 'lognormal' is not supported. Use log(y) with family = 'gaussian' instead.")
} else {
stop("Unsupported family string. Use gaussian, binomial, logit, negbin, ocat, cox, or a clm_<link> family.")
}
}

if (!inherits(family, "family")) {
stop("family must be NULL, a supported string, or a GLM family object.")
}

if (identical(family$family, "gaussian") && identical(family$link, "identity")) {
if (is.null(Z)) {
return(Run_GG(
X = X, y = y, mgcv_model = mgcv_model, crossprodX = crossprodX,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
L.init = L.init,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
if (select_env) {
return(Run_GGE_Select(
X = X, Z = Z, y = y, mgcv_model = mgcv_model, crossprodX = crossprodX,
L.init = L.init,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int, Lenv = L_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
susie_para_env = susie_para_env,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
return(Run_GGE(
X = X, Z = Z, y = y, mgcv_model = mgcv_model, crossprodX = crossprodX,
include_x_squared = include_x_squared,
Lmain = L_main, Lint = L_int, noint_env = noint_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
L.init = L.init,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}

if (is.null(Z)) {
return(Run_GG_GLM(
X = X, y = y, family = family,
include_x_squared = include_x_squared,
mgcv_model = mgcv_model,
Lmain = L_main, Lint = L_int,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
L.init = L.init,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
if (select_env) {
return(Run_GGE_Select_GLM(
X = X, Z = Z, y = y, family = family,
include_x_squared = include_x_squared,
mgcv_model = mgcv_model,
Lmain = L_main, Lint = L_int, Lenv = L_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
susie_para_env = susie_para_env,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
))
}
Run_GGE_GLM(
X = X, Z = Z, y = y, family = family,
include_x_squared = include_x_squared,
mgcv_model = mgcv_model,
Lmain = L_main, Lint = L_int, noint_env = noint_env,
max.iter = max_iter, max.eps = max_eps, min.iter = min_iter,
susie_para_main = susie_para_main,
susie_para_int = susie_para_int,
L.init = L.init,
x_noncs_var = x_noncs_var,
w_noncs_var = w_noncs_var,
noncs_max_abs_cor = noncs_max_abs_cor,
suff_block_size = suff_block_size,
verbose = verbose, n_threads = n_threads,
returnModel = returnModel
)
}
