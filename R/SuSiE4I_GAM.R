#' SuSiE4I with an mgcv GAM null model
#'
#' Fine-maps main effects of `X` and their interactions with the covariates of
#' a GAM null model. The null model is a `gam()`/`bam()` formula over `data`
#' (for example `y ~ s(age, by = sex) + sex`); it is fitted once and its
#' smooths are integrated out of both SuSiE stages by the penalized projection
#' with \eqn{(B^\top W B + S_\lambda/\phi)^{-1}}{(B'WB + S/phi)^-1}, i.e. the
#' GAM's \eqn{V_p/\phi}{Vp/phi}. The main stage uses the interaction fit as an
#' offset and the interaction stage the main fit; only the null model is
#' projected. Otherwise the algorithm is that of the GLM path of [SuSiE4I()].
#'
#' Every `s()` uses the adaptive Matern basis of mgcv.taps whatever `bs` was
#' written (with a warning); it is built once and reused. A factor `by` becomes
#' a common smooth plus one varying-coefficient smooth per centered contrast; a
#' numeric `by` gives a varying coefficient whose constant column is dropped.
#' `te()`, `ti()` and `t2()` are not supported. All covariates (factor
#' contrasts included) and `X` are centered.
#'
#' Interactions follow strong heredity: they are built from main-effect
#' credible sets only, with every covariate `z` (`Main_CS * z`) and, for
#' covariates with a smooth, with its fitted null contribution
#' (`Main_CS * f(z)`, all smooths in `z` including `by` terms).
#'
#' @param formula Null-model formula with univariate `s()` terms.
#' @param data Data frame with the response and the null-model covariates.
#' @param X Numeric n by p genotype matrix.
#' @param family GLM or mgcv family (default `gaussian()`).
#' @param mgcv_model `"gam"` (REML) or `"bam"` (fREML, `discrete = TRUE`).
#' @param k Basis dimension of each smooth unless set inside `s()`.
#' @param noint_vars Names of formula variables that do not interact with `X`.
#' @inheritParams SuSiE4I
#' @return As [SuSiE4I()], plus the null fit `fitNull`.
#' @export
SuSiE4I_GAM <- function(formula, data, X, family = gaussian(), mgcv_model = "gam", k = 10L,
                        scale_data = TRUE, n_threads = 4, L_main = 10, L_int = 5, noint_vars = NULL,
                        int_suggested_coverage = NULL, include_x_squared = FALSE,
                        susie_para_main = NULL, susie_para_int = NULL,
                        max_iter = 10, max_eps = 1e-5, min_iter = 2,
                        x_noncs_var = 0.1, w_noncs_var = 0.1, noncs_max_abs_cor = 0.9,
                        suff_block_size = 10000L, verbose = TRUE, returnModel = FALSE) {
X <- as.matrix(X)
X <- if (scale_data) large_scale(X) else sweep(X, 2L, colMeans(X))
if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
null <- gam_null_setup(formula, data, family, mgcv_model, k, noint_vars)
Run_GAM(X = X, null = null, family = family, mgcv_model = mgcv_model, Lmain = L_main, Lint = L_int,
        max.iter = max_iter, min.iter = min_iter, max.eps = max_eps,
        susie_para_main = .resolve_susie_para(susie_para_main, "susie_para_main"),
        susie_para_int = .resolve_susie_para(susie_para_int, "susie_para_int"),
        verbose = verbose, n_threads = n_threads, x_noncs_var = x_noncs_var, w_noncs_var = w_noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor, include_x_squared = include_x_squared,
        suff_block_size = suff_block_size, int_suggested_coverage = int_suggested_coverage,
        returnModel = returnModel)
}

# Rewrite the null formula (every s() to the frozen AMatern basis), center the
# covariates, fit the null model once, and build the interaction covariates
# (centered z, and f(z) for smoothed z).
gam_null_setup <- function(formula, data, family, mgcv_model, k, noint_vars) {
ig <- mgcv::interpret.gam(formula)
specs <- ig$smooth.spec
if (any(vapply(specs, inherits, TRUE, what = c("tensor.smooth.spec", "t2.smooth.spec")))) {
stop("te(), ti() and t2() are not supported; use univariate s() terms.")
}
if (!all(vapply(specs, inherits, TRUE, what = "AMatern.smooth.spec"))) {
warning("All smooths use the adaptive Matern (AMatern) basis; other bs choices were replaced.", call. = FALSE)
}
par_vars <- attr(stats::terms(ig$pf), "term.labels")
smooth_vars <- unique(vapply(specs, function(sp) sp$term, ""))
by_vars <- setdiff(unique(vapply(specs, function(sp) sp$by, "")), "NA")
all_vars <- unique(c(par_vars, smooth_vars, by_vars))

# Center everything; a factor becomes its centered treatment contrasts.
new <- data.frame(row.names = seq_len(nrow(data)))
new[[ig$response]] <- data[[ig$response]]
cols <- list()
for (v in all_vars) {
if (is.numeric(data[[v]])) {
new[[v]] <- data[[v]] - mean(data[[v]])
cols[[v]] <- v
} else {
f <- as.factor(data[[v]])
D <- stats::model.matrix(~ f)[, -1L, drop = FALSE]
colnames(D) <- make.names(paste0(v, levels(f)[-1L]), unique = TRUE)
for (j in colnames(D)) new[[j]] <- D[, j] - mean(D[, j])
cols[[v]] <- colnames(D)
}
}

env <- new.env(parent = environment(formula))
env$.s4i_bases <- list()
sterm <- function(z, by, kz) {
if (is.null(env$.s4i_bases[[z]])) env$.s4i_bases[[z]] <- gam_amatern_base(new[[z]], k = kz)
sprintf("s(%s%s, bs = \"s4iAM\", xt = list(base = .s4i_bases[[\"%s\"]]))", z,
        if (is.null(by)) "" else paste0(", by = ", by), z)
}
rhs <- unlist(cols[par_vars], use.names = FALSE)
for (sp in specs) {
kz <- if (sp$bs.dim > 0) sp$bs.dim else k
if (identical(sp$by, "NA") || !is.numeric(data[[sp$by]])) rhs <- c(rhs, sterm(sp$term, NULL, kz))
if (!identical(sp$by, "NA")) {
rhs <- c(rhs, cols[[sp$by]], vapply(cols[[sp$by]], function(b) sterm(sp$term, b, kz), ""))
}
}
fml <- stats::as.formula(paste(ig$response, "~", paste(unique(rhs), collapse = " + ")), env = env)

fit <- if (identical(mgcv_model, "bam")) {
mgcv::bam(fml, data = new, family = family, method = "fREML", discrete = TRUE)
} else mgcv::gam(fml, data = new, family = family, method = "REML")

int_vars <- setdiff(all_vars, noint_vars)
Zint <- as.matrix(new[, unlist(cols[int_vars], use.names = FALSE), drop = FALSE])
tt <- stats::predict(fit, type = "terms")
for (z in intersect(smooth_vars, int_vars)) {
lab <- vapply(fit$smooth, function(sm) if (sm$term == z) sm$label else NA_character_, "")
fz <- rowSums(tt[, stats::na.omit(lab), drop = FALSE])
Zint <- cbind(Zint, fz - mean(fz))
colnames(Zint)[ncol(Zint)] <- paste0("f(", z, ")")
}
list(fit = fit, formula = fml, data = new, response = ig$response,
     B = stats::predict(fit, type = "lpmatrix"), Zint = Zint)
}

# S_lambda in the null-model coefficient layout, smoothing parameters taken
# from the latest fit (smooth labels are shared with the joint refit).
gam_null_penalty <- function(fit_null, fit) {
q <- length(stats::coef(fit_null))
S <- matrix(0, q, q)
for (sm in fit_null$smooth) {
ii <- sm$first.para:sm$last.para
S[ii, ii] <- S[ii, ii] + fit$sp[[sm$label]] * sm$S[[1L]]
}
S
}
