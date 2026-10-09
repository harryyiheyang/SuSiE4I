#' SuSiE4I with an mgcv GAM null model
#'
#' Fine-maps main effects of `X` and their interactions with the covariates of
#' a GAM null model. The null model is a `gam()`/`bam()` formula over `data`
#' (for example `y ~ s(age, by = sex) + sex`); it is fitted once and its
#' smooths are integrated out of both SuSiE stages by the penalized projection
#' with \eqn{(B^\top W B + S_\lambda/\phi)^{-1}}{(B'WB + S/phi)^-1}, i.e. the
#' GAM's \eqn{V_p/\phi}{Vp/phi}. The two stages project each other's credible
#' sets as well (the main stage `[B, Int_CS]`, the interaction stage
#' `[B, Main_CS]`, credible-set columns with ridge \eqn{1/V}{1/V}); no offsets
#' are used. Otherwise the algorithm is that of the GLM path of [SuSiE4I()].
#'
#' By default `s()` uses the truncated linear basis (`bs = "truncated"`,
#' Ruppert, Wand and Carroll): 1, z and \eqn{(z-\tau_j)_+}{(z - tau_j)_+} at the
#' \eqn{j/k}{j/k} quantiles, \eqn{j = 1, \dots, k-1}{j = 1, ..., k - 1}, with
#' (1, z) unpenalized and the truncated lines iid, i.e. binned regression with
#' continuity. `bs = "AMatern"` gives the adaptive Matern basis of mgcv.taps
#' (null space (1, z), built once and reused); other `bs` choices become
#' `"truncated"` (with a warning). `bs = "re"` (mgcv's iid random effect) and
#' `bs = "rw1"` (one coefficient per level of a factor, first-difference penalty
#' in level order) are kept. A factor `by` becomes a common smooth plus one
#' varying-coefficient smooth per centered contrast; a numeric `by` gives a
#' varying coefficient whose constant column is dropped. `te()`, `ti()` and
#' `t2()` are not supported. All covariates (factor contrasts included) and `X`
#' are centered.
#'
#' Interactions follow strong heredity: they are built from main-effect
#' credible sets only. A linear covariate `z` gives the single standardized
#' candidate `Main_CS * z`. A group is one single effect (prior
#' \eqn{V I/d}{V*I/d}, an `lfsr` per column): the centered contrasts of a
#' `bs = "re"` / `"rw1"` factor times `Main_CS`, and, for a smoothed `z`
#' (whichever main-effect basis), `Main_CS * M(z)`, where `M(z)` is the adaptive
#' Matern basis with the intercept as its only null space, whitened by its
#' penalty \eqn{\Omega}{Omega} so that the prior is \eqn{Vc\Omega^{-1}}{V*c*Omega^-1}.
#' In the joint refit a group credible set enters by level as that block with
#' the fixed penalty \eqn{\phi I/(cV)}{phi*I/(c*V)} (`c` scales the block to one
#' unit-variance column) and gets the Wald P value of the block.
#'
#' @param formula Null-model formula with univariate `s()` terms.
#' @param data Data frame with the response and the null-model covariates.
#' @param X Numeric n by p genotype matrix, a `geno` object, or a list of
#'   arguments to `geno_open()` (`bedfile` or `pgenfile`, and optionally
#'   `snp_vec`, `sample_vec`, `impute`); the rows of `data` must follow
#'   `sample_vec` (or the .fam/.psam order).
#' @param family GLM or mgcv family (default `gaussian()`).
#' @param mgcv_model `"gam"` (REML) or `"bam"` (fREML, `discrete = TRUE`).
#' @param k Basis dimension of each smooth unless set inside `s()`.
#' @param noint_vars Names of formula variables that do not interact with `X`.
#' @inheritParams SuSiE4I
#' @return As [SuSiE4I()], plus the null fit `fitNull`.
#' @export
SuSiE4I_GAM <- function(formula, data, X, family = gaussian(), mgcv_model = "gam", k = 10L,
                        scale_data = TRUE, n_threads = 4, L_main = 10, L_int = 5, noint_vars = NULL,
                        coverage_nonkilled = NULL, include_x_squared = FALSE,
                        susie_para_main = NULL, susie_para_int = NULL,
                        max_iter = 10, max_eps = 1e-5, min_iter = 2,
                        x_noncs_var = 0.1, w_noncs_var = 0.1, noncs_max_abs_cor = 0.9,
                        suff_block_size = 10000L, verbose = TRUE, returnModel = FALSE) {
if (is.list(X) && !inherits(X, "geno") && !is.data.frame(X)) {
X <- do.call(geno_open, c(X, list(scale = TRUE, threads = n_threads)))
}
if (inherits(X, "geno")) {
if (!scale_data) X$sd[] <- 1   # centre only, as sweep(X, 2, colMeans(X))
} else {
X <- as.matrix(X)
X <- if (scale_data) large_scale(X) else sweep(X, 2L, colMeans(X))
}
if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
null <- gam_null_setup(formula, data, family, mgcv_model, k, noint_vars)
Run_GAM(X = X, null = null, family = family, mgcv_model = mgcv_model, Lmain = L_main, Lint = L_int,
        max.iter = max_iter, min.iter = min_iter, max.eps = max_eps,
        susie_para_main = .resolve_susie_para(susie_para_main, "susie_para_main"),
        susie_para_int = .resolve_susie_para(susie_para_int, "susie_para_int"),
        verbose = verbose, n_threads = n_threads, x_noncs_var = x_noncs_var, w_noncs_var = w_noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor, include_x_squared = include_x_squared,
        suff_block_size = suff_block_size, coverage_nonkilled = coverage_nonkilled,
        returnModel = returnModel)
}

# Rewrite the null formula (s() to the truncated linear basis, or to the frozen
# AMatern basis when asked), center the covariates, fit the null model once, and
# build the interaction covariates (centered z; M(z) for a smoothed z).
gam_null_setup <- function(formula, data, family, mgcv_model, k, noint_vars) {
ig <- mgcv::interpret.gam(formula)
specs <- ig$smooth.spec
if (any(vapply(specs, inherits, TRUE, what = c("tensor.smooth.spec", "t2.smooth.spec")))) {
stop("te(), ti() and t2() are not supported; use univariate s() terms.")
}
keep <- vapply(specs, inherits, TRUE, what = c("re.smooth.spec", "rw1.smooth.spec"))
am <- vapply(specs, inherits, TRUE, what = "AMatern.smooth.spec")
if (!all(keep | am | vapply(specs, inherits, TRUE, what = "truncated.smooth.spec"))) {
warning("Smooths other than bs = \"truncated\", \"AMatern\", \"re\" and \"rw1\" use the truncated linear basis.", call. = FALSE)
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
if (is.numeric(data[[v]]) && !(v %in% vapply(specs[keep], function(sp) sp$term, ""))) {
new[[v]] <- data[[v]] - mean(data[[v]])
cols[[v]] <- v
} else {
f <- as.factor(data[[v]])
D <- stats::model.matrix(~ f)[, -1L, drop = FALSE]
colnames(D) <- make.names(paste0(v, levels(f)[-1L]), unique = TRUE)
for (j in colnames(D)) new[[j]] <- D[, j] - mean(D[, j])
new[[v]] <- f
cols[[v]] <- colnames(D)
}
}

env <- new.env(parent = environment(formula))
env$.s4i_bases <- list()
sterm <- function(z, by, kz, am) {
if (!am) return(sprintf("s(%s%s, bs = \"truncated\", k = %d)", z, if (is.null(by)) "" else paste0(", by = ", by), kz))
if (is.null(env$.s4i_bases[[z]])) env$.s4i_bases[[z]] <- gam_amatern_base(new[[z]], k = kz)
sprintf("s(%s%s, bs = \"s4iAM\", xt = list(base = .s4i_bases[[\"%s\"]]))", z,
        if (is.null(by)) "" else paste0(", by = ", by), z)
}
rhs <- unlist(cols[par_vars], use.names = FALSE)
for (sp in specs[keep]) rhs <- c(rhs, sprintf("s(%s, bs = \"%s\")", sp$term, if (inherits(sp, "re.smooth.spec")) "re" else "rw1"))
kzs <- list()
for (i in which(!keep)) {
sp <- specs[[i]]
kz <- if (sp$bs.dim > 0) sp$bs.dim else k
kzs[[sp$term]] <- kz
if (identical(sp$by, "NA") || !is.numeric(data[[sp$by]])) rhs <- c(rhs, sterm(sp$term, NULL, kz, am[i]))
if (!identical(sp$by, "NA")) {
rhs <- c(rhs, cols[[sp$by]], vapply(cols[[sp$by]], function(b) sterm(sp$term, b, kz, am[i]), ""))
}
}
fml <- stats::as.formula(paste(ig$response, "~", paste(unique(rhs), collapse = " + ")), env = env)

fit <- if (identical(mgcv_model, "bam")) {
mgcv::bam(fml, data = new, family = family, method = "fREML", discrete = TRUE)
} else mgcv::gam(fml, data = new, family = family, method = "REML")

int_vars <- setdiff(all_vars, noint_vars)
zs <- intersect(names(kzs), int_vars)
Zint <- scale(as.matrix(new[, unlist(cols[setdiff(int_vars, zs)], use.names = FALSE), drop = FALSE]))
# The contrasts of a bs = "re" / "rw1" factor interact as one group single effect.
zgroup <- rep(NA_character_, ncol(Zint))
for (v in intersect(vapply(specs[keep], function(sp) sp$term, ""), int_vars)) zgroup[colnames(Zint) %in% cols[[v]]] <- v
# A smoothed z interacts only through M(z): the AMatern part with the intercept as
# null space, whitened by its penalty Omega so that the group prior V I/d is
# V c Omega^-1 (c = n d / tr(M'M) puts it on the scale of d unit-variance columns).
for (z in zs) {
b1 <- gam_amatern_base(new[[z]], k = kzs[[z]], m = 1L)
e <- eigen(b1$S[-1L, -1L], symmetric = TRUE)
M <- scale(gam_amatern_predict(b1, new[[z]])[, -1L] %*% e$vectors %*% diag(1 / sqrt(e$values)), scale = FALSE)
M <- M / sqrt(mean(M^2))
colnames(M) <- paste0("M(", z, ")", seq_len(ncol(M)))
Zint <- cbind(Zint, M)
zgroup <- c(zgroup, rep(z, ncol(M)))
}
list(fit = fit, formula = fml, data = new, response = ig$response,
     B = stats::predict(fit, type = "lpmatrix"), Zint = Zint, zgroup = zgroup)
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
