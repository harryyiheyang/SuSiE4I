#' SuSiE4I with an mgcv GAM null model
#'
#' Fine-maps main effects of `X` and their interactions with the covariates of
#' a GAM null model. The null model is a `gam()`/`bam()` formula over `data`
#' (for example `y ~ s(age, by = sex) + sex`); it is fitted once and its
#' smooths are integrated out of every SuSiE stage by the penalized projection
#' \eqn{I - B (B^\top W B + S_\lambda/\phi)^{-1} B^\top W}{I - B (B'WB + S/phi)^-1 B'W},
#' the Gaussian random-effect form of the smooth (so the projection matrix is
#' \eqn{V_p/\phi}{Vp/phi} of the GAM). The main stage uses the current
#' interaction fit as an offset and the interaction stage uses the current
#' main fit as an offset; only the null model is projected.
#'
#' Every `s()` term is fitted with the adaptive Matern basis of mgcv.taps
#' (null space `1, z` unpenalized, Matern part projected orthogonal to it),
#' whatever `bs` was written; the basis is built once and reused. `by`
#' variables are supported: a factor `by` becomes a common smooth plus one
#' varying-coefficient smooth per centered contrast, and a numeric `by` gives a
#' varying coefficient without its constant column (its main effect enters as a
#' linear term). `te()`, `ti()` and `t2()` are not supported. All covariates,
#' including factor contrasts, and `X` are centered before fitting.
#'
#' Interactions follow strong heredity: candidates are built only from
#' main-effect credible sets. For every covariate `z` they are `Main_CS * z`
#' (linear modification) and, when `z` has a smooth, `Main_CS * f(z)` with
#' `f(z)` the fitted null-model contribution of all smooths in `z` (including
#' `by` terms). The final joint refit is the null formula plus the
#' main and interaction credible-set terms, the latter ridge-penalized with
#' prior variances from SuSiE.
#'
#' @param formula Null-model formula with response, parametric terms and
#'   univariate `s()` terms.
#' @param data Data frame holding the response and every variable in `formula`.
#' @param X Numeric n by p genotype matrix.
#' @param family A GLM or mgcv family supported by IRLS (default `gaussian()`).
#' @param mgcv_model `"gam"` (REML) or `"bam"` (fREML with `discrete = TRUE`).
#' @param k Basis dimension of every smooth (`k` written inside `s()` wins).
#' @param L_main,L_int Number of SuSiE components for the main and interaction
#'   stages.
#' @param noint_vars Names of formula variables that do not interact with `X`.
#' @param scale_X Whether to standardize `X` (it is always centered).
#' @param susie_para_main,susie_para_int Named `susieR::susie_ss()` options,
#'   as in [SuSiE4I()].
#' @param max_iter,min_iter,max_eps Outer-loop controls, as in [SuSiE4I()].
#' @param n_threads Threads for cross-products.
#' @param suff_block_size Row-block size for cross-products.
#' @param verbose Whether to print iteration diagnostics.
#'
#' @return A list with `fitNull`, `fitX`, `fitW`, the final joint mgcv fit
#'   `fitJoint`, `main_discoveries`, `interaction_discoveries`, `report`,
#'   the rewritten null `formula`, and `diagnostics`.
#' @export
SuSiE4I_GAM <- function(formula, data, X, family = gaussian(),
                        mgcv_model = c("gam", "bam"), k = 10L,
                        L_main = 10, L_int = 5, noint_vars = NULL,
                        scale_X = TRUE,
                        susie_para_main = NULL, susie_para_int = NULL,
                        max_iter = 10, min_iter = 2, max_eps = 1e-5,
                        n_threads = 1, suff_block_size = 10000L,
                        verbose = TRUE) {
run_start <- proc.time()[["elapsed"]]
mgcv_model <- match.arg(mgcv_model)
validate_mgcv_irls_family(family)
family <- mgcv_patch_family_environment(family)
susie_para_main <- .resolve_susie_para(susie_para_main, "susie_para_main")
susie_para_int <- .resolve_susie_para(susie_para_int, "susie_para_int")
suff_block_size <- validate_suff_block_size(suff_block_size)
if (!is.data.frame(data)) stop("data must be a data frame.")
X <- as.matrix(X)
if (!is.numeric(X) || nrow(X) != nrow(data)) {
stop("X must be a numeric matrix with nrow(X) == nrow(data).")
}
if (anyNA(X)) stop("X must not contain missing values.")
if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
X <- if (scale_X) large_scale(X) else sweep(X, 2L, colMeans(X))
n <- nrow(X)
gaussian_id <- identical(family$family, "gaussian") && identical(family$link, "identity")

null <- gam_null_setup(formula, data, k = k, noint_vars = noint_vars)
fit_null <- gam_fit_engine(null$formula, null$data, family, mgcv_model)
B <- stats::predict(fit_null, type = "lpmatrix")
f_hat <- gam_smooth_fits(fit_null, null$smooth_vars)
Zint <- cbind(null$int_linear, f_hat)
Zint <- Zint[, apply(Zint, 2L, stats::var) > 1e-8, drop = FALSE]

fit_final <- fit_null
etaX <- 0
etaW <- 0
beta <- rep(0, ncol(X))
fitX <- NULL
fitW <- NULL
W <- NULL
XCS <- NULL
WCS <- NULL
g <- numeric(0)
for (iter in seq_len(max_iter)) {
beta_prev <- beta
work <- extract_mgcv_working(fit_final, weight_cutoff = 0.0025, eta_clip_range = c(-50, 50))
phi <- work$phi0
P <- gam_null_penalty(fit_null, fit_final) / phi
n_susie <- if (gaussian_id) n else max(0.95 * n, work$n_eff)
res_var <- if (gaussian_id) 1 else work$phi0

ssX <- gam_projected_suffstats(X, work$pseudo_response - etaW, B, work$weights, P,
                               n_threads = n_threads, block_size = suff_block_size)
fitX <- .fit_susie_stage(structural = list(XtX = ssX$XtX, Xty = ssX$Xty, yty = ssX$yty, n = n_susie, L = L_main),
                         susie_para = susie_para_main, stage = "main", iter = iter, min.iter = min_iter,
                         gaussian = gaussian_id, residual_variance = res_var)
beta <- colSums(fitX$alpha * fitX$mu)
XCS <- build_cs_design_from_fit(X, fitX, "Main_CS")$design

W <- NULL
fitW <- NULL
WCS <- NULL
if (!is.null(XCS)) {
# Offset for the interaction stage: main effects from a refit without interactions.
fit_main <- gam_refit(null, XCS, NULL, fitX, NULL, family, mgcv_model, phi)
etaX <- gam_block_eta(fit_main, XCS)
W <- gam_interactions(XCS, Zint)
if (!is.null(W)) {
ssW <- gam_projected_suffstats(W, work$pseudo_response - etaX, B, work$weights, P,
                               n_threads = n_threads, block_size = suff_block_size)
fitW <- .fit_susie_stage(structural = list(XtX = ssW$XtX, Xty = ssW$Xty, yty = ssW$yty, n = n_susie, L = L_int),
                         susie_para = susie_para_int, stage = "int", iter = iter, min.iter = min_iter,
                         gaussian = gaussian_id, residual_variance = res_var)
WCS <- build_cs_design_from_fit(W, fitW, "Int_CS")$design
}
}
fit_final <- if (is.null(XCS)) fit_null else gam_refit(null, XCS, WCS, fitX, fitW, family, mgcv_model, phi)
etaX <- gam_block_eta(fit_final, XCS)
etaW <- gam_block_eta(fit_final, WCS)

err <- sqrt(mean((beta - beta_prev)^2))
phi_new <- mgcv_refit_dispersion(fit_final)
err <- max(err, abs(phi_new - phi) / max(1, phi_new, phi))
g[iter] <- err
if (verbose) cat(sprintf("Iteration %d: err = %.3e\n", iter, err))
if (err < max_eps && iter > min_iter) {
if (verbose) cat("Converged!\n")
break
}
}

G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
MainIndex <- safe_add_p(Identifying_MainEffect(fitX, colnames(X)), G)
IntIndex <- if (is.null(W)) NULL else safe_add_p(Identifying_IntEffect(fitW, colnames(W)), G)
out <- list(diagnostics = make_diagnostics(length(g), g, run_start),
            fitNull = fit_null, fitX = fitX, fitW = fitW, fitJoint = fit_final,
            main_discoveries = MainIndex, interaction_discoveries = IntIndex,
            formula = null$formula)
out$report <- extract_direction_table(out, G)
out
}

# Parse the null formula, rewrite every s() to the frozen AMatern basis and
# build the centered data. Returns the new formula/data, the variables with
# smooths, and the centered linear interaction columns.
gam_null_setup <- function(formula, data, k = 10L, noint_vars = NULL) {
ig <- mgcv::interpret.gam(formula)
response <- ig$response
if (is.null(response) || !response %in% names(data)) {
stop("The formula response must be a column of data.")
}
par_vars <- attr(stats::terms(ig$pf), "term.labels")
bad <- par_vars[!par_vars %in% names(data)]
if (length(bad)) {
stop("Parametric terms must be plain columns of data: ", paste(bad, collapse = ", "), ".")
}
specs <- ig$smooth.spec
redirected <- character(0)
for (sp in specs) {
if (inherits(sp, c("tensor.smooth.spec", "t2.smooth.spec"))) {
stop("te(), ti() and t2() are not supported; use univariate s() terms.")
}
if (length(sp$term) != 1L) stop("Only univariate s() terms are supported: ", sp$label, ".")
if (!sp$term %in% names(data)) stop("Smooth variable not in data: ", sp$term, ".")
if (!identical(sp$by, "NA") && !sp$by %in% names(data)) stop("by variable not in data: ", sp$by, ".")
if (!inherits(sp, "AMatern.smooth.spec")) redirected <- c(redirected, sp$label)
}
if (length(redirected)) {
warning("All smooths use the adaptive Matern (AMatern) basis; replaced the basis of ",
        paste(unique(redirected), collapse = ", "), ".", call. = FALSE)
}
smooth_vars <- unique(vapply(specs, function(sp) sp$term, ""))
by_vars <- unique(vapply(specs, function(sp) sp$by, ""))
by_vars <- by_vars[by_vars != "NA"]
all_vars <- unique(c(par_vars, smooth_vars, by_vars))

is_fac <- function(v) is.factor(data[[v]]) || is.character(data[[v]]) || is.logical(data[[v]])
new <- data.frame(row.names = seq_len(nrow(data)))
new[[response]] <- data[[response]]
contrasts <- list()
for (v in all_vars) {
if (is_fac(v)) {
f <- droplevels(as.factor(data[[v]]))
if (nlevels(f) < 2L) stop("Factor ", v, " has fewer than two levels.")
D <- stats::model.matrix(~ f)[, -1L, drop = FALSE]
colnames(D) <- make.names(paste0(v, levels(f)[-1L]), unique = TRUE)
D <- sweep(D, 2L, colMeans(D))
for (j in colnames(D)) new[[j]] <- D[, j]
contrasts[[v]] <- colnames(D)
} else {
z <- as.numeric(data[[v]])
if (anyNA(z)) stop("Variable ", v, " has missing values.")
new[[v]] <- z - mean(z)
contrasts[[v]] <- v
}
}

env <- new.env(parent = environment(formula))
env$.s4i_bases <- list()
smooth_term <- function(z, by = NULL, kz) {
if (is.null(env$.s4i_bases[[z]])) env$.s4i_bases[[z]] <- gam_amatern_base(new[[z]], k = kz)
by_txt <- if (is.null(by)) "" else paste0(", by = ", by)
sprintf("s(%s%s, bs = \"s4iAM\", xt = list(base = .s4i_bases[[\"%s\"]]))", z, by_txt, z)
}
rhs <- unlist(contrasts[par_vars], use.names = FALSE)
for (sp in specs) {
z <- sp$term
kz <- if (is.numeric(sp$bs.dim) && sp$bs.dim > 0) sp$bs.dim else k
if (identical(sp$by, "NA")) {
rhs <- c(rhs, smooth_term(z, NULL, kz))
} else {
# factor by: common smooth + one varying coefficient per centered contrast;
# numeric by: varying coefficient. The by-variable's main effect is kept.
by_cols <- contrasts[[sp$by]]
if (is_fac(sp$by)) rhs <- c(rhs, smooth_term(z, NULL, kz))
rhs <- c(rhs, by_cols, vapply(by_cols, function(b) smooth_term(z, b, kz), ""))
}
}
rhs <- unique(rhs)
fml <- stats::as.formula(paste(response, "~", paste(rhs, collapse = " + ")), env = env)

int_vars <- setdiff(all_vars, noint_vars)
int_cols <- unlist(contrasts[int_vars], use.names = FALSE)
int_linear <- as.matrix(new[, int_cols, drop = FALSE])
list(formula = fml, data = new, response = response, smooth_vars = intersect(smooth_vars, int_vars),
     int_linear = int_linear)
}

gam_fit_engine <- function(formula, data, family, mgcv_model, paraPen = NULL) {
if (identical(mgcv_model, "gam")) {
mgcv::gam(formula, data = data, family = family, method = "REML", paraPen = paraPen)
} else {
mgcv::bam(formula, data = data, family = family, method = "fREML", discrete = TRUE, paraPen = paraPen)
}
}

# Centered fitted contribution of all smooths in each variable (by terms included).
gam_smooth_fits <- function(fit, vars) {
if (!length(vars)) return(NULL)
tt <- stats::predict(fit, type = "terms")
term_var <- vapply(fit$smooth, function(sm) sm$term[1L], "")
labels <- vapply(fit$smooth, function(sm) sm$label, "")
out <- vapply(vars, function(v) {
cols <- intersect(labels[term_var == v], colnames(tt))
f <- rowSums(tt[, cols, drop = FALSE])
f - mean(f)
}, numeric(nrow(tt)))
out <- matrix(out, nrow(tt))
colnames(out) <- paste0("f(", vars, ")")
out
}

# S_lambda in the null-model coefficient layout, with smoothing parameters
# taken from the current fit (the null fit or the latest joint refit).
gam_null_penalty <- function(fit_null, fit_current) {
q <- length(stats::coef(fit_null))
S <- matrix(0, q, q)
sp <- fit_current$sp
for (sm in fit_null$smooth) {
lam <- if (!is.null(sp) && sm$label %in% names(sp)) sp[[sm$label]] else fit_null$sp[[sm$label]]
ii <- sm$first.para:sm$last.para
S[ii, ii] <- S[ii, ii] + lam * sm$S[[1L]]
}
S
}

# Sufficient statistics of X after integrating out the null smooths:
# X'MX, X'My, y'My with M = W - W B (B'WB + P)^{-1} B'W.
gam_projected_suffstats <- function(X, y, B, weights, P, n_threads = 1, block_size = 10000L) {
y <- as.numeric(y)
weights <- as.numeric(weights)
weights[!is.finite(weights) | weights < 0] <- 0
wy <- weights * y
Bw <- B * weights
wc <- weighted_crossprod(X, weights, cbind(Bw, wy), n_threads = n_threads, block_size = block_size)
q <- ncol(B)
BtWB <- crossprod(B, Bw) + P
BtX <- t(wc$XtM[, seq_len(q), drop = FALSE])
Bty <- as.numeric(crossprod(B, wy))
R <- chol((BtWB + t(BtWB)) / 2 + diag(1e-10, q))
U <- backsolve(R, BtX, transpose = TRUE)
u <- backsolve(R, Bty, transpose = TRUE)
XtX <- wc$XtWX - crossprod(U)
XtX <- (XtX + t(XtX)) / 2
dimnames(XtX) <- list(colnames(X), colnames(X))
Xty <- stats::setNames(as.numeric(wc$XtM[, q + 1L]) - as.numeric(crossprod(U, u)), colnames(X))
yty <- max(sum(weights * y^2) - sum(u^2), 0)
list(XtX = XtX, Xty = Xty, yty = yty)
}

gam_interactions <- function(XCS, Zint) {
if (is.null(XCS) || is.null(Zint) || !ncol(Zint)) return(NULL)
W <- do.call(cbind, lapply(colnames(Zint), function(z) XCS * Zint[, z]))
colnames(W) <- as.vector(outer(colnames(XCS), colnames(Zint), paste, sep = "*"))
W[, apply(W, 2L, stats::var) > 1e-8, drop = FALSE]
}

# Null formula + ridge-penalized credible-set terms (prior variances from SuSiE).
gam_refit <- function(null, XCS, WCS, fitX, fitW, family, mgcv_model, dispersion) {
Xp <- cbind(XCS, WCS)
V <- refit_penalty_variance(fitX, fitW, colnames(Xp))
dat <- as.list(null$data)
dat$X_pen <- Xp
fml <- stats::update(null$formula, . ~ . + X_pen)
environment(fml) <- environment(null$formula)
PP <- list(X_pen = list(diag(dispersion / V, nrow = length(V)), sp = 1))
fit <- gam_fit_engine(fml, dat, family, mgcv_model, paraPen = PP)
ii <- grep("^X_pen", names(fit$coefficients))
if (length(ii) != ncol(Xp)) stop("mgcv did not keep every credible-set coefficient.")
names(fit$coefficients)[ii] <- colnames(Xp)
fit
}

gam_block_eta <- function(fit, D) {
if (is.null(D)) return(0)
drop(D %*% stats::coef(fit)[colnames(D)])
}
