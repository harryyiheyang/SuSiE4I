refit_penalty_variance <- function(fitX, fitW = NULL, term_names,
                                   fitZ = NULL) {
term_names <- as.character(term_names)
if (!length(term_names)) return(stats::setNames(numeric(0), character(0)))

map_stage <- function(fit, prefix, noncs_name = NULL, terms,
                      noncs_V_default = 2) {
V <- as.numeric(fit$V)
if (!length(V) && any(terms != noncs_name)) {
  stop("The SuSiE fit does not contain prior variances V.")
}
out <- stats::setNames(rep(NA_real_, length(terms)), terms)
is_cs <- grepl(paste0("^", prefix, "[0-9]+$"), terms)
if (any(is_cs)) {
idx <- as.integer(sub(prefix, "", terms[is_cs], fixed = TRUE))
if (any(idx < 1L | idx > length(V))) {
stop("A refit credible-set term does not map to a SuSiE component V.")
}
out[is_cs] <- V[idx]
}
is_noncs <- if (is.null(noncs_name)) rep(FALSE, length(terms)) else
  terms == noncs_name
if (any(is_noncs)) {
positive <- V[is.finite(V) & V > 0]
if (!length(positive)) positive <- noncs_V_default
# Use the least-shrunk fitted component to retain the non-CS rescue term.
out[is_noncs] <- max(positive)
}
out
}

main_terms <- grepl("^Main_CS[0-9]+$", term_names) | term_names == "Main_noncs_res"
int_terms <- grepl("^Int_CS[0-9]+$", term_names) | term_names == "Int_noncs_res"
env_terms <- grepl("^Env_CS[0-9]+$", term_names) | term_names == "Z_noncs_res"
unknown <- !main_terms & !int_terms & !env_terms
if (any(unknown)) {
stop("Unknown penalized refit term: ", paste(term_names[unknown], collapse = ", "), ".")
}
out <- stats::setNames(rep(NA_real_, length(term_names)), term_names)
if (any(main_terms)) {
out[main_terms] <- map_stage(
fitX, "Main_CS", "Main_noncs_res", term_names[main_terms]
)
}
if (any(int_terms)) {
if (is.null(fitW)) stop("Interaction refit terms require fitW.")
out[int_terms] <- map_stage(
fitW, "Int_CS", "Int_noncs_res", term_names[int_terms]
)
}
if (any(env_terms)) {
if (is.null(fitZ)) stop("Environment refit terms require fitZ.")
out[env_terms] <- map_stage(
fitZ, "Env_CS", "Z_noncs_res", term_names[env_terms]
)
}
if (any(!is.finite(out) | out <= 0)) {
stop("Every penalized refit term requires a positive finite SuSiE V.")
}
out
}

refit_penalty_terms <- function(term_names) {
term_names[grepl(
"^(Main_CS[0-9]+|Main_noncs_res|Int_CS[0-9]+|Int_noncs_res|Env_CS[0-9]+|Z_noncs_res)$",
term_names, perl = TRUE
)]
}

projection_penalty_precision <- function(ZI, fitX = NULL, fitW = NULL,
                                         fitZ = NULL) {
if (is.null(ZI)) return(stats::setNames(numeric(0), character(0)))
term_names <- refit_penalty_terms(colnames(as.matrix(ZI)))
if (!length(term_names)) {
return(stats::setNames(numeric(0), character(0)))
}
penalty_V <- refit_penalty_variance(
fitX, fitW, term_names, fitZ = fitZ
)
stats::setNames(1 / as.numeric(penalty_V), names(penalty_V))
}

align_projection_precision <- function(ZI, nuisance_precision) {
q <- if (is.null(ZI)) 0L else ncol(as.matrix(ZI))
if (!q) return(numeric(0))

term_names <- colnames(as.matrix(ZI))
if (is.null(term_names)) {
if (length(nuisance_precision)) {
stop("A penalized nuisance projection design must have column names.")
}
return(numeric(q))
}
required <- refit_penalty_terms(term_names)
precision_names <- names(nuisance_precision)
if ((any(!nzchar(term_names)) || anyDuplicated(term_names)) &&
    (length(nuisance_precision) || length(required))) {
stop("A penalized nuisance projection design must have unique non-empty column names.")
}
if (length(nuisance_precision)) {
if (is.null(precision_names) || any(!nzchar(precision_names)) ||
    anyDuplicated(precision_names)) {
stop("nuisance_precision must be a uniquely named vector.")
}
if (any(!is.finite(nuisance_precision) | nuisance_precision <= 0)) {
stop("nuisance_precision must contain positive finite values.")
}
unknown <- setdiff(precision_names, term_names)
if (length(unknown)) {
stop("nuisance_precision names are absent from the projection design: ",
     paste(unknown, collapse = ", "), ".")
}
if (any(precision_names %in% c("Intercept", "(Intercept)"))) {
stop("The projection intercept cannot be penalized.")
}
}
missing_required <- setdiff(required, precision_names)
if (length(missing_required)) {
stop(
"Every penalized refit nuisance term requires projection precision: ",
paste(missing_required, collapse = ", "), "."
)
}

out <- stats::setNames(numeric(q), term_names)
if (length(nuisance_precision)) {
out[precision_names] <- as.numeric(nuisance_precision)
}
out
}

mgcv_refit_dispersion <- function(fit) {
phi <- as.numeric(fit$sig2)[1L]
if (!length(phi) || !is.finite(phi) || phi <= 0) {
phi <- as.numeric(summary(fit)$dispersion)[1L]
}
if (!length(phi) || !is.finite(phi) || phi <= 0) {
stop("The previous mgcv fit must provide a positive finite dispersion.")
}
phi
}

gaussian_ridge_refit <- function(y, pred, penalty_V, dispersion,
                                 offset = NULL, fixed_point = FALSE,
                                 tol = 1e-10, max_iter = 100L,
                                 n_threads = 1L,
                                 block_size = NULL) {
y <- as.numeric(y)
pred <- as.data.frame(pred)
n <- length(y)
if (nrow(pred) != n) stop("nrow(pred) must equal length(y).")
if (is.null(offset)) offset <- rep(0, n)
offset <- as.numeric(offset)
if (length(offset) != n || any(!is.finite(offset))) {
stop("offset must contain one finite value per observation.")
}
if (!is.numeric(dispersion) || length(dispersion) != 1L ||
    !is.finite(dispersion) || dispersion <= 0) {
stop("dispersion must be a positive finite scalar.")
}

G <- cbind("(Intercept)" = 1, as.matrix(pred))
y0 <- y - offset
ymat <- matrix(y0, ncol = 1L)
K <- blockwise_crossprod(
G, n_threads = n_threads, block_size = block_size
)
g <- as.numeric(blockwise_crossprod(
G, ymat, n_threads = n_threads, block_size = block_size
))
yty <- as.numeric(CppMatrix::matrixMultiply(
ymat, ymat, transA = TRUE
))

precision <- numeric(ncol(G))
names(precision) <- colnames(G)
if (length(penalty_V)) {
if (is.null(names(penalty_V)) ||
    !all(names(penalty_V) %in% colnames(pred))) {
stop("penalty_V names must identify columns of pred.")
}
precision[names(penalty_V)] <- 1 / as.numeric(penalty_V)
}

phi <- as.numeric(dispersion)
n_update <- if (isTRUE(fixed_point)) as.integer(max_iter) else 1L
for (i in seq_len(n_update)) {
P <- diag(phi * precision, nrow = ncol(G), ncol = ncol(G))
A <- K + P
b <- as.numeric(CppMatrix::matrixSolve(A, matrix(g, ncol = 1L)))
AiK <- CppMatrix::matrixSolve(A, K)
edf <- sum(diag(AiK))
Kb <- as.numeric(CppMatrix::matrixVectorMultiply(K, b))
rss <- yty - 2 * sum(b * g) + sum(b * Kb)
phi_new <- rss / (n - edf)
if (!is.finite(phi_new) || phi_new <= 0) {
stop("Gaussian ridge refit produced a non-positive dispersion.")
}
if (!isTRUE(fixed_point) ||
    abs(phi_new - phi) <= tol * max(1, phi, phi_new)) break
phi <- phi_new
}

names(b) <- colnames(G)
fit <- list(
coefficients = b, dispersion = phi_new,
sig2 = phi_new, edf = edf, rss = rss,
penalty_dispersion = phi, fixed_point_iterations = i
)
class(fit) <- "gaussian_ridge_refit"
attr(fit, "refit_penalty") <- list(
V = penalty_V, precision = 1 / penalty_V,
working_precision = phi * (1 / penalty_V),
dispersion = phi, sp = 1
)
fit
}

mgcv_fit_fixed_ridge <- function(response, rhs, data, family, penalty_V,
                                 dispersion = 1, offset = NULL,
                                 mgcv_model = NULL) {
penalty_names <- names(penalty_V)
penalty_V <- as.numeric(penalty_V)
names(penalty_V) <- penalty_names
if (!length(penalty_V)) {
return(mgcv_fit_explicit(response, rhs, data, family, offset = offset,
                         mgcv_model = mgcv_model))
}
if (is.null(penalty_names) || any(!nzchar(penalty_names)) ||
    anyDuplicated(penalty_names)) {
stop("penalty_V must be a uniquely named vector.")
}
if (any(!is.finite(penalty_V) | penalty_V <= 0)) {
stop("penalty_V must contain positive finite SuSiE variances.")
}
if (!is.numeric(dispersion) || length(dispersion) != 1L ||
    !is.finite(dispersion) || dispersion <= 0) {
stop("dispersion must be a positive finite scalar.")
}
if (!all(penalty_names %in% rhs) || !all(penalty_names %in% names(data))) {
stop("Every penalized refit term must occur in both rhs and data.")
}

X_pen <- as.matrix(data[, penalty_names, drop = FALSE])
dat <- data[, setdiff(names(data), penalty_names), drop = FALSE]
dat$X_pen <- I(X_pen)
rhs <- c(setdiff(rhs, penalty_names), "X_pen")
PP <- list(X_pen = list(
diag(dispersion / penalty_V, nrow = length(penalty_V),
     ncol = length(penalty_V)),
sp = 1
))

family <- mgcv_patch_family_environment(family)
fml <- mgcv_explicit_formula(response = response, rhs = rhs, offset = offset)
model <- mgcv_model_name(mgcv_model, nrow(dat))
fit_fun <- if (identical(model, "gam")) mgcv::gam else mgcv::bam
method <- if (identical(model, "gam")) "REML" else "fREML"
fit <- fit_fun(fml, data = dat, family = family, method = method, paraPen = PP)

bundled <- grep("^X_pen", names(stats::coef(fit)))
if (length(bundled) != length(penalty_names)) {
stop("mgcv did not preserve every penalized refit coefficient.")
}
names(fit$coefficients)[bundled] <- penalty_names
attr(fit, "refit_penalty") <- list(
V = stats::setNames(penalty_V, penalty_names),
precision = stats::setNames(1 / penalty_V, penalty_names),
working_precision = stats::setNames(dispersion / penalty_V, penalty_names),
dispersion = dispersion, sp = 1
)
fit
}

cox_fit_fixed_ridge <- function(y, status, data, penalty_V,
                                offset = NULL) {
data <- as.data.frame(data)
penalty_names <- names(penalty_V)
penalty_V <- as.numeric(penalty_V)
names(penalty_V) <- penalty_names
if (length(penalty_V) &&
    (is.null(penalty_names) || any(!nzchar(penalty_names)) ||
     any(!is.finite(penalty_V) | penalty_V <= 0))) {
stop("Cox penalty_V must be a named vector of positive finite variances.")
}
if (!all(penalty_names %in% names(data))) {
stop("Every penalized Cox refit term must occur in the model data.")
}
if (!is.null(offset) && !offset %in% names(data)) {
stop("The Cox offset must occur in the model data.")
}

ordinary <- setdiff(names(data), c(penalty_names, offset))
rhs <- if (length(ordinary)) mgcv_backtick(ordinary) else character(0)
if (length(penalty_names)) {
ridge_rhs <- vapply(seq_along(penalty_names), function(i) {
paste0(
"survival::ridge(", mgcv_backtick(penalty_names[i]),
", theta = ", format(1 / penalty_V[i], digits = 17, scientific = TRUE),
", scale = FALSE)"
)
}, character(1))
rhs <- c(rhs, ridge_rhs)
}
if (!is.null(offset)) rhs <- c(rhs, paste0("offset(", offset, ")"))
rhs_text <- if (length(rhs)) paste(rhs, collapse = " + ") else "1"
data$.time <- as.numeric(y)
data$.status <- as.integer(status)
form <- stats::as.formula(paste(
"survival::Surv(.time, .status) ~", rhs_text
))
fit <- survival::coxph(form, data = data, ties = "breslow")

if (length(penalty_names)) {
penalized <- utils::tail(seq_along(stats::coef(fit)), length(penalty_names))
names(fit$coefficients)[penalized] <- penalty_names
attr(fit, "refit_penalty") <- list(
V = stats::setNames(penalty_V, penalty_names),
precision = stats::setNames(1 / penalty_V, penalty_names),
theta = stats::setNames(1 / penalty_V, penalty_names),
scale = FALSE
)
}
fit
}

cox_coef_table <- function(fit) {
G0 <- summary(fit)$coefficients
if (is.null(G0) || !length(G0)) return(NULL)
est <- G0[, "coef"]
se <- G0[, "se(coef)"]
z <- est / se
G <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
           `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
rownames(G) <- rownames(G0)
G
}
