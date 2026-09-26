get_pairwise_interactions <- function(W, Z = NULL, noint_env = NULL,
                                      include_x_squared = FALSE) {
if (is.null(W)) return(NULL)
W <- as.matrix(W)
if (!is.logical(include_x_squared) || length(include_x_squared) != 1L ||
    is.na(include_x_squared)) {
stop("include_x_squared must be TRUE or FALSE.")
}
n <- nrow(W); p <- ncol(W)
colnames_W <- colnames(W)
if (is.null(colnames_W)) colnames_W <- paste0("W", seq_len(p))
q <- 0L
Z_mat <- NULL
colnames_Z <- character(0)
if (!is.null(Z)) {
Z0 <- as.matrix(Z)
if (nrow(Z0) != n) stop("nrow(Z) must equal nrow(W).")
if (is.null(noint_env)) noint_env <- integer(0)
noint_env <- intersect(unique(as.integer(noint_env)), seq_len(ncol(Z0)))
indz <- setdiff(seq_len(ncol(Z0)), noint_env)
if (length(indz) > 0L) {
Z_mat <- Z0[, indz, drop = FALSE]
q <- ncol(Z_mat)
colnames_Z <- colnames(Z_mat)
if (is.null(colnames_Z)) colnames_Z <- paste0("Z", seq_len(q))
} else {
Z_mat <- NULL
q <- 0L
}
}
# --- W x W ---
ww_cols <- if (isTRUE(include_x_squared)) {
p * (p + 1L) / 2L
} else {
p * (p - 1L) / 2L
}
WW <- matrix(NA_real_, n, ww_cols)
colnames_WW <- character(ww_cols)
col_idx <- 1L
for (i in seq_len(p)) {
j0 <- if (isTRUE(include_x_squared)) i else i + 1L
if (j0 > p) next
for (j in j0:p) {
WW[, col_idx] <- W[, i] * W[, j]
colnames_WW[col_idx] <- paste0(colnames_W[i], "*", colnames_W[j])
col_idx <- col_idx + 1L
}
}
if (q > 0L) {
zw_cols <- q * p
ZW <- matrix(NA_real_, n, zw_cols)
colnames_ZW <- character(zw_cols)
col_idx <- 1L
for (i in seq_len(q)) {
for (j in seq_len(p)) {
ZW[, col_idx] <- Z_mat[, i] * W[, j]
colnames_ZW[col_idx] <- paste0(colnames_Z[i], "*", colnames_W[j])
col_idx <- col_idx + 1L
}
}
out <- cbind(WW, ZW)
colnames(out) <- c(colnames_WW, colnames_ZW)
} else {
out <- WW
colnames(out) <- colnames_WW
}
out
}

susie_cs_list <- function(fit) {
if (is.null(fit) || is.null(fit$sets$cs) || !length(fit$sets$cs)) {
return(list(index = integer(0), vars = list()))
}
index <- as.integer(fit$sets$cs_index)
vars <- lapply(fit$sets$cs, as.integer)
ord <- order(index)
list(index = index[ord], vars = vars[ord])
}

Identifying_CSEffect <- function(fit, nam, prefix) {
cs <- susie_cs_list(fit)
if (!length(cs$index)) return(NULL)
S <- lapply(seq_along(cs$index), function(k) {
i <- cs$index[k]
a <- cs$vars[[k]]
a <- a[order(fit$pip[a], decreasing = TRUE)]
data.frame(Index = a, Variable = unname(nam[a]), CS = paste0(prefix, i),
           logBF = unname(fit$lbf[i]), PIP = unname(fit$pip[a]))
})
out <- do.call(rbind, S)
rownames(out) <- NULL
out
}

Identifying_MainEffect <- function(fit, nam) {
Identifying_CSEffect(fit, nam, "Main_CS")
}
###############################################################################
Identifying_EnvEffect <- function(fit, nam) {
Identifying_CSEffect(fit, nam, "Env_CS")
}
###############################################################################
Identifying_IntEffect <- function(fitW, namW) {
if (is.null(namW)) return(NULL)
Identifying_CSEffect(fitW, namW, "Int_CS")
}

filter_noncs_interactions <- function(IntIndex) {
if (is.null(IntIndex) || nrow(IntIndex) == 0) return(IntIndex)
if (!("Variable" %in% names(IntIndex))) return(IntIndex)
IntIndex[!grepl("Main_noncs_res", IntIndex$Variable, fixed = TRUE), , drop = FALSE]
}

filter_noncs_joint_coef <- function(G) {
if (is.null(G) || is.null(rownames(G))) return(G)
keep <- !rownames(G) %in% c("Main_noncs_res", "Int_noncs_res")
G[keep, , drop = FALSE]
}

###############################################################################
init_k_from_L <- function(L.init, p) {
if (!is.numeric(L.init) || length(L.init) != 1L ||
    !is.finite(L.init) || L.init < 1) {
stop("L.init must be a positive numeric scalar.")
}
min(as.integer(ceiling(L.init)), as.integer(p))
}

validate_suff_block_size <- function(suff_block_size) {
if (is.null(suff_block_size)) return(NULL)
if (!is.numeric(suff_block_size) || length(suff_block_size) != 1L ||
    !is.finite(suff_block_size) || suff_block_size < 1) {
stop("suff_block_size must be a positive numeric scalar.")
}
as.integer(suff_block_size)
}

residual_variance_from_lm <- function(fit) {
v <- stats::sigma(fit)^2
if (!is.finite(v) || v <= 0) {
r <- stats::residuals(fit)
v <- mean(r^2)
}
if (!is.finite(v) || v <= 0) {
stop("The Gaussian refit produced an invalid residual variance.")
}
as.numeric(v)
}

solve_with_ridge <- function(A, B = NULL, ridge = 1e-8) {
A <- as.matrix(A)
if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
if (is.finite(ridge) && ridge > 0) {
diag(A) <- diag(A) + ridge
}
if (is.null(B)) solve(A) else CppMatrix::matrixSolve(A, as.matrix(B))
}

make_init_data <- function(y = NULL, Z = NULL, X = NULL, selected = integer(0)) {
if (!is.null(y)) {
n <- length(y)
Data <- data.frame(y = y)
} else if (!is.null(Z) && ncol(Z) > 0) {
n <- nrow(Z)
Data <- data.frame(row.names = seq_len(n))
} else {
n <- nrow(X)
Data <- data.frame(row.names = seq_len(n))
}

if (!is.null(Z) && ncol(Z) > 0) {
Zdf <- as.data.frame(Z)
colnames(Zdf) <- paste0("Z", seq_len(ncol(Z)))
Data <- cbind(Data, Zdf)
}

if (length(selected) > 0) {
Xdf <- as.data.frame(X[, selected, drop = FALSE])
colnames(Xdf) <- paste0("InitX", seq_along(selected))
Data <- cbind(Data, Xdf)
}

Data
}

fit_init_lm <- function(X, y, Z, selected) {
Data <- make_init_data(y = y, Z = Z, X = X, selected = selected)
rhs_n <- ncol(Data) - 1L
if (rhs_n == 0L) {
stats::lm(y ~ 1, data = Data, model = FALSE, x = FALSE, y = FALSE)
} else {
stats::lm(y ~ ., data = Data, model = FALSE, x = FALSE, y = FALSE)
}
}

fit_init_cox <- function(X, y, status, Z, selected) {
surv_y <- survival::Surv(y, status)
Data <- make_init_data(Z = Z, X = X, selected = selected)
if (ncol(Data) == 0L) {
survival::coxph(surv_y ~ 1, ties = "breslow", model = TRUE)
} else {
survival::coxph(surv_y ~ ., data = Data, ties = "breslow", model = TRUE)
}
}

select_by_residual_score <- function(X, residual, available) {
if (!any(available)) return(NA_integer_)

r <- as.numeric(residual)
ok <- is.finite(r)
if (!any(ok)) return(NA_integer_)
r[!ok] <- 0
scores <- as.numeric(CppMatrix::matrixMultiply(
X, matrix(r, ncol = 1), transA = TRUE
))
scores[!available] <- NA_real_
scores[!is.finite(scores)] <- NA_real_
if (all(is.na(scores))) return(NA_integer_)
which.max(abs(scores))
}

greedy_lm_warm_start <- function(X, y, Z = NULL, L.init = 1) {
p <- ncol(X)
k_init <- init_k_from_L(L.init, p)
selected <- integer(0)
available <- rep(TRUE, p)
fit <- fit_init_lm(X = X, y = y, Z = Z, selected = selected)

for (step in seq_len(k_init)) {
r <- stats::residuals(fit)
j <- select_by_residual_score(X = X, residual = r, available = available)
if (is.na(j)) break
selected <- c(selected, j)
available[j] <- FALSE
fit <- fit_init_lm(X = X, y = y, Z = Z, selected = selected)
}

list(fit = fit, selected = selected)
}

greedy_cox_warm_start <- function(X, y, status, Z = NULL, L.init = 1) {
p <- ncol(X)
k_init <- init_k_from_L(L.init, p)
selected <- integer(0)
available <- rep(TRUE, p)
fit <- fit_init_cox(X = X, y = y, status = status, Z = Z, selected = selected)

for (step in seq_len(k_init)) {
r <- stats::residuals(fit, type = "martingale")
j <- select_by_residual_score(X = X, residual = r, available = available)
if (is.na(j)) break
selected <- c(selected, j)
available[j] <- FALSE
fit <- fit_init_cox(X = X, y = y, status = status, Z = Z, selected = selected)
}

list(fit = fit, selected = selected)
}

clean_coef <- function(x) {
x <- as.numeric(x)
x[!is.finite(x)] <- 0
x
}

mgcv_backtick <- function(x) {
paste0("`", gsub("`", "``", x, fixed = TRUE), "`")
}

mgcv_explicit_formula <- function(response, rhs, offset = NULL) {
terms <- character(0)
if (length(rhs)) terms <- c(terms, mgcv_backtick(rhs))
if (!is.null(offset)) terms <- c(terms, paste0("offset(", offset, ")"))
if (!length(terms)) return(stats::as.formula(paste(response, "~ 1")))
stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

validate_mgcv_irls_family <- function(family) {
if (!inherits(family, "family")) {
stop("family must be a GLM or mgcv family object.")
}
fam_text <- tolower(paste(c(family$family, class(family)), collapse = " "))
blocked <- c("zero inflated", "zero-inflated", "zip", "ordered",
             "categorical", "cox", "censored")
if (any(vapply(blocked, grepl, logical(1), x = fam_text, fixed = TRUE))) {
stop("Unsupported family for the mgcv IRLS path: ", family$family)
}
if (identical(family$family, "gaussian") && identical(family$link, "log")) {
stop("Lognormal responses are not supported by this IRLS path. Use log(y) with gaussian identity instead.")
}
if (!is.function(family$variance) || !is.function(family$mu.eta)) {
stop("family must provide variance() and mu.eta() for working IRLS.")
}
invisible(TRUE)
}

mgcv_patch_family_environment <- function(family) {
fam_name <- tolower(paste(family$family, collapse = " "))
if (grepl("tweedie", fam_name, fixed = TRUE)) {
ld <- get("ldTweedie", envir = asNamespace("mgcv"))
for (nm in names(family)) {
if (is.function(family[[nm]])) {
env <- environment(family[[nm]])
if (!environmentIsLocked(env) && !exists("ldTweedie", envir = env, inherits = TRUE)) {
assign("ldTweedie", ld, envir = env)
}
}
}
}
family
}

mgcv_model_name <- function(mgcv_model = NULL, n) {
if (is.null(mgcv_model)) {
return("gam")
} else {
if (!is.character(mgcv_model) || length(mgcv_model) != 1L || is.na(mgcv_model)) {
stop("mgcv_model must be NULL, 'gam', or 'bam'.")
}
mgcv_model <- tolower(mgcv_model)
if (!mgcv_model %in% c("gam", "bam")) {
stop("mgcv_model must be NULL, 'gam', or 'bam'.")
}
}
mgcv_model
}

mgcv_fit_explicit <- function(response, rhs, data, family, offset = NULL,
                              mgcv_model = NULL) {
family <- mgcv_patch_family_environment(family)
fml <- mgcv_explicit_formula(response = response, rhs = rhs, offset = offset)
mgcv_model <- mgcv_model_name(mgcv_model, nrow(data))
mgcv_fit <- if (identical(mgcv_model, "gam")) mgcv::gam else mgcv::bam
method <- if (identical(mgcv_model, "gam")) "REML" else "fREML"
mgcv_fit(fml, data = data, family = family, method = method)
}

mgcv_prepare_response <- function(y, family) {
is_binom <- identical(family$family, "binomial") ||
identical(family$family, "quasibinomial")
if (is.matrix(y)) {
if (!is_binom || ncol(y) != 2L) {
stop("Matrix y is only supported for two-column binomial responses.")
}
return(list(
data = data.frame(y_success = y[, 1], y_failure = y[, 2]),
response = "cbind(y_success, y_failure)",
n = nrow(y)
))
}

y <- as.numeric(y)
if (is_binom) {
ymax <- max(y, na.rm = TRUE)
is_count <- all(is.finite(y)) && all(y >= 0) &&
all(abs(y - round(y)) < sqrt(.Machine$double.eps))
if (is_count && ymax > 1) {
trials <- as.integer(ymax)
if (any(y > trials)) stop("Binomial counts exceed inferred trial size.")
return(list(
data = data.frame(y_success = y, y_failure = trials - y),
response = "cbind(y_success, y_failure)",
n = length(y)
))
}
}

list(data = data.frame(y = y), response = "y", n = length(y))
}

mgcv_predictor_data <- function(Z = NULL, Xextra = NULL, n = NULL) {
if (is.null(n)) n <- if (!is.null(Z)) nrow(Z) else nrow(Xextra)
out <- data.frame(row.names = seq_len(n))
if (!is.null(Z) && ncol(Z) > 0L) {
Zdf <- as.data.frame(Z)
colnames(Zdf) <- paste0("Z", seq_len(ncol(Z)))
out <- cbind(out, Zdf)
}
if (!is.null(Xextra) && ncol(as.matrix(Xextra)) > 0L) {
Xdf <- as.data.frame(Xextra)
out <- cbind(out, Xdf)
}
out
}

mgcv_fit_init <- function(X, response_info, Z, selected, family,
                          mgcv_model = NULL) {
Xinit <- NULL
if (length(selected) > 0L) {
Xinit <- X[, selected, drop = FALSE]
colnames(Xinit) <- paste0("InitX", seq_along(selected))
}
pred <- mgcv_predictor_data(Z = Z, Xextra = Xinit, n = response_info$n)
dat <- cbind(response_info$data, pred)
mgcv_fit_explicit(response_info$response, colnames(pred), dat, family,
                  mgcv_model = mgcv_model)
}

mgcv_greedy_warm_start <- function(X, response_info, Z = NULL, family,
                                   L.init = 1,
                                   mgcv_model = NULL) {
p <- ncol(X)
k_init <- init_k_from_L(L.init, p)
selected <- integer(0)
available <- rep(TRUE, p)
fit <- mgcv_fit_init(X = X, response_info = response_info, Z = Z,
                     selected = selected, family = family,
                     mgcv_model = mgcv_model)

for (step in seq_len(k_init)) {
r <- stats::residuals(fit, type = "response")
j <- select_by_residual_score(X = X, residual = r, available = available)
if (is.na(j)) break
selected <- c(selected, j)
available[j] <- FALSE
fit <- mgcv_fit_init(X = X, response_info = response_info, Z = Z,
                     selected = selected, family = family,
                     mgcv_model = mgcv_model)
}

list(fit = fit, selected = selected)
}

robust_weight <- function(w, cutoff = 0.005) {
w <- as.numeric(w)
n <- length(w)
if (!n) return(w)
w[!is.finite(w) | w < 0] <- NA_real_
na <- is.na(w)
n_eff <- sum(!na)
if (n_eff <= 1L) {
w[na] <- 0
return(w)
}
cutoff <- as.numeric(cutoff)[1L]
if (!is.finite(cutoff) || cutoff <= 0) cutoff <- 1e-6
if (cutoff >= 0.05) cutoff <- 0.049
if (n_eff < 1 / cutoff) {
lo <- min(w, na.rm = TRUE)
hi <- max(w, na.rm = TRUE)
} else {
lo <- stats::quantile(w, probs = cutoff, na.rm = TRUE, names = FALSE, type = 7)
hi <- stats::quantile(w, probs = 1 - cutoff, na.rm = TRUE, names = FALSE, type = 7)
}
w <- pmin(pmax(w, lo), hi)
w[na] <- 0
w
}

extract_mgcv_working <- function(fit, weight_cutoff = 0.005,
                                 eta_clip_range = c(-50, 50)) {
eta <- as.numeric(fit$linear.predictors)
eta <- pmin(pmax(eta, eta_clip_range[1]), eta_clip_range[2])
y_work <- as.numeric(fit$y)
family <- fit$family
mu <- as.numeric(family$linkinv(eta))
mu_eta <- family$mu.eta(eta)
var_mu <- family$variance(mu)
prior_w <- fit$prior.weights
if (is.null(prior_w)) prior_w <- rep(1, length(mu))

pseudo_response <- eta + (y_work - mu) / mu_eta
W_diag <- as.numeric(prior_w) * (mu_eta^2) / var_mu
bad <- !is.finite(pseudo_response) | !is.finite(W_diag) | W_diag <= 0
if (mean(bad) > 0.9) stop("Too many invalid working observations.")
if (any(bad)) {
W_diag[bad] <- 0
pseudo_response[bad] <- 0
}
W_diag <- robust_weight(W_diag, cutoff = weight_cutoff)
weight_denom <- sum(W_diag^2)
if (!is.finite(weight_denom) || weight_denom <= 0) {
stop("All working weights are zero.")
}
phi0 <- tryCatch(summary(fit)$dispersion, error = function(e) NA_real_)
if (!is.finite(phi0) || phi0 <= 0) phi0 <- 1

list(
pseudo_response = pseudo_response,
W_diag = W_diag,
phi0 = phi0,
weights = W_diag / phi0,
n_eff = (sum(W_diag)^2) / weight_denom
)
}

safe_add_p <- function(idx, Coefmat) {
if (is.null(idx)) return(NULL)
if (is.data.frame(idx) && nrow(idx) == 0) return(idx)
if (!("CS" %in% names(idx))) return(idx)

cs <- as.character(idx$CS)
pos <- match(cs, rownames(Coefmat))
p   <- rep(NA_real_, length(cs))
p[!is.na(pos)] <- Coefmat[pos[!is.na(pos)], 4]

idx$Pvalue <- p
idx
}

build_noncs_residual <- function(eta, projection_design, cor_design = NULL,
                                 noncs_var = 0.1,
                                 min_eta_var = 1e-7,
                                 noncs_max_abs_cor = 0.9) {
eta <- as.numeric(eta)
if (!length(eta) || any(!is.finite(eta))) return(NULL)
eta_var <- stats::var(eta)
if (!is.finite(eta_var) || eta_var < min_eta_var) return(NULL)

if (is.null(projection_design)) return(NULL)
projection_design <- as.matrix(projection_design)
if (ncol(projection_design) == 0L) return(NULL)
if (nrow(projection_design) != length(eta) || any(!is.finite(projection_design))) return(NULL)

X_full <- cbind(Intercept = 1, projection_design)
if (qr(X_full, tol = 1e-10)$rank < ncol(X_full)) return(NULL)

resid <- as.numeric(stats::lm.fit(x = X_full, y = eta)$residuals)
if (any(!is.finite(resid))) return(NULL)
resid_var_ratio <- stats::var(resid) / eta_var
if (!is.finite(resid_var_ratio) || resid_var_ratio < noncs_var) return(NULL)

cor_mat <- projection_design
if (!is.null(cor_design)) {
cor_design <- as.matrix(cor_design)
if (nrow(cor_design) != length(eta) || any(!is.finite(cor_design))) return(NULL)
cor_mat <- cbind(cor_mat, cor_design)
}
if (ncol(cor_mat) > 0L) {
noncs_max_abs_cor <- min(0.9, noncs_max_abs_cor)
cors <- suppressWarnings(stats::cor(cor_mat, resid))
if (any(is.finite(cors) & abs(cors) >= noncs_max_abs_cor)) return(NULL)
}

resid <- resid - mean(resid)
resid_sd <- stats::sd(resid)
if (!is.finite(resid_sd) || resid_sd <= 0) return(NULL)
as.numeric(resid / resid_sd)
}

build_noncs_refit_term <- function(X, fit, CSdt = NULL, cs_indices = NULL, XCS,
                                    noncs_var = 0.1,
                                    min_eta_var = 1e-7,
                                    noncs_max_abs_cor = 0.9,
                                    cor_design = NULL) {
if (is.null(fit)) return(NULL)

beta_total <- clean_coef(coef.susie(fit)[-1L])
if (!length(beta_total) || length(beta_total) != ncol(X)) return(NULL)

eta_total <- as.numeric(matrixVectorMultiply(X, beta_total))
build_noncs_residual(
eta = eta_total,
projection_design = XCS,
cor_design = cor_design,
noncs_var = noncs_var,
min_eta_var = min_eta_var,
noncs_max_abs_cor = noncs_max_abs_cor
)
}

build_full_noncs_refit_term <- function(X, fit) {
if (is.null(fit)) return(NULL)
if (!length(fit$V)) return(NULL)

beta_total <- clean_coef(coef.susie(fit)[-1L])
if (length(beta_total) != ncol(X)) {
stop("The SuSiE coefficient vector does not match ncol(X).")
}

eta_total <- as.numeric(matrixVectorMultiply(X, beta_total))
eta_total[!is.finite(eta_total)] <- 0

eta_total
}

append_noncs_refit_term <- function(design, term, name,
                                   corr_threshold = 0.9) {
if (is.null(term)) return(design)
term <- as.numeric(term)
if (!length(term) || any(!is.finite(term))) return(design)
if (!is.finite(stats::sd(term)) || stats::sd(term) <= 1e-3) return(design)

design_mat <- as.matrix(design)
ok <- is.finite(term)
if (ncol(design_mat) > 0L && sum(ok) > 2L) {
corr_threshold <- min(0.9, corr_threshold)
cors <- suppressWarnings(stats::cor(design_mat[ok, , drop = FALSE], term[ok]))
if (any(is.finite(cors) & abs(cors) >= corr_threshold)) return(design)
}

out <- cbind(design_mat, term)
colnames(out)[ncol(out)] <- name
out
}

passes_noncs_correlation <- function(term, design, max_abs_cor) {
if (is.null(design)) return(TRUE)
design <- as.matrix(design)
if (ncol(design) == 0L) return(TRUE)
if (nrow(design) != length(term) || any(!is.finite(design))) return(FALSE)

cors <- suppressWarnings(stats::cor(design, term))
max_abs_cor <- min(0.9, max_abs_cor)
!any(is.finite(cors) & abs(cors) >= max_abs_cor)
}

build_w_noncs_refit_term <- function(W, fitW, WCS, etaX, XCS, Z = NULL,
                                     w_noncs_var = 0.1,
                                     min_etaW_var = 1e-7,
                                     noncs_max_abs_cor = 0.9) {
if (is.null(fitW) || is.null(W) || ncol(W) == 0L) return(NULL)

beta_total <- clean_coef(coef.susie(fitW)[-1L])
if (!length(beta_total) || length(beta_total) != ncol(W)) return(NULL)

etaW_total <- as.numeric(matrixVectorMultiply(W, beta_total))
etaW_var <- stats::var(etaW_total)
if (!is.finite(etaW_var) || etaW_var < min_etaW_var) return(NULL)
w_cor_threshold <- min(0.9, noncs_max_abs_cor)
cor_design <- cbind(Z, XCS)

if (is.null(WCS) || ncol(as.matrix(WCS)) == 0L) {
etaX_var <- stats::var(as.numeric(etaX))
if (!is.finite(etaX_var) || etaX_var <= 0) return(NULL)
if (etaW_var / etaX_var < w_noncs_var) return(NULL)
if (!passes_noncs_correlation(etaW_total, cor_design, w_cor_threshold)) return(NULL)
return(etaW_total)
}

build_noncs_residual(
eta = etaW_total,
projection_design = WCS,
cor_design = cor_design,
noncs_var = w_noncs_var,
min_eta_var = min_etaW_var,
noncs_max_abs_cor = w_cor_threshold
)
}

interaction_design_available <- function(W, iter, min_iter,
                                         allow_empty = FALSE) {
if (!is.null(W) && ncol(W) > 0L) return(TRUE)
if (isTRUE(allow_empty)) return(FALSE)
if (iter > min_iter) {
stop(
"No interaction candidates remain after the main-effect CS step. ",
"With one main CS, set include_x_squared = TRUE or provide an interacting environmental covariate.",
call. = FALSE
)
}
FALSE
}

build_cs_design_from_fit <- function(X, fit, prefix) {
if (is.null(fit)) return(list(design = NULL, cs_indices = integer(0)))
cs <- susie_cs_list(fit)
cs_indices <- cs$index
if (!length(cs_indices)) return(list(design = NULL, cs_indices = integer(0)))

Alpha_filtered <- fit$alpha * 0
for (k in seq_along(cs_indices)) {
i <- cs_indices[k]
vars_in_cs_i <- cs$vars[[k]]
vars_in_cs_i <- vars_in_cs_i[vars_in_cs_i >= 1L & vars_in_cs_i <= ncol(X)]
if (length(vars_in_cs_i) > 0) Alpha_filtered[i, vars_in_cs_i] <- fit$alpha[i, vars_in_cs_i] / sum(fit$alpha[i, vars_in_cs_i])
}
Alpha_filtered <- Alpha_filtered * sign(fit$mu)
XCS <- matrixMultiply(X, t(as.matrix(Alpha_filtered)))
XCS <- XCS[, cs_indices, drop = FALSE]
if (is.null(dim(XCS))) XCS <- matrix(XCS, ncol = 1)
colnames(XCS) <- paste0(prefix, cs_indices)
list(design = XCS, cs_indices = cs_indices)
}

build_component_design_from_fit <- function(X, fit, prefix) {
build_cs_design_from_fit(X, fit, prefix)
}

build_select_env_terms <- function(Z, fitZ, noncs_env_var = 0.1,
                                   min_eta_env_var = 1e-7,
                                   noncs_max_abs_cor = 0.9) {
env_cs <- build_cs_design_from_fit(Z, fitZ, "Env_CS")
if (!length(env_cs$cs_indices)) {
return(list(ZCS = NULL, ZCS_refit = NULL, cs_indices = integer(0)))
}

z_noncs <- build_noncs_refit_term(
Z, fitZ, cs_indices = env_cs$cs_indices, XCS = env_cs$design,
noncs_var = noncs_env_var, min_eta_var = min_eta_env_var,
noncs_max_abs_cor = noncs_max_abs_cor
)
ZCS_refit <- append_noncs_refit_term(
env_cs$design, z_noncs, "Z_noncs_res",
corr_threshold = noncs_max_abs_cor
)
list(ZCS = env_cs$design, ZCS_refit = ZCS_refit,
     cs_indices = env_cs$cs_indices)
}

initial_continuous_env <- function(Z, y, Lenv, n_threads = 1,
                                   block_size = NULL) {
y <- as.numeric(y)
if (length(y) != nrow(Z) || any(!is.finite(y))) {
stop("Continuous environment initialization requires one finite outcome per row of Z.",
call. = FALSE)
}
y_sd <- stats::sd(y)
if (!is.finite(y_sd) || y_sd <= 0) {
stop("Continuous environment initialization requires a non-constant outcome.",
call. = FALSE)
}
y <- (y - mean(y)) / y_sd
ZtZ <- blockwise_crossprod(
X = Z, n_threads = n_threads, block_size = block_size
)
Zty <- as.numeric(CppMatrix::matrixMultiply(
Z, matrix(y, ncol = 1), transA = TRUE
))
fitZ <- .fit_susie_stage(
structural = list(
XtX = ZtZ, Xty = Zty, yty = sum(y^2), n = nrow(Z), L = Lenv
),
susie_para = list(), stage = "env", iter = 1, min.iter = 1
)
env_terms <- build_select_env_terms(Z, fitZ)
if (is.null(env_terms$ZCS)) stop_select_env_no_cs()
list(fitZ = fitZ, ZCS = env_terms$ZCS, ZCS_refit = env_terms$ZCS_refit)
}

stop_select_env_no_cs <- function() {
stop("No environmental credible set detected in select_env path; please try setting Z = NULL.",
     call. = FALSE)
}

select_initial_cox_env <- function(Z, y, beta_1se, beta_min) {
Z_index <- which(beta_1se != 0)

if (!length(Z_index)) {
Z_index <- which(beta_min != 0)
}

if (!length(Z_index)) {
z_cor <- abs(as.numeric(stats::cor(Z, y)))
z_cor[!is.finite(z_cor)] <- -Inf
if (all(z_cor == -Inf)) {
stop("Cox initialization failed because all cor(Z, observed_time) values are undefined.",
call. = FALSE)
}
Z_index <- which.max(z_cor)
}

Z_index
}

weighted_projected_suffstats <- function(X, y, ZI, weights,
                                         nuisance_precision,
                                         n_threads = 1,
                                         ridge = 1e-8,
                                         block_size = NULL) {
if (missing(nuisance_precision)) {
stop("nuisance_precision must be supplied explicitly for every projection.")
}
if (!is.null(ZI)) ZI <- as.matrix(ZI)
q <- if (is.null(ZI)) 0L else ncol(ZI)
projection_precision <- align_projection_precision(ZI, nuisance_precision)
if (!is.null(block_size)) block_size <- max(1L, as.integer(block_size))

y <- as.numeric(y)
weights <- as.numeric(weights)
weights[!is.finite(weights) | weights < 0] <- 0

wy <- weights * y
Zw <- if (q > 0L) ZI * weights else NULL
wc <- weighted_crossprod(X, weights, cbind(Zw, wy),
                         n_threads = n_threads, block_size = block_size)
XtX <- wc$XtWX
Xty <- as.numeric(wc$XtM[, q + 1L])
yty <- sum(weights * y^2)
yty_raw <- yty
if (q > 0L) {
ZtZ <- matrixMultiply(ZI, Zw, transA = TRUE)
diag(ZtZ) <- diag(ZtZ) + projection_precision
ZtX <- t(wc$XtM[, seq_len(q), drop = FALSE])
Zty <- as.numeric(matrixMultiply(ZI, matrix(wy, ncol = 1), transA = TRUE))
rm(Zw)

Zinv_ZtX <- solve_with_ridge(ZtZ, ZtX, ridge = ridge)
Zinv_Zty <- solve_with_ridge(ZtZ, matrix(Zty, ncol = 1), ridge = ridge)

XtX <- XtX - matrixMultiply(ZtX, Zinv_ZtX, transA = TRUE)
Xty <- Xty - as.numeric(matrixMultiply(ZtX, Zinv_Zty, transA = TRUE))
yty <- yty - as.numeric(crossprod(Zty, Zinv_Zty))
}

XtX <- (XtX + t(XtX)) / 2
if (is.finite(yty) && yty < 0 &&
    yty > -sqrt(.Machine$double.eps) * max(1, abs(yty_raw))) {
yty <- 0
}
dimnames(XtX) <- list(colnames(X), colnames(X))
names(Xty) <- colnames(X)
list(XtX = XtX, Xty = Xty, yty = as.numeric(yty))
}

.susie_default_para <- function(stage = c("main", "int", "env"),
                                 gaussian = FALSE,
                                 residual_variance = NULL) {
stage <- match.arg(stage)
native_formals <- formals(susieR::susie_ss)
coverage <- eval(native_formals[["coverage"]], envir = asNamespace("susieR"))
min_abs_corr <- eval(native_formals[["min_abs_corr"]],
                     envir = asNamespace("susieR"))
out <- list(
standardize = FALSE,
scaled_prior_variance = 2,
estimate_prior_variance = TRUE,
estimate_prior_method = "optim",
max_iter = 300,
coverage = coverage,
min_abs_corr = min_abs_corr
)

if (isTRUE(gaussian)) {
out$estimate_residual_variance <- FALSE
if (!is.null(residual_variance)) {
out$residual_variance <- as.numeric(residual_variance)
}
} else {
out$estimate_residual_variance <- TRUE
out$residual_variance <- 0.5
out$residual_variance_lowerbound <- 0.1
out$residual_variance_upperbound <- 1.01
}
out
}

.validate_susie_para <- function(x, arg) {
if (is.null(x)) return(list())
if (!is.list(x)) stop(arg, " must be NULL or a named list.")
if (!length(x)) return(list())
nm <- names(x)
if (is.null(nm) || any(!nzchar(nm))) stop(arg, " must be a named list.")
if (anyDuplicated(nm)) stop(arg, " must not contain duplicate names.")

protected <- c("XtX", "Xty", "yty", "n", "L")
blocked <- intersect(nm, protected)
if (length(blocked)) {
stop(arg, " cannot set structural SuSiE inputs: ",
     paste(blocked, collapse = ", "), ".")
}

if (all(c("prior_variance", "scaled_prior_variance") %in% nm)) {
stop(arg, " cannot contain both prior_variance and scaled_prior_variance.")
}
if ("prior_variance" %in% nm) {
V <- x$prior_variance
if (!is.numeric(V) || length(V) != 1L || !is.finite(V) || V <= 0) {
stop(arg, "$prior_variance must be a positive finite numeric scalar.")
}
}
if ("scaled_prior_variance" %in% nm) {
V <- x$scaled_prior_variance
if (!is.numeric(V) || length(V) != 1L || !is.finite(V) || V <= 0) {
stop(arg, "$scaled_prior_variance must be a positive finite numeric scalar.")
}
warning(
"scaled_prior_variance is interpreted as an absolute coefficient prior variance and re-scaled at every IRLS iteration; use prior_variance instead.",
call. = FALSE
)
x$prior_variance <- V
x$scaled_prior_variance <- NULL
nm <- names(x)
}

valid <- c(names(formals(susieR::susie_ss)), "prior_variance")
unknown <- setdiff(nm, valid)
if (length(unknown)) {
stop("Unknown susieR::susie_ss parameter in ", arg, ": ",
     paste(unknown, collapse = ", "), ".")
}
x
}

.resolve_susie_para <- function(susie_para = NULL, arg) {
.validate_susie_para(susie_para, arg)
}

.susie_iteration_args <- function(susie_para, structural, stage,
                                  iter, min.iter, gaussian = FALSE,
                                  residual_variance = NULL) {
args <- .susie_default_para(
stage = stage, gaussian = gaussian,
residual_variance = residual_variance
)
overrides <- susie_para[!vapply(susie_para, is.null, logical(1))]
if (iter <= min.iter && length(overrides)) {
warm_V_controls <- c(
"prior_variance", "scaled_prior_variance", "estimate_prior_variance"
)
overrides <- overrides[setdiff(names(overrides), warm_V_controls)]
}
if (length(overrides)) {
args[names(overrides)] <- overrides
}
if (iter <= min.iter) {
args$estimate_prior_variance <- FALSE
}

epv <- args$estimate_prior_variance
if (!is.logical(epv) || length(epv) != 1L || is.na(epv)) {
stop("estimate_prior_variance must be TRUE or FALSE.")
}
if ("prior_variance" %in% names(args)) {
y_scale <- structural$yty / (structural$n - 1)
if (!is.numeric(y_scale) || length(y_scale) != 1L ||
    !is.finite(y_scale) || y_scale <= 0) {
stop("structural$yty / (structural$n - 1) must be positive and finite to convert prior_variance.")
}
args$scaled_prior_variance <- args$prior_variance / y_scale
args$prior_variance <- NULL
}
args[names(structural)] <- structural
args
}

.fit_susie_stage <- function(structural, susie_para, stage,
                             iter, min.iter, gaussian = FALSE,
                             residual_variance = NULL) {
args <- .susie_iteration_args(
susie_para = susie_para, structural = structural, stage = stage,
iter = iter, min.iter = min.iter, gaussian = gaussian,
residual_variance = residual_variance
)
fit <- do.call(susieR::susie_ss, args)
fit$cs_config <- list(
coverage = as.numeric(args$coverage),
min_abs_corr = as.numeric(args$min_abs_corr)
)
if (is.null(fit$sets$requested_coverage) ||
    !isTRUE(all.equal(as.numeric(fit$sets$requested_coverage),
                     fit$cs_config$coverage))) {
stop("The fitted SuSiE requested coverage does not match its effective CS configuration.")
}
fit
}

make_diagnostics <- function(iterations, eps, start_time) {
final_eps <- if (length(eps)) as.numeric(utils::tail(eps, 1L)) else NA_real_
data.frame(
iterations = as.integer(iterations),
eps = final_eps,
runtime_seconds = unname(proc.time()[["elapsed"]] - start_time)
)
}

extract_direction_table <- function(fit, joint_coef) {
  int_idx <- fit$interaction_discoveries
  if (is.null(int_idx) || nrow(int_idx) == 0) { return(NULL) }
  JC <- filter_noncs_joint_coef(joint_coef)
  results <- vector("list", nrow(int_idx))
  for (i in seq_len(nrow(int_idx))) {
    var_name <- as.character(int_idx$Variable[i])
    int_cs   <- as.character(int_idx$CS[i])
    parts        <- strsplit(var_name, "\\*")[[1]]
    is_gxg       <- all(grepl("^Main_", parts))
    is_quadratic <- is_gxg && (parts[1] == parts[2])
    if (is_quadratic) {
      main_term <- parts[1]
      main_z    <- if (main_term %in% rownames(JC)) JC[main_term, 3] else NA_real_
      int_z     <- if (int_cs    %in% rownames(JC)) JC[int_cs,    3] else NA_real_
      cat_id <- NA_character_
      anno   <- NA_character_
      if (!is.na(main_z) && !is.na(int_z)) {
        if (sign(main_z) != sign(int_z)) { cat_id <- "5"; anno <- "Overdominance"
        } else { cat_id <- "6"; anno <- "Dominance" }
      }
      results[[i]] <- data.frame(Main_Term=main_term, Environment_Term=NA_character_, Interaction_Zscore=int_z, Main_Zscore=main_z, Environment_Zscore=NA_real_, Category=cat_id, Annotation=anno, stringsAsFactors=FALSE)
      next
    }
    if (is_gxg) {
      z1 <- if (parts[1] %in% rownames(JC)) JC[parts[1], 3] else NA_real_
      z2 <- if (parts[2] %in% rownames(JC)) JC[parts[2], 3] else NA_real_
      if (!is.na(z1) && !is.na(z2) && abs(z1) >= abs(z2)) { main_term <- parts[1]; env_term <- parts[2]; main_z <- z1; env_z <- z2
      } else if (!is.na(z2)) { main_term <- parts[2]; env_term <- parts[1]; main_z <- z2; env_z <- z1
      } else { main_term <- parts[1]; env_term <- parts[2]; main_z <- z1; env_z <- z2 }
    } else {
      main_pos  <- grep("^Main_", parts)
      main_term <- parts[main_pos]
      env_term  <- parts[-main_pos]
      main_z    <- if (main_term %in% rownames(JC)) JC[main_term, 3] else NA_real_
      env_z     <- if (env_term  %in% rownames(JC)) JC[env_term,  3] else NA_real_
    }
    int_z <- if (int_cs %in% rownames(JC)) JC[int_cs, 3] else NA_real_
    cat_id <- NA_character_
    anno   <- NA_character_
    if (!is.na(main_z) && !is.na(env_z) && !is.na(int_z)) {
      s_main <- sign(main_z); s_env <- sign(env_z); s_int <- sign(int_z)
      e_str <- if (is_gxg) "G" else "E"
      if (s_main > 0 && s_env > 0 && s_int > 0) { cat_id <- "1"; anno <- sprintf("Harmful %s amplifies Harmful G", e_str) }
      if (s_main < 0 && s_env < 0 && s_int < 0) { cat_id <- "1"; anno <- sprintf("Protective %s amplifies Protective G", e_str) }
      if (s_main > 0 && s_env > 0 && s_int < 0) { cat_id <- "2"; anno <- sprintf("Harmful %s buffers Harmful G", e_str) }
      if (s_main < 0 && s_env < 0 && s_int > 0) { cat_id <- "2"; anno <- sprintf("Protective %s buffers Protective G", e_str) }
      if (s_main > 0 && s_env < 0 && s_int > 0) { cat_id <- "3"; anno <- sprintf("Protective %s eroded by Harmful G", e_str) }
      if (s_main < 0 && s_env > 0 && s_int < 0) { cat_id <- "3"; anno <- sprintf("Harmful %s eroded by Protective G", e_str) }
      if (s_main > 0 && s_env < 0 && s_int < 0) { cat_id <- "4"; anno <- sprintf("Protective %s amplified by Harmful G", e_str) }
      if (s_main < 0 && s_env > 0 && s_int > 0) { cat_id <- "4"; anno <- sprintf("Harmful %s amplified by Protective G", e_str) }
    }
    results[[i]] <- data.frame(Main_Term=main_term, Environment_Term=env_term, Interaction_Zscore=int_z, Main_Zscore=main_z, Environment_Zscore=env_z, Category=cat_id, Annotation=anno, stringsAsFactors=FALSE)
  }
  do.call(rbind, results)
}
