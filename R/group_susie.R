###############################################################################
# Group SuSiE on sufficient statistics.
#
# Each single effect selects one predefined group of columns. Within a group
# the effect vector has prior b_g ~ N(0, V I), and the prior inclusion weight
# is 1 / G over the G groups. Singleton groups reduce to the usual SuSiE SER.
# The fit is used for the interaction stage when Z contains factor-coded
# variables (haplotypes) whose level-wise interactions must enter together.
###############################################################################

is_group_fit <- function(fit) {
inherits(fit, "gsusie")
}

susie_stage_coef <- function(fit) {
if (is_group_fit(fit)) {
return(clean_coef(colSums(fit$alpha[, fit$groups, drop = FALSE] * fit$mu)))
}
clean_coef(coef.susie(fit)[-1L])
}

.group_lbf <- function(V, sigma2, ctil2, d, gid, G) {
if (V <= 0) return(numeric(G))
term <- 0.5 * V * ctil2 / (sigma2 * (sigma2 + V * d)) -
  0.5 * log1p(V * d / sigma2)
out <- numeric(G)
s <- rowsum(term, gid, reorder = TRUE)
out[as.integer(rownames(s))] <- s[, 1]
out
}

.log_sum_exp <- function(x) {
m <- max(x)
m + log(sum(exp(x - m)))
}

.group_ser <- function(r, eig, sigma2, V, log_pi, estimate_V) {
G <- length(eig$idx)
ctil <- numeric(length(eig$d))
for (g in seq_len(G)) {
pos <- eig$pos[[g]]
ctil[pos] <- crossprod(eig$U[[g]], r[eig$idx[[g]]])
}
ctil2 <- ctil^2
lbf_model_at <- function(v) {
.log_sum_exp(log_pi + .group_lbf(v, sigma2, ctil2, eig$d, eig$gid, G))
}
if (isTRUE(estimate_V)) {
opt <- stats::optimize(function(lv) lbf_model_at(exp(lv)),
                       interval = c(-30, 15), maximum = TRUE)
V <- if (opt$objective > 0) exp(opt$maximum) else 0
}
p <- length(r)
if (V <= 0) {
return(list(V = 0, alpha = exp(log_pi), lbf = numeric(G), lbf_model = 0,
            mu = numeric(p), mu2 = numeric(p), EbAb = 0))
}
lbf <- .group_lbf(V, sigma2, ctil2, eig$d, eig$gid, G)
lw <- log_pi + lbf
lbf_model <- .log_sum_exp(lw)
alpha <- exp(lw - lbf_model)
mu <- numeric(p)
mu2 <- numeric(p)
EbAb <- 0
for (g in seq_len(G)) {
pos <- eig$pos[[g]]
j <- eig$idx[[g]]
d <- eig$d[pos]
U <- eig$U[[g]]
shrink <- V / (V * d + sigma2)
mtil <- shrink * ctil[pos]
post_var <- V * sigma2 / (V * d + sigma2)
mu[j] <- as.numeric(U %*% mtil)
mu2[j] <- as.numeric((U^2) %*% post_var) + mu[j]^2
EbAb <- EbAb + alpha[g] * (sum(d * post_var) + sum(d * mtil^2))
}
list(V = V, alpha = alpha, lbf = lbf, lbf_model = lbf_model,
     mu = mu, mu2 = mu2, EbAb = EbAb)
}

.group_eigen <- function(XtX, groups) {
idx <- split(seq_along(groups), groups)
G <- length(idx)
U <- vector("list", G)
d <- vector("list", G)
for (g in seq_len(G)) {
e <- eigen(XtX[idx[[g]], idx[[g]], drop = FALSE], symmetric = TRUE)
U[[g]] <- e$vectors
d[[g]] <- pmax(e$values, 0)
}
sizes <- lengths(idx)
ends <- cumsum(sizes)
pos <- lapply(seq_len(G), function(g) (ends[g] - sizes[g] + 1L):ends[g])
list(idx = idx, U = U, d = unlist(d), pos = pos,
     gid = rep(seq_len(G), sizes))
}

.group_credible_sets <- function(alpha, V, XtX, idx, coverage, min_abs_corr,
                                 prior_tol = 1e-9) {
L <- nrow(alpha)
dXtX <- diag(XtX)
R <- XtX / sqrt(outer(pmax(dXtX, 1e-12), pmax(dXtX, 1e-12)))
cs <- list()
cs_index <- integer(0)
purity <- numeric(0)
for (l in seq_len(L)) {
if (!is.finite(V[l]) || V[l] <= prior_tol) next
o <- order(alpha[l, ], decreasing = TRUE)
k <- which(cumsum(alpha[l, o]) >= coverage)[1L]
if (is.na(k)) k <- length(o)
S <- sort(o[seq_len(k)])
if (any(vapply(cs, identical, logical(1), S))) next
pur <- 1
if (length(S) > 1L) {
pair_cor <- utils::combn(S, 2L, function(gh) {
max(abs(R[idx[[gh[1L]]], idx[[gh[2L]]], drop = FALSE]))
})
pur <- min(pair_cor)
}
if (pur < min_abs_corr) next
cs[[length(cs) + 1L]] <- S
cs_index <- c(cs_index, l)
purity <- c(purity, pur)
}
if (length(cs)) names(cs) <- paste0("L", cs_index)
list(cs = if (length(cs)) cs else NULL,
     cs_index = if (length(cs)) cs_index else NULL,
     purity = if (length(cs)) data.frame(min.abs.corr = purity) else NULL,
     coverage = if (length(cs)) rep(coverage, length(cs)) else NULL,
     requested_coverage = coverage)
}

gsusie_ss <- function(XtX, Xty, yty, n, L, groups, group_names = NULL,
                      scaled_prior_variance = 0.2,
                      residual_variance = NULL,
                      estimate_residual_variance = TRUE,
                      estimate_prior_variance = TRUE,
                      residual_variance_lowerbound = 0,
                      residual_variance_upperbound = Inf,
                      coverage = 0.95, min_abs_corr = 0.5,
                      max_iter = 100, tol = 1e-3, prior_tol = 1e-9, ...) {
XtX <- as.matrix(XtX)
Xty <- as.numeric(Xty)
p <- length(Xty)
if (!identical(dim(XtX), c(p, p))) stop("XtX must be p by p with p = length(Xty).")
if (length(groups) != p) stop("groups must have one entry per column.")
group_levels <- unique(groups)
groups <- match(groups, group_levels)
G <- length(group_levels)
if (is.null(group_names)) group_names <- as.character(group_levels)
if (length(group_names) != G) stop("group_names must have one entry per group.")
L <- max(1L, min(as.integer(L), G))
if (is.null(residual_variance_lowerbound)) residual_variance_lowerbound <- 0
if (is.null(residual_variance_upperbound)) residual_variance_upperbound <- Inf

eig <- .group_eigen(XtX, groups)
var_y <- yty / (n - 1)
sigma2 <- if (is.null(residual_variance)) var_y else as.numeric(residual_variance)
V <- rep(scaled_prior_variance * var_y, L)
log_pi <- rep(-log(G), G)
alpha <- matrix(1 / G, L, G)
mu <- matrix(0, L, p)
mu2 <- matrix(0, L, p)
lbf <- numeric(L)
lbf_variable <- matrix(0, L, G)
EbAb <- numeric(L)
bAb <- numeric(L)
KL <- numeric(L)
XtXb <- numeric(p)
elbo <- -Inf
converged <- FALSE

for (iter in seq_len(max_iter)) {
for (l in seq_len(L)) {
bl <- alpha[l, groups] * mu[l, ]
XtXbl <- as.numeric(XtX %*% bl)
r <- Xty - XtXb + XtXbl
ser <- .group_ser(r, eig, sigma2, V[l], log_pi, estimate_prior_variance)
V[l] <- ser$V
alpha[l, ] <- ser$alpha
mu[l, ] <- ser$mu
mu2[l, ] <- ser$mu2
lbf[l] <- ser$lbf_model
lbf_variable[l, ] <- ser$lbf
EbAb[l] <- ser$EbAb
bl_new <- alpha[l, groups] * mu[l, ]
XtXbl_new <- as.numeric(XtX %*% bl_new)
XtXb <- XtXb - XtXbl + XtXbl_new
bAb[l] <- sum(bl_new * XtXbl_new)
KL[l] <- -ser$lbf_model + (sum(bl_new * r) - 0.5 * ser$EbAb) / sigma2
}
bbar <- colSums(alpha[, groups, drop = FALSE] * mu)
ERSS <- yty - 2 * sum(bbar * Xty) + sum(bbar * XtXb) + sum(EbAb - bAb)
ERSS <- max(ERSS, 1e-12)
elbo_new <- -0.5 * n * log(2 * pi * sigma2) - 0.5 * ERSS / sigma2 - sum(KL)
if (isTRUE(estimate_residual_variance)) {
sigma2 <- min(max(ERSS / n, residual_variance_lowerbound, 1e-12),
              residual_variance_upperbound)
}
if (is.finite(elbo) && abs(elbo_new - elbo) < tol) {
elbo <- elbo_new
converged <- TRUE
break
}
elbo <- elbo_new
}

active <- is.finite(V) & V > prior_tol
pip <- if (any(active)) {
1 - apply(1 - alpha[active, , drop = FALSE], 2L, prod)
} else numeric(G)
names(pip) <- group_names
colnames(alpha) <- group_names
colnames(lbf_variable) <- group_names
colnames(mu) <- colnames(mu2) <- colnames(XtX)
sets <- .group_credible_sets(alpha, V, XtX, eig$idx, coverage, min_abs_corr,
                             prior_tol = prior_tol)
structure(list(
alpha = alpha, mu = mu, mu2 = mu2, V = V, lbf = lbf,
lbf_variable = lbf_variable, pip = pip, sets = sets,
sigma2 = sigma2, elbo = elbo, niter = iter, converged = converged,
groups = groups, group_names = group_names, group_index = eig$idx
), class = c("gsusie", "list"))
}

###############################################################################
# Grouped Z handling
###############################################################################

normalize_z_groups <- function(z_groups, Z) {
if (is.null(z_groups)) return(NULL)
q <- ncol(Z)
nm <- colnames(Z)
out <- rep(NA_character_, q)
if (is.list(z_groups)) {
labels <- names(z_groups)
if (is.null(labels)) labels <- rep("", length(z_groups))
labels[!nzchar(labels)] <- paste0("G", which(!nzchar(labels)))
if (anyDuplicated(labels)) stop("z_groups names must be unique.")
for (k in seq_along(z_groups)) {
cols <- z_groups[[k]]
pos <- if (is.character(cols)) match(cols, nm) else as.integer(cols)
if (!length(pos) || any(is.na(pos)) || any(pos < 1L | pos > q)) {
stop("z_groups element '", labels[k], "' names columns not in Z.")
}
if (any(!is.na(out[pos]))) stop("A Z column belongs to more than one z_groups entry.")
out[pos] <- labels[k]
}
} else {
if (length(z_groups) != q) {
stop("z_groups must be a list of Z columns or a vector with one label per Z column.")
}
lab <- as.character(z_groups)
lab[lab %in% c("", "0")] <- NA_character_
out <- lab
}
if (all(is.na(out))) stop("z_groups does not assign any Z column to a group.")
out
}

z_group_support <- function(Z) {
Z <- as.matrix(Z)
out <- matrix(FALSE, nrow(Z), ncol(Z), dimnames = dimnames(Z))
for (j in seq_len(ncol(Z))) {
z <- Z[, j]
tab <- table(round(z, 8))
base <- as.numeric(names(tab)[which.max(tab)])
out[, j] <- abs(z - base) > 1e-8
}
out
}

check_z_group_coding <- function(Z_raw, z_groups) {
for (f in unique(stats::na.omit(z_groups))) {
cols <- which(z_groups == f)
if (length(cols) < 2L) next
s <- rowSums(as.matrix(Z_raw[, cols, drop = FALSE]))
if (stats::sd(s) < 1e-8) {
warning("The Z columns of group '", f, "' sum to a constant and are collinear with the intercept; drop one reference level.",
        call. = FALSE)
}
}
invisible(NULL)
}

get_grouped_pairwise_interactions <- function(W, Z, noint_env = NULL,
                                              include_x_squared = FALSE,
                                              z_groups,
                                              z_support = NULL,
                                              min_group_int_obs = 100L) {
Z <- as.matrix(Z)
n <- nrow(Z)
nmZ <- colnames(Z)
if (is.null(nmZ)) nmZ <- paste0("Z", seq_len(ncol(Z)))
if (length(z_groups) != ncol(Z)) stop("z_groups must have one entry per Z column.")
if (is.null(z_support)) z_support <- z_group_support(Z)
if (!identical(dim(z_support), dim(Z))) stop("z_support must have the same dimensions as Z.")
if (is.null(noint_env)) noint_env <- integer(0)
interacting <- setdiff(seq_len(ncol(Z)), intersect(unique(as.integer(noint_env)), seq_len(ncol(Z))))
ungrouped <- interacting[is.na(z_groups[interacting])]
grouped <- interacting[!is.na(z_groups[interacting])]
factors <- unique(z_groups[grouped])
n_support <- colSums(z_support)

base <- NULL
if (!is.null(W)) {
W <- as.matrix(W)
if (nrow(W) != n) stop("nrow(Z) must equal nrow(W).")
if (is.null(colnames(W))) colnames(W) <- paste0("W", seq_len(ncol(W)))
base <- get_pairwise_interactions(
W, Z = if (length(ungrouped)) Z[, ungrouped, drop = FALSE] else NULL,
include_x_squared = include_x_squared
)
}

blocks <- list()
block_names <- character(0)
add_block <- function(M, name) {
keep <- apply(M, 2L, function(v) {
s <- stats::sd(v)
is.finite(s) && s > 1e-8
})
if (!any(keep)) return(invisible(NULL))
blocks[[length(blocks) + 1L]] <<- M[, keep, drop = FALSE]
block_names[length(block_names) + 1L] <<- name
invisible(NULL)
}

if (!is.null(W)) {
for (f in factors) {
lv <- grouped[z_groups[grouped] == f]
lv <- lv[n_support[lv] >= min_group_int_obs]
if (!length(lv)) next
for (j in seq_len(ncol(W))) {
M <- Z[, lv, drop = FALSE] * W[, j]
colnames(M) <- paste0(nmZ[lv], "*", colnames(W)[j])
add_block(M, paste0(f, "*", colnames(W)[j]))
}
}
}

if (length(factors) > 1L) {
for (a in seq_len(length(factors) - 1L)) {
for (b in (a + 1L):length(factors)) {
lva <- grouped[z_groups[grouped] == factors[a]]
lvb <- grouped[z_groups[grouped] == factors[b]]
cols <- list()
for (i in lva) {
for (k in lvb) {
if (sum(z_support[, i] & z_support[, k]) < min_group_int_obs) next
v <- Z[, i] * Z[, k]
cols[[paste0(nmZ[i], "*", nmZ[k])]] <- v
}
}
if (!length(cols)) next
M <- do.call(cbind, cols)
colnames(M) <- names(cols)
add_block(M, paste0(factors[a], "*", factors[b]))
}
}
}

n_base <- if (is.null(base)) 0L else ncol(base)
if (!n_base && !length(blocks)) return(NULL)
out <- do.call(cbind, c(if (n_base) list(base) else list(), blocks))
groups <- c(seq_len(n_base),
            n_base + rep(seq_along(blocks), vapply(blocks, ncol, integer(1))))
attr(groups, "group_names") <- c(if (n_base) colnames(base) else character(0), block_names)
attr(out, "groups") <- groups
out
}
