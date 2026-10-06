###############################################################################
# Group single-effect regression on sufficient statistics, used by the
# interaction stage when groupint_ind is given.
#
# group[j] is the group id of column j. All columns of a group (e.g. one main
# CS times every level of a factor) form ONE candidate with prior
# beta_g ~ N(0, V Sigma0_g), Sigma0_g = I / d_g (d_g = number of columns) or the
# trace-1 RW1 shape of an ordinal group (attr "Sigma0" of group), so V is the
# total prior variance of the group effect. Levels inside a group do not compete for the
# lbf or alpha; the level is identified by its lfsr. With singleton groups this
# is susie_ss. The fit is returned at column level (alpha and pip of a column
# are those of its group) so the downstream CS / refit code works unchanged.
###############################################################################

gsusie_ss <- function(XtX, Xty, yty, n, L, group, suggested_coverage,
                      scaled_prior_variance = 0.2,
                      estimate_prior_variance = TRUE,
                      residual_variance = NULL,
                      estimate_residual_variance = TRUE,
                      residual_variance_lowerbound = 0,
                      residual_variance_upperbound = Inf,
                      max_iter = 100, tol = 1e-3, coverage = 0.95,
                      min_abs_corr = 0.5, prior_tol = 1e-9, ...) {
XtX <- as.matrix(XtX)
p <- ncol(XtX)
ix <- split(seq_len(p), group)
G <- length(ix)
gid <- match(group, as.integer(names(ix)))
d <- lengths(ix)
S0 <- attr(group, "Sigma0")
# beta_g = Lg u with u ~ N(0, V I), Lg Lg' = Sigma0_g
Lg <- lapply(seq_len(G), function(g) if (is.null(S0)) diag(1 / sqrt(d[g]), d[g]) else t(chol(S0[[g]])))
eg <- lapply(seq_len(G), function(g) eigen(crossprod(Lg[[g]], XtX[ix[[g]], ix[[g]], drop = FALSE] %*% Lg[[g]]), symmetric = TRUE))
lam <- lapply(eg, function(e) pmax(e$values, 0))
s2 <- if (is.null(residual_variance)) yty / (n - 1) else residual_variance
V0 <- scaled_prior_variance * yty / (n - 1)
logsum <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
# per-group log BF under N(0, V Sigma0_g) and the posterior given the group
ser <- function(V, z, post = FALSE) {
out <- numeric(G)
for (g in seq_len(G)) {
a <- 1 / V + lam[[g]] / s2
qz <- drop(crossprod(eg[[g]]$vectors, crossprod(Lg[[g]], z[ix[[g]]]))) / s2
out[g] <- 0.5 * sum(qz^2 / a) - 0.5 * sum(log(a)) - 0.5 * d[g] * log(V)
if (post) {
Q <- Lg[[g]] %*% eg[[g]]$vectors
m[ix[[g]]] <<- drop(Q %*% (qz / a))
v[ix[[g]]] <<- drop((Q^2) %*% (1 / a))
eq[g] <<- sum(lam[[g]] * (qz / a)^2) + sum(lam[[g]] / a)
}
}
out
}
alpha <- matrix(1 / G, L, G)
mu <- matrix(0, L, p)
mu2 <- matrix(0, L, p)
lbf_g <- matrix(0, L, G)
V <- rep(V0, L)
lbf <- rep(0, L)
EQ <- matrix(0, L, G)  # E[b_g' X_g'X_g b_g | g]
B <- matrix(0, L, p)
for (it in seq_len(max_iter)) {
alpha_old <- alpha
for (l in seq_len(L)) {
z <- Xty - drop(XtX %*% colSums(B[-l, , drop = FALSE]))
if (estimate_prior_variance) {
o <- stats::optimize(function(lv) -logsum(ser(exp(lv), z) - log(G)), c(-30, 15))
V[l] <- exp(o$minimum)
}
m <- v <- numeric(p)
eq <- numeric(G)
lb <- ser(V[l], z, post = TRUE)
lbf[l] <- logsum(lb - log(G))
if (estimate_prior_variance && lbf[l] <= 0) {
V[l] <- 0
lb[] <- 0
lbf[l] <- 0
m[] <- 0
v[] <- 0
eq[] <- 0
}
alpha[l, ] <- exp(lb - log(G) - logsum(lb - log(G)))
lbf_g[l, ] <- lb
mu[l, ] <- m
mu2[l, ] <- m^2 + v
EQ[l, ] <- eq
B[l, ] <- alpha[l, gid] * m
}
if (estimate_residual_variance) {
bbar <- colSums(B)
erss <- yty - 2 * sum(bbar * Xty) + sum(bbar * drop(XtX %*% bbar))
for (l in seq_len(L)) {
erss <- erss + sum(alpha[l, ] * EQ[l, ]) - sum(B[l, ] * drop(XtX %*% B[l, ]))
}
s2 <- min(max(erss / n, residual_variance_lowerbound), residual_variance_upperbound)
}
if (max(abs(alpha - alpha_old)) < tol) break
}
# first canonical correlation between two groups (|r| for single columns)
R <- stats::cov2cor(XtX + diag(1e-12, p))
cc <- function(a, b) {
if (a == b) return(1)
ia <- ix[[a]]
ib <- ix[[b]]
ih <- function(M) { e <- eigen(M, symmetric = TRUE); k <- e$values > 1e-10
  e$vectors[, k, drop = FALSE] %*% (t(e$vectors[, k, drop = FALSE]) / sqrt(e$values[k])) }
max(svd(ih(R[ia, ia, drop = FALSE]) %*% R[ia, ib, drop = FALSE] %*%
  ih(R[ib, ib, drop = FALSE]), 0, 0)$d)
}
purity <- function(S) {
if (length(S) < 2L) return(c(1, 1, 1))
r <- utils::combn(S, 2, function(x) cc(x[1], x[2]))
c(min(r), mean(r), stats::median(r))
}
cs <- list()
cs_index <- integer(0)
cs_cov <- numeric(0)
pur <- NULL
for (l in which(V > prior_tol)) {
o <- order(alpha[l, ], decreasing = TRUE)
S <- o[seq_len(which(cumsum(alpha[l, o]) >= coverage)[1L])]
if (any(vapply(cs, function(s) setequal(s, S), TRUE))) next
pu <- purity(S)
if (pu[1] < min_abs_corr) next
cs[[length(cs) + 1L]] <- S
cs_index <- c(cs_index, l)
cs_cov <- c(cs_cov, sum(alpha[l, S]))
pur <- rbind(pur, pu)
}
# non-CS components that reach suggested_coverage after purification
sug <- list(index = integer(0), vars = list(), coverage = numeric(0))
taken <- unlist(cs)
for (l in setdiff(which(V > prior_tol), cs_index)) {
o <- order(alpha[l, ], decreasing = TRUE)
if (o[1L] %in% taken) next
S <- o[seq_len(which(cumsum(alpha[l, o]) >= suggested_coverage)[1L])]
S <- S[!(S %in% taken)]
P <- S[vapply(S, function(g) cc(o[1L], g), 0) >= min_abs_corr]
if (sum(alpha[l, P]) < suggested_coverage) next
sug$index <- c(sug$index, l)
sug$vars <- c(sug$vars, list(unlist(ix[P], use.names = FALSE)))
sug$coverage <- c(sug$coverage, sum(alpha[l, P]))
taken <- c(taken, P)
}
alpha_c <- alpha[, gid, drop = FALSE]
sdv <- sqrt(pmax(mu2 - mu^2, 1e-300))
nrm <- t(apply(mu, 1, function(x) sqrt(rowsum(x^2, gid)[gid, 1])))
fit <- list(
alpha = alpha_c, mu = mu, mu2 = mu2, V = V, lbf = lbf,
lbf_variable = lbf_g[, gid, drop = FALSE], sigma2 = s2,
pip = 1 - apply(1 - alpha_c[V > prior_tol, , drop = FALSE], 2, prod),
lfsr = 1 - alpha_c * pmax(stats::pnorm(mu / sdv), stats::pnorm(-mu / sdv)),
unit_mu = ifelse(nrm > 0, mu / pmax(nrm, 1e-300), 0),
group = gid, null_index = 0, intercept = 0,
X_column_scale_factors = rep(1, p), niter = it,
sets = list(
cs = if (length(cs)) lapply(cs, function(S) unlist(ix[S], use.names = FALSE)) else NULL,
cs_index = if (length(cs)) cs_index else NULL,
coverage = if (length(cs)) cs_cov else NULL,
purity = if (length(cs)) data.frame(min.abs.corr = pur[, 1], mean.abs.corr = pur[, 2],
  median.abs.corr = pur[, 3]) else NULL,
requested_coverage = coverage
),
suggested = sug
)
class(fit) <- "susie"
fit
}
