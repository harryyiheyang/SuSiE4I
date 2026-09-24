pkgload::load_all(".", quiet = TRUE)

scale_var <- function(x, v) {
  x <- as.numeric(x) - mean(x)
  x * sqrt(v / stats::var(x))
}

set.seed(20260713)
n <- 600L
p <- 20L
q <- 3L

X <- scale(matrix(rnorm(n * p), n, p))
Z <- scale(matrix(rnorm(n * q), n, q))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- paste0("Z", seq_len(q))

etaZ <- scale_var(Z[, 1] - 0.7 * Z[, 2], 0.30)
etaX <- scale_var(X[, 1] + 0.8 * X[, 2], 0.15)
etaI <- scale_var(X[, 1] * Z[, 1], 0.05)
eta <- etaZ + etaX + etaI

signal_variance <- c(
  etaZ = stats::var(etaZ),
  etaX = stats::var(etaX),
  etaI = stats::var(etaI)
)
stopifnot(max(abs(signal_variance - c(0.30, 0.15, 0.05))) < 1e-10)

# Gaussian
y_gaussian <- eta + stats::rnorm(n)
fit_gaussian <- SuSiE4I(
  X = X, Z = Z, y = y_gaussian, family = gaussian(),
  scale_data = FALSE, n_threads = 1,
  L_main = 3, L_int = 2, min_iter = 2, max_iter = 4,
  susie_iter = 100, verbose = FALSE
)

# Binary, exact 1:2 case-control ratio
case <- sample(seq_len(n), n / 3, prob = exp(eta), replace = FALSE)
y_binary <- integer(n)
y_binary[case] <- 1L
case_control <- table(y_binary)
stopifnot(unname(case_control["1"]) * 2 == unname(case_control["0"]))
fit_binary <- SuSiE4I(
  X = X, Z = Z, y = y_binary, family = binomial("logit"),
  scale_data = FALSE, n_threads = 1,
  L_main = 3, L_int = 2, min_iter = 2, max_iter = 4,
  susie_iter = 100, verbose = FALSE
)

# Cox, exact censoring rate = 0.60
hazard <- exp(eta)
t_event <- rexp(n, rate = hazard)
cox_cutoff <- sort(t_event)[as.integer(0.40 * n)]
t_censor <- rep(cox_cutoff, n)
y_cox <- pmin(t_event, t_censor)
status <- as.integer(t_event <= t_censor)
censoring_rate <- mean(status == 0L)
stopifnot(censoring_rate == 0.60)
fit_cox <- SuSiE4I(
  X = X, Z = Z, y = y_cox, status = status, family = "cox",
  scale_data = FALSE, n_threads = 1,
  L_main = 3, L_int = 2, min_iter = 2, max_iter = 4,
  susie_iter = 100, verbose = FALSE
)

# Negative binomial, theta = 5
y_nb <- rnbinom(n, mu = exp(eta), size = 5)
fit_nb <- SuSiE4I(
  X = X, Z = Z, y = y_nb, family = mgcv::nb(theta = 5),
  scale_data = FALSE, n_threads = 1,
  L_main = 3, L_int = 2, min_iter = 2, max_iter = 4,
  susie_iter = 100, verbose = FALSE
)

# Tweedie, power = 1.5
y_tw <- mgcv::rTweedie(mu = exp(eta), p = 1.5, phi = 1)
fit_tw <- SuSiE4I(
  X = X, Z = Z, y = y_tw,
  family = mgcv::tw(theta = 1.5, link = "log"),
  scale_data = FALSE, n_threads = 1,
  L_main = 3, L_int = 2, min_iter = 2, max_iter = 4,
  susie_iter = 100, verbose = FALSE
)

# Ordinal, four equally sized categories
latent <- eta + stats::rlogis(n)
cutoff <- stats::quantile(latent, probs = 1:3 / 4, names = FALSE)
y_ordinal <- cut(
  latent, breaks = c(-Inf, cutoff, Inf), labels = 1:4,
  ordered_result = TRUE
)
ordinal_counts <- table(y_ordinal)
stopifnot(max(ordinal_counts) - min(ordinal_counts) <= 1L)
fit_ordinal <- SuSiE4I(
  X = X, Z = Z, y = y_ordinal, family = "clm_logit",
  scale_data = FALSE, n_threads = 1,
  L_main = 3, L_int = 2, min_iter = 2, max_iter = 4,
  susie_iter = 100, verbose = FALSE
)
