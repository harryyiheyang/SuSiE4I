library(SuSiE4I)
library(mgcv)

set.seed(44)
n <- 160L
p <- 10L
X <- scale(matrix(rnorm(n * p), n, p))
Z <- scale(matrix(rnorm(n), n, 1))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- "Z1"
eta <- 0.2 + 1.0 * X[, 2] + 0.6 * Z[, 1] + 0.8 * X[, 2] * Z[, 1]

y_nb <- rnbinom(n, mu = exp(eta), size = 2.75)
fit_nb_gam <- SuSiE4I(
  X = X, Z = Z, y = y_nb, family = mgcv::nb(theta = NULL),
  mgcv_model = "gam", scale_data = FALSE, n_threads = 1,
  L_main = 2, L_int = 2, coverage_main = 0.5, coverage_int = 0.5,
  min_iter = 1, max_iter = 2, susie_iter = 80,
  verbose = FALSE, refit_noncs = TRUE
)
fit_nb_bam <- SuSiE4I(
  X = X, Z = Z, y = y_nb, family = mgcv::nb(theta = 2.75),
  mgcv_model = "bam", scale_data = FALSE, n_threads = 1,
  L_main = 2, L_int = 2, coverage_main = 0.5, coverage_int = 0.5,
  min_iter = 1, max_iter = 2, susie_iter = 80,
  verbose = FALSE, refit_noncs = TRUE
)

y_tw <- mgcv::rTweedie(mu = exp(eta), p = 1.6, phi = 0.8)
fit_tw_gam <- SuSiE4I(
  X = X, Z = Z, y = y_tw,
  family = mgcv::tw(theta = NULL, link = "log"),
  mgcv_model = "gam", scale_data = FALSE, n_threads = 1,
  L_main = 2, L_int = 2, coverage_main = 0.5, coverage_int = 0.5,
  min_iter = 1, max_iter = 2, susie_iter = 80,
  verbose = FALSE, refit_noncs = TRUE
)
fit_tw_bam <- SuSiE4I(
  X = X, Z = Z, y = y_tw,
  family = mgcv::tw(theta = 1.6, link = "log"),
  mgcv_model = "bam", scale_data = FALSE, n_threads = 1,
  L_main = 2, L_int = 2, coverage_main = 0.5, coverage_int = 0.5,
  min_iter = 1, max_iter = 2, susie_iter = 80,
  verbose = FALSE, refit_noncs = TRUE
)

stopifnot(
  inherits(fit_nb_gam$fitJoint, "gam"),
  inherits(fit_nb_bam$fitJoint, "bam"),
  identical(fit_nb_bam$fitJoint$family$getTheta(TRUE), 2.75),
  inherits(fit_tw_gam$fitJoint, "gam"),
  inherits(fit_tw_bam$fitJoint, "bam"),
  identical(fit_tw_bam$fitJoint$family$getTheta(TRUE), 1.6)
)
