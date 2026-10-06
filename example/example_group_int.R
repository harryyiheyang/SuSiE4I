pkgload::load_all(".", quiet = TRUE)

# Two haplotype factors in Z (reference coding), SNP genotypes in X.
# Truth: level a3 of HapA interacts with snp10; level a1 of HapA interacts
# with level b2 of HapB.
set.seed(11)
n <- 5000L
p <- 200L
X <- sapply(runif(p, 0.1, 0.4), function(f) rbinom(n, 2, f))
colnames(X) <- paste0("snp", seq_len(p))
hA <- sample(c("a0", "a1", "a2", "a3", "a4"), n, TRUE,
             prob = c(0.35, 0.25, 0.2, 0.18, 0.015))
hB <- sample(c("b0", "b1", "b2"), n, TRUE, prob = c(0.5, 0.3, 0.2))
DA <- sapply(c("a1", "a2", "a3", "a4"), function(l) as.numeric(hA == l))
DB <- sapply(c("b1", "b2"), function(l) as.numeric(hB == l))
colnames(DA) <- paste0("HA_", colnames(DA))
colnames(DB) <- paste0("HB_", colnames(DB))
age <- rnorm(n)
Z <- cbind(DA, DB, age = age)

xs <- scale(X)
y <- 0.15 * xs[, 10] - 0.12 * xs[, 50] + 0.2 * DA[, "HA_a2"] + 0.1 * age +
  0.35 * DA[, "HA_a3"] * xs[, 10] + 0.5 * DA[, "HA_a1"] * DB[, "HB_b2"] +
  rnorm(n)

# HapA and HapB are marked cross, so HapA x HapB interactions are also built.
gi <- list(HapA = colnames(DA), HapB = colnames(DB))
attr(gi$HapA, "cross") <- TRUE
attr(gi$HapB, "cross") <- TRUE
fit <- SuSiE4I(
  X = X, Z = Z, y = y, family = "gaussian", groupint_ind = gi,
  L_main = 5, verbose = FALSE
)

fit$main_discoveries
# Each selected group lists its levels; lfsr identifies the interacting levels.
# InCS = FALSE marks refit components without a credible set; Coverage is the CS
# coverage or the purified coverage at coverage_nonkilled.
fit$interaction_discoveries[, c("Pair", "CS", "PIP", "lfsr", "Coverage", "Pvalue", "InCS")]

# An ordinal factor (levels in order, reference lowest) gets the random-walk prior.
attr(gi$HapA, "type") <- "ordinal"
fit_ord <- SuSiE4I(X = X, Z = Z, y = y, family = "gaussian", groupint_ind = gi,
                   L_main = 5, verbose = FALSE)
fit_ord$interaction_discoveries[, c("Pair", "CS", "PIP", "lfsr", "Pvalue", "InCS")]
