pkgload::load_all(".", quiet = TRUE)

# A drinking-category factor in Z (baseline "never" dropped), SNP genotypes in X.
# Truth: the effect of snp10 differs across drinking levels (spread over levels).
set.seed(11)
n <- 5000L
p <- 200L
X <- sapply(runif(p, 0.1, 0.4), function(f) rbinom(n, 2, f))
colnames(X) <- paste0("snp", seq_len(p))
drink <- sample(c("never", "rare", "monthly", "weekly", "daily"), n, TRUE,
                prob = c(0.3, 0.25, 0.2, 0.15, 0.1))
D <- sapply(c("rare", "monthly", "weekly", "daily"), function(l) as.numeric(drink == l))
colnames(D) <- paste0("dr_", colnames(D))
age <- rnorm(n)
Z <- cbind(D, age = age)

xs <- scale(X)
y <- 0.15 * xs[, 10] - 0.12 * xs[, 50] + 0.1 * D[, "dr_weekly"] + 0.1 * age +
  drop(D %*% c(0.15, -0.1, 0.2, 0.25)) * xs[, 10] + rnorm(n)

fit <- SuSiE4I(
  X = X, Z = Z, y = y, family = "gaussian",
  groupint_ind = list(Drink = colnames(D)),
  L_main = 5, verbose = FALSE
)

fit$main_discoveries
# Each selected group lists its levels with the lfsr of that single effect.
# InCS = FALSE marks refit components without a credible set; Coverage is the CS
# coverage or the purified coverage at coverage_nonkilled.
fit$interaction_discoveries[, c("Pair", "CS", "PIP", "lfsr", "Coverage", "Pvalue", "InCS")]
