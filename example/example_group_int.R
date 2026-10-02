pkgload::load_all(".", quiet = TRUE)

# Two haplotype factors in Z (reference coding), SNP genotypes in X.
# Truth: level a3 of HapA interacts with snp10; level a1 of HapA interacts
# with level b2 of HapB. Level a4 is rare and is dropped by min_group_int_obs.
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

fit <- SuSiE4I(
  X = X, Z = Z, y = y, family = "gaussian",
  z_groups = list(HapA = colnames(DA), HapB = colnames(DB)),
  min_group_int_obs = 100,
  L_main = 5, L_int = 5, verbose = FALSE
)

fit$main_discoveries
# One row per level column of each selected group; PostMean/PostSD show
# which levels carry the interaction.
fit$interaction_discoveries
