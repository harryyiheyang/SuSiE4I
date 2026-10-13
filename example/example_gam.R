library(SuSiE4I)

set.seed(3)
n <- 3000L
p <- 200L
age <- runif(n, 40, 75)
sex <- factor(sample(c("F", "M"), n, replace = TRUE))
X <- scale(matrix(rnorm(n * p), n, p))
colnames(X) <- paste0("rs", seq_len(p))

u <- ecdf(age)(age)
male <- as.numeric(sex == "M")
f_age <- scale(plogis(6 * (u - 0.4)) * (1 + 0.5 * male))[, 1] * sqrt(0.3)
y <- f_age + 0.3 * male +
  0.1 * (X[, 5] + X[, 105] + X[, 155]) +
  0.06 * X[, 105] * scale(plogis(6 * (u - 0.4)))[, 1] +   # rs105 modified by age
  0.06 * X[, 5] * scale(male)[, 1] +                       # rs5 modified by sex
  rnorm(n, sd = sqrt(0.6))
dat <- data.frame(y = y, age = age, sex = sex)

# Null model: a sex-specific age curve. Any bs is replaced by AMatern (with a warning).
fit <- SuSiE4I_GAM(y ~ s(age, by = sex) + sex, data = dat, X = X, mgcv_model = "bam")
fit$main_discoveries
fit$interaction_discoveries

# Binary outcome
yb <- rbinom(n, 1, plogis(-0.5 + 2 * (y - mean(y))))
fit_b <- SuSiE4I_GAM(yb ~ s(age, by = sex) + sex, data = transform(dat, yb = yb), X = X,
                     family = binomial(), mgcv_model = "bam", verbose = FALSE)
fit_b$main_discoveries
