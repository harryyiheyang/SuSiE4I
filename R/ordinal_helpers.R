ocat_is_family <- function(family) {
  inherits(family, "family") &&
    grepl("ordered categorical", tolower(paste(family$family, collapse = " ")),
          fixed = TRUE)
}

ocat_family_ncat <- function(family) {
  if (!ocat_is_family(family)) return(NULL)
  if (!is.function(family$getTheta)) return(NULL)
  alpha <- family$getTheta(TRUE)
  if (is.null(alpha) || !length(alpha)) return(NULL)
  length(alpha) + 1L
}

ocat_prepare_response <- function(y, family = NULL) {
  ncat <- ocat_family_ncat(family)
  if (is.factor(y)) {
    if (is.null(ncat)) ncat <- nlevels(y)
    if (nlevels(y) != ncat) {
      stop("Number of y levels does not match the ordered-categorical family.")
    }
    y_ord <- ordered(y, levels = levels(y))
    y_int <- as.integer(y_ord)
  } else {
    y_num <- as.numeric(y)
    if (any(!is.finite(y_num))) stop("Ordinal y must contain finite class labels.")
    if (any(abs(y_num - round(y_num)) > sqrt(.Machine$double.eps))) {
      stop("Ordinal y must be integer class labels or an ordered factor.")
    }
    y_int <- as.integer(y_num)
    if (is.null(ncat)) ncat <- max(y_int)
    if (any(y_int < 1L) || any(y_int > ncat)) {
      stop("Ordinal y class labels must be in 1, ..., R.")
    }
    y_ord <- ordered(y_int, levels = seq_len(ncat))
  }
  if (ncat < 2L) stop("Ordinal y must have at least two categories.")
  list(y = y_ord, y_int = y_int, ncat = ncat)
}

ocat_backtick <- function(x) {
  paste0("`", gsub("`", "``", x, fixed = TRUE), "`")
}

ocat_formula <- function(rhs, offset = NULL) {
  terms <- character(0)
  if (length(rhs)) terms <- c(terms, ocat_backtick(rhs))
  if (!is.null(offset)) terms <- c(terms, paste0("offset(", offset, ")"))
  if (!length(terms)) return(stats::as.formula("y ~ 1"))
  stats::as.formula(paste("y ~", paste(terms, collapse = " + ")))
}

ocat_predictor_data <- function(Z = NULL, Xextra = NULL, n = NULL) {
  if (is.null(n)) {
    n <- if (!is.null(Z)) nrow(Z) else nrow(Xextra)
  }
  out <- data.frame(row.names = seq_len(n))
  if (!is.null(Z) && ncol(as.matrix(Z)) > 0L) {
    Zdf <- as.data.frame(Z)
    out <- cbind(out, Zdf)
  }
  if (!is.null(Xextra) && ncol(as.matrix(Xextra)) > 0L) {
    Xdf <- as.data.frame(Xextra)
    out <- cbind(out, Xdf)
  }
  out
}

ocat_validate_link <- function(clm_link) {
  match.arg(clm_link, c("logit", "probit", "cloglog", "loglog", "cauchit"))
}

ocat_link_parts <- function(x, clm_link) {
  clm_link <- ocat_validate_link(clm_link)
  if (identical(clm_link, "logit")) {
    F <- stats::plogis(x)
    f <- stats::dlogis(x)
    fp <- f * (1 - 2 * F)
  } else if (identical(clm_link, "probit")) {
    F <- stats::pnorm(x)
    f <- stats::dnorm(x)
    fp <- -x * f
  } else if (identical(clm_link, "cloglog")) {
    ex <- exp(pmin(x, 700))
    F <- -expm1(-ex)
    f <- exp(x - ex)
    fp <- f * (1 - ex)
  } else if (identical(clm_link, "loglog")) {
    emx <- exp(pmin(-x, 700))
    F <- exp(-emx)
    f <- exp(-x - emx)
    fp <- f * (emx - 1)
  } else {
    F <- stats::pcauchy(x)
    f <- stats::dcauchy(x)
    fp <- f * (-2 * x / (1 + x^2))
  }
  left <- x == -Inf
  right <- x == Inf
  if (any(left)) {
    F[left] <- 0
    f[left] <- 0
    fp[left] <- 0
  }
  if (any(right)) {
    F[right] <- 1
    f[right] <- 0
    fp[right] <- 0
  }
  list(F = F, f = f, fp = fp)
}

ocat_fit_explicit <- function(y, pred, alpha_start = NULL, offset = NULL,
                              clm_link, family = NULL,
                              penalty_V = numeric(0)) {
  clm_link <- ocat_validate_link(clm_link)
  if (identical(clm_link, "logit")) {
    stop("clm_logit must use the mgcv::ocat runner.")
  }
  if (!is.null(family) || length(penalty_V)) {
    stop("CLM refits do not accept an mgcv family or a ridge penalty.")
  }
  dat <- cbind(data.frame(y = y), pred)
  rhs <- colnames(pred)
  offset_name <- NULL
  if (!is.null(offset)) {
    offset_name <- "offset_eta"
    dat[[offset_name]] <- as.numeric(offset)
  }
  start <- NULL
  if (!is.null(alpha_start)) {
    start <- c(alpha_start, rep(0, length(rhs)))
  }
  ordinal::clm(
    ocat_formula(rhs, offset = offset_name), data = dat, link = clm_link,
    threshold = "flexible", start = start, model = TRUE
  )
}

ocat_penalty_variance <- function(clm_link, fitX, fitW, term_names,
                                  fitZ = NULL) {
  stats::setNames(numeric(0), character(0))
}

ocat_linear_predictor <- function(fit, pred, n, offset = NULL) {
  eta <- if (is.null(offset)) rep(0, n) else as.numeric(offset)
  if (is.null(fit$beta) || !length(fit$beta)) return(eta)
  beta <- as.numeric(fit$beta)
  names(beta) <- names(fit$beta)
  Xp <- as.matrix(pred[, names(beta), drop = FALSE])
  eta + as.numeric(CppMatrix::matrixVectorMultiply(Xp, beta))
}

ocat_thresholds <- function(fit) {
  alpha <- fit$alpha
  if (is.null(alpha) || any(!is.finite(alpha))) {
    stop("Could not extract finite ordered-categorical thresholds.")
  }
  as.numeric(alpha)
}

ocat_nuisance_coef <- function(fit, n_keep = 0L) {
  bb <- if (n_keep > 0L && length(fit$beta) >= n_keep) {
    clean_coef(fit$beta[seq_len(n_keep)])
  } else {
    numeric(0)
  }
  c(clean_coef(fit$alpha), bb)
}

ocat_coef_table <- function(fit) {
  G <- summary(fit)$coefficients
  if (is.null(G) || is.null(dim(G))) return(NULL)
  G
}

ocat_prob_parts <- function(y_int, eta, alpha, clm_link = "logit",
                            eps = 1e-12) {
  n <- length(y_int)
  K <- length(alpha)
  cuts <- c(-Inf, alpha, Inf)
  lo <- cuts[y_int]
  hi <- cuts[y_int + 1L]
  tl <- lo - eta
  tu <- hi - eta

  lp_lo <- ocat_link_parts(tl, clm_link = clm_link)
  lp_hi <- ocat_link_parts(tu, clm_link = clm_link)
  Fl <- lp_lo$F
  Fu <- lp_hi$F
  fl <- lp_lo$f
  fu <- lp_hi$f
  fpl <- lp_lo$fp
  fpu <- lp_hi$fp

  pr <- pmax(Fu - Fl, eps)
  A <- fl - fu
  B <- fpu - fpl
  u_eta <- A / pr
  h_eta <- A^2 / pr^2 - B / pr

  D <- matrix(0, n, K)
  C <- matrix(0, n, K)
  E <- matrix(0, n, K)

  ii <- which(y_int <= K)
  if (length(ii)) {
    jj <- y_int[ii]
    D[cbind(ii, jj)] <- fu[ii]
    C[cbind(ii, jj)] <- -fpu[ii]
    E[cbind(ii, jj)] <- fpu[ii]
  }

  ii <- which(y_int > 1L)
  if (length(ii)) {
    jj <- y_int[ii] - 1L
    D[cbind(ii, jj)] <- D[cbind(ii, jj)] - fl[ii]
    C[cbind(ii, jj)] <- C[cbind(ii, jj)] + fpl[ii]
    E[cbind(ii, jj)] <- E[cbind(ii, jj)] - fpl[ii]
  }

  h_eta_alpha <- sweep(D, 1L, A / pr^2, "*") -
    sweep(C, 1L, 1 / pr, "*")

  list(
    u_eta = u_eta, h_eta = h_eta, h_eta_alpha = h_eta_alpha,
    D = D, E = E, pr = pr
  )
}

ocat_suffstat_block <- function(Xblk, y_int, eta, Znui, alpha,
                                clm_link = "logit",
                                n_threads = 1, ridge = 1e-6,
                                block_size = 10000L) {
  Xblk <- as.matrix(Xblk)
  n <- nrow(Xblk)
  p <- ncol(Xblk)
  if (is.null(Znui)) {
    Znui <- matrix(nrow = n, ncol = 0)
  } else {
    Znui <- as.matrix(Znui)
    if (nrow(Znui) != n) stop("nrow(Znui) must equal nrow(Xblk).")
  }
  q <- ncol(Znui)
  K <- length(alpha)

  pp <- ocat_prob_parts(y_int = y_int, eta = eta, alpha = alpha,
                        clm_link = clm_link)
  h <- pmax(as.numeric(pp$h_eta), 1e-8)
  sw <- sqrt(h)

  hZ <- if (q > 0L) Znui * h else NULL
  wc <- weighted_crossprod(
    Xblk, h, cbind(h * eta, pp$u_eta, pp$h_eta_alpha, hZ),
    n_threads = n_threads, block_size = block_size
  )
  HXX <- wc$XtWX
  HXeta_eta <- as.numeric(wc$XtM[, 1L])
  UX <- as.numeric(wc$XtM[, 2L])
  HXT <- wc$XtM[, 2L + seq_len(K), drop = FALSE]
  HTeta_eta <- CppMatrix::matrixMultiply(
    matrix(eta, ncol = 1), pp$h_eta_alpha, transA = TRUE
  )
  UT <- colSums(pp$D / pp$pr)

  Dp <- pp$D / pp$pr
  TTT <- crossprod(Dp) - diag(colSums(pp$E / pp$pr), K)

  if (q > 0L) {
    Zh <- Znui * sw
    HXZ <- wc$XtM[, 2L + K + seq_len(q), drop = FALSE]
    HZZ <- CppMatrix::matrixMultiply(Zh, Zh, transA = TRUE)
    HZT <- CppMatrix::matrixMultiply(Znui, pp$h_eta_alpha, transA = TRUE)
    HZeta_eta <- as.numeric(CppMatrix::matrixMultiply(
      Znui, matrix(h * eta, ncol = 1), transA = TRUE
    ))
    UZ <- as.numeric(CppMatrix::matrixMultiply(
      Znui, matrix(pp$u_eta, ncol = 1), transA = TRUE
    ))
    HNN <- rbind(
      cbind(HZZ, HZT),
      cbind(t(HZT), TTT)
    )
    HXN <- cbind(HXZ, HXT)
    HNeta_eta <- c(HZeta_eta, as.numeric(HTeta_eta))
    UN <- c(UZ, UT)
  } else {
    HNN <- TTT
    HXN <- HXT
    HNeta_eta <- as.numeric(HTeta_eta)
    UN <- UT
  }

  HNN_inv_HNX <- solve_with_ridge(HNN, t(HXN), ridge = ridge)
  HNN_inv_rhs <- solve_with_ridge(
    HNN, matrix(HNeta_eta + UN, ncol = 1), ridge = ridge
  )

  XtX <- HXX - CppMatrix::matrixMultiply(HXN, HNN_inv_HNX)
  Xty <- HXeta_eta + UX -
    as.numeric(CppMatrix::matrixMultiply(HXN, HNN_inv_rhs))
  XtX <- (XtX + t(XtX)) / 2
  XtX_pre_ridge <- XtX
  diag(XtX) <- diag(XtX) + ridge

  dimnames(XtX) <- list(colnames(Xblk), colnames(Xblk))
  dimnames(XtX_pre_ridge) <- dimnames(XtX)
  names(Xty) <- colnames(Xblk)
  list(
    XtX = XtX, XtX_pre_ridge = XtX_pre_ridge,
    Xty = Xty, yty = n - 1, n_eff = n,
    min_pr = min(pp$pr), min_h = min(h),
    med_h = stats::median(h), max_h = max(h)
  )
}
