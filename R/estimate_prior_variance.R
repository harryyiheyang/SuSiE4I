#' Estimate the SuSiE prior variance from selected variants
#'
#' Fits the outcome once using `Z` and the all-chromosome PRS, then constructs
#' chromosome-specific LOCO sufficient statistics for the columns of `X`.
#' The columns of `X` must already be standardized.
#'
#' @param X Standardized numeric genotype matrix with SNP IDs as column names.
#' @param y Outcome accepted by [SuSiE4I()]. A Cox outcome may be supplied as a
#'   `survival::Surv` object.
#' @param Z Non-empty numeric covariate matrix that does not contain PRS
#'   columns. `NULL` is not supported because the nuisance predictor must be
#'   fixed before projection.
#' @param PRS Numeric matrix with 22 or 23 chromosome-specific PRS columns, or
#'   `NULL`. Column names must exactly match `SNPInfo$CHR` values used by `X`.
#' @param SNPInfo Data frame containing `SNP` and `CHR`. Its SNPs must match
#'   `colnames(X)` exactly.
#' @param family Outcome family accepted by [SuSiE4I()].
#' @param status Optional event indicator for Cox models when `y` is a time
#'   vector.
#' @param mgcv_model `NULL`, `"gam"`, or `"bam"` for mgcv nuisance fits.
#' @param n_threads Number of threads used for matrix operations.
#' @param suff_block_size Row-block size for sufficient-statistic cross-products.
#' @param r_thres Absolute within-chromosome genotype-correlation threshold.
#'   The default `0` uses the independent-variant estimator. Positive values
#'   use chromosome-block sparse sufficient statistics.
#' @param max_iter Maximum number of EM updates for the prior variance.
#' @param tol Relative EM convergence tolerance.
#'
#' @return A list with `prior_variance`, `scaled_prior_variance`, and a one-row
#'   `diagnostics` data frame. The recommended package use is
#'   `susie_para_main = list(prior_variance = out$prior_variance,
#'   estimate_prior_variance = FALSE)`. The scaled value is retained for direct
#'   `susieR::susie_ss()` use with this estimator's sufficient-statistic
#'   convention. Do not rerun this estimator inside the outer iterations.
#'
#' @export
estimate_prior_variance <- function(X, y, Z, PRS, SNPInfo,
                                    family = NULL, status = NULL,
                                    mgcv_model = NULL,
                                    n_threads = 4,
                                    suff_block_size = 10000L,
                                    r_thres = 0,
                                    max_iter = 100L,
                                    tol = 1e-6) {
  X <- as.matrix(X)
  has_prs <- !is.null(PRS)
  if (has_prs) PRS <- as.matrix(PRS)
  if (!is.numeric(X) || !nrow(X) || !ncol(X)) {
    stop("X must be a non-empty numeric matrix.")
  }
  if (is.null(colnames(X)) || any(!nzchar(colnames(X))) ||
      anyDuplicated(colnames(X))) {
    stop("X must have unique, non-empty SNP column names.")
  }
  n <- nrow(X)
  p <- ncol(X)
  if (has_prs) {
    if (!is.numeric(PRS) || nrow(PRS) != n ||
        !ncol(PRS) %in% c(22L, 23L)) {
      stop("PRS must be NULL or a numeric matrix with nrow(X) rows and 22 or 23 columns.")
    }
    if (is.null(colnames(PRS)) || any(!nzchar(colnames(PRS))) ||
        anyDuplicated(colnames(PRS))) {
      stop("PRS must have unique, non-empty chromosome column names.")
    }
  }
  if (!is.data.frame(SNPInfo) || !all(c("SNP", "CHR") %in% names(SNPInfo))) {
    stop("SNPInfo must be a data frame containing SNP and CHR.")
  }

  snp <- as.character(SNPInfo$SNP)
  chr <- as.character(SNPInfo$CHR)
  if (anyNA(snp) || anyNA(chr) || any(!nzchar(snp)) || any(!nzchar(chr)) ||
      anyDuplicated(snp)) {
    stop("SNPInfo$SNP and SNPInfo$CHR must be complete, with unique SNP IDs.")
  }
  pos <- match(colnames(X), snp)
  if (anyNA(pos) || length(snp) != p || any(!snp %in% colnames(X))) {
    stop("SNPInfo$SNP must match colnames(X) exactly.")
  }
  chr <- chr[pos]
  chr_levels <- unique(chr)
  if (has_prs && !all(chr_levels %in% colnames(PRS))) {
    stop("Every SNPInfo$CHR value used by X must exactly match a PRS column name.")
  }

  if (is.null(Z)) stop("Z must be a non-empty numeric matrix; NULL is not supported.")
  Z <- as.matrix(Z)
  if (!is.numeric(Z) || nrow(Z) != n || !ncol(Z)) {
    stop("Z must be a non-empty numeric matrix with nrow(X) rows.")
  }
  if (any(!is.finite(X)) || any(!is.finite(Z)) ||
      (has_prs && any(!is.finite(PRS)))) {
    stop("X, Z, and PRS must contain only finite values.")
  }
  colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  suff_block_size <- validate_suff_block_size(suff_block_size)
  if (!is.numeric(r_thres) || length(r_thres) != 1L ||
      !is.finite(r_thres) || r_thres < 0 || r_thres >= 1) {
    stop("r_thres must be a finite scalar in [0, 1).")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1) {
    stop("max_iter must be a positive numeric scalar.")
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite numeric scalar.")
  }

  is_surv_y <- inherits(y, "Surv")
  if (is_surv_y) {
    status <- as.integer(y[, 2])
    y <- as.numeric(y[, 1])
  }
  if (NROW(y) != n) stop("NROW(y) must equal nrow(X).")

  family_string <- NULL
  if (is.character(family)) {
    if (length(family) != 1L) stop("family must be a scalar string or a family object.")
    family_string <- tolower(family)
  }
  if (is.null(family)) {
    if (is_surv_y || !is.null(status)) {
      family_string <- "cox"
    } else if (is.ordered(y)) {
      stop("Specify an ordinal family, for example family = 'clm_probit'.")
    } else {
      yy <- unique(stats::na.omit(as.numeric(y)))
      family <- if (length(yy) <= 2L && all(yy %in% c(0, 1))) {
        stats::binomial(link = "logit")
      } else {
        stats::gaussian()
      }
    }
  }

  path <- "glm"
  clm_link <- NULL
  if (!is.null(family_string) &&
      family_string %in% c("cox", "coxph", "survival")) {
    path <- "cox"
  } else if (!is.null(family_string) && startsWith(family_string, "clm_")) {
    clm_link <- sub("^clm_", "", family_string)
    clm_link <- ocat_validate_link(clm_link)
    path <- if (identical(clm_link, "logit")) "ocat" else "ordinal"
  } else if ((!is.null(family_string) &&
              family_string %in% c("ocat", "ordinal", "ordered",
                                   "ordered.categorical", "ordered_categorical")) ||
             ocat_is_family(family)) {
    path <- "ocat"
  } else {
    if (!is.null(family_string)) {
      if (family_string %in% c("gaussian", "normal", "linear")) {
        family <- stats::gaussian()
      } else if (family_string %in% c("binomial", "logit", "logistic")) {
        family <- stats::binomial(link = "logit")
      } else if (family_string %in%
                 c("negbin", "nb", "negative.binomial", "negative_binomial")) {
        family <- mgcv::nb(theta = NULL)
      } else {
        stop("Unsupported family.")
      }
    }
    if (!inherits(family, "family")) {
      stop("family must be NULL, a supported string, or a family object.")
    }
    if (identical(family$family, "gaussian") &&
        identical(family$link, "identity")) {
      path <- "gaussian"
    } else {
      validate_mgcv_irls_family(family)
    }
  }

  PRS_all <- if (has_prs) rowSums(PRS) else NULL
  Q_all <- if (has_prs) cbind(Z, PRS_All = PRS_all) else Z
  s <- numeric(p)
  d <- numeric(p)
  var_y_chr <- numeric(if (has_prs) length(chr_levels) else 1L)
  p_chr <- if (has_prs) integer(length(chr_levels)) else p
  XtX_blocks <- if (r_thres > 0) vector("list", length(chr_levels)) else NULL
  Xty_blocks <- if (r_thres > 0) vector("list", length(chr_levels)) else NULL
  block_jitter <- numeric(length(chr_levels))
  sigma2 <- 1

  if (identical(path, "gaussian")) {
    yy <- as.numeric(y)
    Q_fit <- cbind(Intercept = 1, Q_all)
    fit <- stats::lm.fit(Q_fit, yy)
    if (fit$rank < ncol(Q_fit)) stop("Z and PRS_All are rank deficient.")
    sigma2 <- sum(fit$residuals^2) / (n - fit$rank)
    if (!is.finite(sigma2) || sigma2 <= 0) {
      stop("The Gaussian nuisance fit produced an invalid residual variance.")
    }
    work_y <- yy
    weights <- rep(1, n)
    n_ss <- n
  } else if (identical(path, "glm")) {
    response_info <- mgcv_prepare_response(y, family)
    pred <- mgcv_predictor_data(Z = Q_all, n = n)
    dat <- cbind(response_info$data, pred)
    fit <- mgcv_fit_explicit(
      response_info$response, colnames(pred), dat, family,
      mgcv_model = mgcv_model
    )
    work <- extract_mgcv_working(fit, weight_cutoff = 0.0025,
                                 eta_clip_range = c(-50, 50))
    work_y <- work$pseudo_response
    weights <- work$weights
    n_ss <- max(0.95 * n, work$n_eff)
  } else if (identical(path, "ocat")) {
    y_info <- ocat_prepare_response(y, family = if (ocat_is_family(family)) family else NULL)
    if (!ocat_is_family(family)) family <- mgcv::ocat(R = y_info$ncat)
    pred <- mgcv_predictor_data(Z = Q_all, n = n)
    dat <- cbind(data.frame(y = as.integer(y_info$y)), pred)
    fit <- mgcv_fit_explicit("y", colnames(pred), dat, family,
                             mgcv_model = mgcv_model)
    eta_work <- pmin(pmax(as.numeric(fit$linear.predictors), -50), 50)
    weights <- fit$working.weights
    if (is.null(weights)) weights <- fit$weights
    weights <- as.numeric(weights)
    work_y <- eta_work + as.numeric(fit$residuals)
    bad <- !is.finite(work_y) | !is.finite(weights) | weights <= 0
    if (mean(bad) > 0.9) stop("Too many invalid OCAT working observations.")
    work_y[bad] <- 0
    weights[bad] <- 0
    weights <- robust_weight(weights, cutoff = 0.0025)
    n_ss <- max(0.95 * n, sum(weights)^2 / sum(weights^2))
  } else if (identical(path, "cox")) {
    if (is.null(status)) stop("status is required for Cox models unless y is a Surv object.")
    status <- as.integer(status)
    if (length(status) != n || any(!status %in% c(0L, 1L))) {
      stop("status must contain one 0/1 value per row of X.")
    }
    surv_y <- survival::Surv(as.numeric(y), status)
    dat <- as.data.frame(Q_all)
    fit <- survival::coxph(surv_y ~ ., data = dat, ties = "breslow",
                           model = FALSE, x = FALSE, y = FALSE)
    eta <- pmin(pmax(as.numeric(fit$linear.predictors), -50), 50)
    n_ss <- n
  } else {
    y_info <- ocat_prepare_response(y, family = if (ocat_is_family(family)) family else NULL)
    pred <- ocat_predictor_data(Z = Q_all, n = n)
    fit <- ocat_fit_explicit(y_info$y, pred, clm_link = clm_link)
    eta <- ocat_linear_predictor(fit, pred, n = n)
    n_ss <- n
  }

  project_block <- function(X_block, Q_loco) {
    if (path %in% c("gaussian", "glm", "ocat")) {
      ZI <- cbind(Intercept = 1, Q_loco)
      ss <- weighted_projected_suffstats(
        X = X_block, y = work_y, ZI = ZI,
        weights = weights, nuisance_precision = numeric(0),
        n_threads = n_threads,
        block_size = suff_block_size
      )
    } else if (identical(path, "cox")) {
      ss <- cox_suffstat_block(
        Xblk = X_block, eta = eta, Znui = Q_loco,
        surv_time = as.numeric(y), surv_status = status,
        nuisance_precision = numeric(0),
        n_threads = n_threads, block_size = suff_block_size
      )
      ss$yty <- n - 1
    } else {
      ss <- ocat_suffstat_block(
        Xblk = X_block, y_int = y_info$y_int,
        eta = eta, Znui = Q_loco, alpha = fit$alpha,
        clm_link = clm_link, n_threads = n_threads,
        block_size = suff_block_size
      )
    }
    ss
  }

  sparsify_block <- function(X_block, B, threshold) {
    LD <- crossprod(X_block)
    LD_diag <- diag(LD)
    if (any(!is.finite(LD_diag)) || any(LD_diag <= 0)) {
      stop("The genotype LD matrix is invalid.")
    }
    LD_cor <- LD / sqrt(tcrossprod(LD_diag))
    B <- as.matrix(B)
    B_diag <- diag(B)
    B[abs(LD_cor) < threshold] <- 0
    diag(B) <- B_diag
    B <- (B + t(B)) / 2

    ev_min <- min(eigen(B, symmetric = TRUE, only.values = TRUE)$values)
    ev_floor <- 1e-10 * max(B_diag)
    jitter <- max(0, ev_floor - ev_min)
    if (jitter > 0) diag(B) <- diag(B) + jitter
    list(
      XtX = Matrix::forceSymmetric(Matrix::Matrix(B, sparse = TRUE)),
      jitter = jitter
    )
  }

  if (has_prs) {
    for (i in seq_along(chr_levels)) {
      cc <- chr_levels[i]
      idx <- which(chr == cc)
      X_block <- X[, idx, drop = FALSE]
      p_chr[i] <- length(idx)
      PRS_loco <- PRS_all - PRS[, cc]
      Q_loco <- cbind(Z, PRS_LOCO = PRS_loco)
      ss <- project_block(X_block, Q_loco)

      s[idx] <- as.numeric(ss$Xty)
      d[idx] <- diag(ss$XtX)
      var_y_chr[i] <- ss$yty / (n_ss - 1)

      if (r_thres > 0) {
        sparse_block <- sparsify_block(X_block, ss$XtX, r_thres)
        XtX_blocks[[i]] <- sparse_block$XtX
        Xty_blocks[[i]] <- as.numeric(ss$Xty)
        block_jitter[i] <- sparse_block$jitter
      }
    }
  } else {
    ss <- project_block(X, Z)
    s <- as.numeric(ss$Xty)
    d <- diag(ss$XtX)
    var_y_chr[1] <- ss$yty / (n_ss - 1)

    if (r_thres > 0) {
      for (i in seq_along(chr_levels)) {
        idx <- which(chr == chr_levels[i])
        sparse_block <- sparsify_block(
          X[, idx, drop = FALSE], ss$XtX[idx, idx, drop = FALSE], r_thres
        )
        XtX_blocks[[i]] <- sparse_block$XtX
        Xty_blocks[[i]] <- as.numeric(ss$Xty[idx])
        block_jitter[i] <- sparse_block$jitter
      }
    }
  }

  if (any(!is.finite(s)) || any(!is.finite(d)) || any(d <= 0) ||
      any(!is.finite(var_y_chr)) || any(var_y_chr <= 0)) {
    stop("The projected sufficient statistics are invalid.")
  }

  Vhat <- 2
  converged <- FALSE
  xtx_density <- NA_real_
  if (r_thres == 0) {
    for (iter in seq_len(max_iter)) {
      den <- d + sigma2 / Vhat
      m <- s / den
      u <- sigma2 / den
      V1 <- mean(m^2 + u)
      if (abs(V1 - Vhat) <= tol * max(Vhat, V1, 1e-12)) {
        converged <- TRUE
        Vhat <- V1
        break
      }
      Vhat <- V1
    }
  } else {
    XtX_sparse <- Matrix::bdiag(XtX_blocks)
    Xty_sparse <- unlist(Xty_blocks, use.names = FALSE)
    xtx_density <- Matrix::nnzero(XtX_sparse) / length(XtX_sparse)
    for (iter in seq_len(max_iter)) {
      P <- Matrix::forceSymmetric(
        XtX_sparse + Matrix::Diagonal(p, x = sigma2 / Vhat)
      )
      F <- Matrix::Cholesky(P, LDL = FALSE, perm = TRUE)
      m <- as.numeric(Matrix::solve(F, Xty_sparse))
      Pinv <- Matrix::solve(F, Matrix::Diagonal(p))
      V1 <- (sum(m^2) + sigma2 * sum(Matrix::diag(Pinv))) / p
      if (abs(V1 - Vhat) <= tol * max(Vhat, V1, 1e-12)) {
        converged <- TRUE
        Vhat <- V1
        break
      }
      Vhat <- V1
    }
  }

  var_y_ss <- stats::weighted.mean(var_y_chr, p_chr)
  diagnostics <- data.frame(
    n = n,
    p = p,
    n_chr = length(chr_levels),
    n_ss = n_ss,
    residual_variance = sigma2,
    var_y_ss = var_y_ss,
    r_thres = r_thres,
    xtx_density = xtx_density,
    max_block_jitter = max(block_jitter),
    iterations = iter,
    converged = converged
  )

  list(
    prior_variance = as.numeric(Vhat),
    scaled_prior_variance = as.numeric(Vhat / var_y_ss),
    diagnostics = diagnostics
  )
}
