#' Cox sufficient statistics for one design block, projected onto the nuisance complement
#'
#' Internal helper shared by the Cox runners. Given a design block Xblk to fine-map,
#' the current linear predictor eta, and nuisance covariates Znui (an intercept
#' is always appended for the projection), it builds the Breslow observed information
#' and score statistic for the Xblk block after projecting out the nuisance columns,
#' then returns the projected XtX and Xty with n_eff equal to the number of events.
#'
#' @keywords internal
#' @noRd
cox_suffstat_block <- function(Xblk, eta, Znui, surv_time, surv_status,
                               nuisance_precision,
                               n_threads = 1, ridge = 1e-6,
                               block_size = 10000L) {
  if (missing(nuisance_precision)) {
    stop("nuisance_precision must be supplied explicitly for every projection.")
  }
  Xblk <- as.matrix(Xblk)
  n <- nrow(Xblk)
  p <- ncol(Xblk)

  if (is.null(Znui)) {
    ZI <- matrix(1, n, 1)
    colnames(ZI) <- "Intercept"
  } else {
    Znui <- as.matrix(Znui)
    ZI <- cbind(Intercept = 1, Znui)
  }
  q <- ncol(ZI)
  projection_precision <- align_projection_precision(ZI, nuisance_precision)

  N <- cbind(eta, ZI)
  k <- ncol(N)
  status <- as.integer(surv_status)

  rsX <- cox_riskset(X = Xblk, eta = eta, time = surv_time,
                     status = status, n_threads = n_threads)
  rsN <- cox_riskset(X = N, eta = eta, time = surv_time,
                     status = status, n_threads = 1L)
  a <- as.numeric(rsX$a)
  M <- as.numeric(rsX$M)
  dev <- as.numeric(rsX$dev)
  n_eff <- rsX$d

  # [X N]' diag(a) [X N] - B' diag(dev) B, split into X and N blocks.
  AX <- weighted_crossprod(Xblk, a, cbind(N * a, M),
                           n_threads = n_threads, block_size = block_size)
  BX <- weighted_crossprod(rsX$B, dev, rsN$B * dev,
                           n_threads = n_threads, block_size = block_size)
  XX <- AX$XtWX - BX$XtWX
  XN <- AX$XtM[, seq_len(k), drop = FALSE] - BX$XtM
  NN <- crossprod(N, N * a) - crossprod(rsN$B, rsN$B * dev)
  NN <- (NN + t(NN)) / 2

  XtX <- (XX + t(XX)) / 2
  dimnames(XtX) <- list(colnames(Xblk), colnames(Xblk))
  XtE <- XN[, 1L, drop = FALSE]
  XtZ <- XN[, 1L + seq_len(q), drop = FALSE]
  ZtZ <- NN[1L + seq_len(q), 1L + seq_len(q), drop = FALSE]
  diag(ZtZ) <- diag(ZtZ) + projection_precision
  ZtX <- t(XtZ)
  ZtE <- NN[1L + seq_len(q), 1L, drop = FALSE]
  XtM <- as.numeric(AX$XtM[, k + 1L])
  ZtM <- as.numeric(crossprod(ZI, M))

  Zinv_ZtX <- solve_with_ridge(ZtZ, ZtX, ridge = ridge)
  Zinv_ZtE_score <- solve_with_ridge(
    ZtZ, matrix(ZtE + ZtM, ncol = 1L), ridge = ridge
  )

  XtX_proj <- XtX - matrixMultiply(XtZ, Zinv_ZtX)
  XtE_proj <- as.vector(
    XtE + XtM - matrixVectorMultiply(XtZ, Zinv_ZtE_score)
  )

  Xty <- XtE_proj
  XtX <- (XtX_proj + t(XtX_proj)) / 2
  XtX_pre_ridge <- XtX
  diag(XtX) <- diag(XtX) + ridge

  dXtX <- diag(XtX)
  zhat <- Xty / sqrt(dXtX)
  R <- cov2cor(XtX)
  R <- (R + t(R)) / 2

  list(
    z = zhat, R = R, n_eff = n_eff, XtX = XtX,
    XtX_pre_ridge = XtX_pre_ridge, Xty = Xty
  )
}
