Run_GG_CLM <- function(X, y, clm_link,
                         Lmain, Lint, max.iter, min.iter, max.eps,
                         susie_para_main, susie_para_int,
                         verbose = TRUE, n_threads = 1,
                         x_noncs_var = 0.1,
                         w_noncs_var = 0.1,
                         noncs_max_abs_cor = 0.9,
                         include_x_squared = FALSE,
                         suff_block_size = 10000L,
                         returnModel = FALSE) {
  run_start <- proc.time()[["elapsed"]]
  n <- NROW(y)
  p <- ncol(X)
  ridge <- 1e-8
  eta_clip_range <- c(-50, 50)
  min_etaX_var <- 1e-7
  min_etaW_var <- 1e-7
  clm_link <- ocat_validate_link(clm_link)
  if (identical(clm_link, "logit")) {
    stop("clm_logit must use the mgcv::ocat runner.")
  }
  y_info <- ocat_prepare_response(y)
  y <- y_info$y
  y_int <- y_info$y_int
  nZ <- 0L
  suff_block_size <- validate_suff_block_size(suff_block_size)

  pred_final <- ocat_predictor_data(n = n)
  fit_final <- ocat_fit_explicit(
    y, pred_final, clm_link = clm_link, family = NULL
  )
  theta <- ocat_nuisance_coef(fit_final, n_keep = nZ)

  g <- c()
  beta <- rep(0, p)
  XCS <- NULL
  XCS_refit <- NULL
  WCS <- NULL
  WCS_refit <- NULL
  W <- NULL
  fitX <- NULL
  fitW <- NULL
  ssX <- NULL

  fitX_no_cs_streak <- 0L
  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    theta_prev <- theta

    eta <- ocat_linear_predictor(fit_final, pred_final, n = n)
    eta <- pmin(pmax(eta, eta_clip_range[1]), eta_clip_range[2])
    Znui_main <- cbind(NULL, WCS_refit)
    ssX <- ocat_suffstat_block(
      X, y_int = y_int, eta = eta, Znui = Znui_main,
      alpha = ocat_thresholds(fit_final), clm_link = clm_link,
      n_threads = n_threads, ridge = ridge,
      block_size = suff_block_size
    )
    fitX <- .fit_susie_stage(
      structural = list(
        XtX = ssX$XtX, Xty = ssX$Xty, yty = ssX$yty,
        n = n, L = Lmain
      ),
      susie_para = susie_para_main, stage = "main",
      iter = iter, min.iter = min.iter
    )
    beta <- coef.susie(fitX)[-1]

    CSdt <- summary(fitX)$vars
    x_component <- build_component_design_from_fit(
      X, fitX, "Main_CS"
    )
    cs_indices <- x_component$cs_indices
    fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L
    main_no_cs <- is.null(x_component$design)
    if (is.null(x_component$design)) {
        noncs_main <- build_full_noncs_refit_term(X, fitX)
        XCS <- NULL
        if (is.null(noncs_main)) {
          XCS_refit <- NULL
        } else {
          XCS_refit <- matrix(noncs_main, ncol = 1)
          colnames(XCS_refit) <- "Main_noncs_res"
        }
        W <- NULL
        WCS <- NULL
        WCS_refit <- NULL
        fitW <- NULL
    } else {
      XCS <- x_component$design
      XCS_refit <- XCS
      {
        main_noncs <- build_noncs_refit_term(
          X = X, fit = fitX, CSdt = CSdt, cs_indices = cs_indices,
          XCS = XCS, noncs_var = x_noncs_var,
          min_eta_var = min_etaX_var,
          noncs_max_abs_cor = noncs_max_abs_cor,
          cor_design = NULL
        )
        XCS_refit <- append_noncs_refit_term(
          XCS_refit, main_noncs, "Main_noncs_res",
          corr_threshold = noncs_max_abs_cor
        )
      }
    }

      pred_final <- ocat_predictor_data(
        Z = NULL,
        Xextra = cbind(XCS_refit, WCS_refit), n = n
      )
      penalty_names <- refit_penalty_terms(colnames(pred_final))
      penalty_V <- ocat_penalty_variance(
        clm_link, fitX, fitW, penalty_names
      )
      fit_final <- ocat_fit_explicit(
        y, pred_final, alpha_start = ocat_thresholds(fit_final),
        clm_link = clm_link, family = NULL, penalty_V = penalty_V
      )
      eta <- ocat_linear_predictor(fit_final, pred_final, n = n)
      eta <- pmin(pmax(eta, eta_clip_range[1]), eta_clip_range[2])

      XCS_W <- XCS_refit
      W <- get_pairwise_interactions(
        XCS_W, Z = NULL,
        include_x_squared = if (main_no_cs) FALSE else include_x_squared
      )
      if (!interaction_design_available(W, iter, min.iter, allow_empty = main_no_cs)) {
        W <- NULL
        fitW <- NULL
        WCS <- NULL
        WCS_refit <- NULL
      } else {
      Znui_int <- cbind(NULL, XCS_refit)
      ssW <- ocat_suffstat_block(
        W, y_int = y_int, eta = eta, Znui = Znui_int,
        alpha = ocat_thresholds(fit_final), clm_link = clm_link,
        n_threads = n_threads, ridge = ridge,
        block_size = suff_block_size
      )
      fitW <- .fit_susie_stage(
        structural = list(
          XtX = ssW$XtX, Xty = ssW$Xty, yty = ssW$yty,
          n = n, L = Lint
        ),
        susie_para = susie_para_int, stage = "int",
        iter = iter, min.iter = min.iter
      )

      CSdt_w <- summary(fitW)$vars
      cs_w <- sort(unique(CSdt_w$cs[CSdt_w$cs > 0]))
      w_component <- build_component_design_from_fit(
        W, fitW, "Int_CS"
      )
      WCS <- w_component$design
      WCS_refit <- WCS

      {
      w_noncs <- build_w_noncs_refit_term(
        W = W, fitW = fitW, WCS = WCS,
        etaX = matrixVectorMultiply(X, beta),
        XCS = XCS, Z = NULL,
        w_noncs_var = w_noncs_var, min_etaW_var = min_etaW_var,
        noncs_max_abs_cor = noncs_max_abs_cor
      )
      if (!is.null(w_noncs)) {
        if (is.null(WCS_refit)) {
          WCS_refit <- matrix(w_noncs, ncol = 1)
          colnames(WCS_refit) <- "Int_noncs_res"
        } else {
          WCS_refit <- append_noncs_refit_term(
            WCS_refit, w_noncs, "Int_noncs_res",
            corr_threshold = noncs_max_abs_cor
          )
        }
      }
      }
      }

    pred_final <- ocat_predictor_data(
      Z = NULL,
      Xextra = cbind(XCS_refit, WCS_refit), n = n
    )
    penalty_names <- refit_penalty_terms(colnames(pred_final))
    penalty_V <- ocat_penalty_variance(
      clm_link, fitX, fitW, penalty_names
    )
    fit_final <- ocat_fit_explicit(
      y, pred_final, alpha_start = ocat_thresholds(fit_final),
      clm_link = clm_link, family = NULL, penalty_V = penalty_V
    )
    theta <- ocat_nuisance_coef(fit_final, n_keep = nZ)

    errX <- sqrt(mean((beta - beta_prev)^2))
    errT <- sqrt(mean((theta - theta_prev)^2))
    err <- max(errX, errT)
    g[iter] <- err
    if (verbose) {
      cat(sprintf(
        "Iteration %d: err = %.3e, n = %d, min_pr = %.3g, med_h = %.3g\n",
        iter, err, ssX$n_eff, ssX$min_pr, ssX$med_h
      ))
    }
    if (fitX_no_cs_streak >= 3L) {
      if (verbose) cat("No main credible set detected in 3 consecutive iterations; stopping.\n")
      break
    }
    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  pred_final <- ocat_predictor_data(
    Z = NULL,
    Xextra = cbind(XCS_refit, WCS_refit), n = n
  )
  penalty_names <- refit_penalty_terms(colnames(pred_final))
  penalty_V <- ocat_penalty_variance(
    clm_link, fitX, fitW, penalty_names
  )
  fit_final <- ocat_fit_explicit(y, pred_final, alpha_start = ocat_thresholds(fit_final),
                                 clm_link = clm_link, family = NULL,
                                 penalty_V = penalty_V)
  fit_final$n_eff <- ssX$n_eff
  G <- ocat_coef_table(fit_final)
  MainIndex <- Identifying_MainEffect(fitX, colnames(X))
  MainIndex <- safe_add_p(MainIndex, G)
  IntIndex <- Identifying_IntEffect(fitW, colnames(W))
  IntIndex <- filter_noncs_interactions(IntIndex)
  IntIndex <- safe_add_p(IntIndex, G)

  if (verbose) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration", ylab = "Max Parameter Change",
         main = paste0("Convergence Trace (Ordered Categorical, ", clm_link, ")"))
    for (i in seq_along(g)) {
      text(x = i, y = g[i], labels = formatC(g[i], format = "e", digits = 1),
           pos = 3, cex = 0.7, col = "red")
    }
  }

  AA <- list(
    diagnostics = make_diagnostics(iter, g, run_start),
    fitX = fitX,
    fitW = fitW,
    fitJoint = fit_final,
    main_discoveries = MainIndex,
    interaction_discoveries = IntIndex
  )
  if (returnModel) AA$FinalModel <- pred_final
  AA$report <- extract_direction_table(AA, G)
  return(AA)
}
