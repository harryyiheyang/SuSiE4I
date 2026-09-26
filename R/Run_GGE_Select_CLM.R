Run_GGE_Select_CLM <- function(X, Z, y, clm_link,
                                Lmain, Lenv, Lint,
                                max.iter, min.iter, max.eps,
                                susie_para_main, susie_para_int, susie_para_env,
                                verbose = TRUE, n_threads = 1,
                                x_noncs_var = 0.1,
                                w_noncs_var = 0.1,
                                noncs_max_abs_cor = 0.9,
                                include_x_squared = FALSE,
                                suff_block_size = 10000L,
                                returnModel = FALSE) {
  run_start <- proc.time()[["elapsed"]]
  ridge <- 1e-8
  eta_clip_range <- c(-50, 50)
  min_etaX_var <- 1e-7
  min_etaW_var <- 1e-7
  n <- NROW(y)
  p <- ncol(X)
  q <- ncol(Z)
  clm_link <- ocat_validate_link(clm_link)
  if (identical(clm_link, "logit")) {
    stop("clm_logit must use the mgcv::ocat runner.")
  }
  y_info <- ocat_prepare_response(y)
  y <- y_info$y
  y_int <- y_info$y_int
  if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(q))
  suff_block_size <- validate_suff_block_size(suff_block_size)

  env_init <- initial_continuous_env(
    Z, y_int, Lenv, n_threads = n_threads, block_size = suff_block_size
  )
  ZCS <- env_init$ZCS
  ZCS_refit <- env_init$ZCS_refit
  pred_initial <- ocat_predictor_data(Xextra = ZCS_refit, n = n)
  fit_final <- ocat_fit_explicit(
    y, pred_initial, clm_link = clm_link, family = NULL
  )

  g <- c()
  beta <- rep(0, p)
  alpha_env <- rep(0, q)
  etaZ <- rep(0, n)
  theta <- clean_coef(ocat_thresholds(fit_final))
  XCS <- NULL
  XCS_refit <- NULL
  WCS <- NULL
  WCS_refit <- NULL
  W <- NULL
  fitZ <- NULL
  fitX <- NULL
  fitW <- NULL
  ssX <- NULL
  pred_offset <- pred_initial
  eta_offset <- NULL

  fitX_no_cs_streak <- 0L
  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_env_prev <- alpha_env
    theta_prev <- theta

    eta <- ocat_linear_predictor(fit_final, pred_offset, n = n, offset = eta_offset)
    eta <- pmin(pmax(eta, eta_clip_range[1]), eta_clip_range[2])

    ssZ <- ocat_suffstat_block(
      Z, y_int = y_int, eta = eta,
      Znui = cbind(XCS_refit, WCS_refit),
      alpha = ocat_thresholds(fit_final), clm_link = clm_link,
      n_threads = n_threads, ridge = ridge,
      block_size = suff_block_size
    )
    fitZ <- .fit_susie_stage(
      structural = list(
        XtX = ssZ$XtX, Xty = ssZ$Xty, yty = ssZ$yty,
        n = n, L = Lenv
      ),
      susie_para = susie_para_env, stage = "env",
      iter = iter, min.iter = min.iter
    )
    alpha_env <- coef.susie(fitZ)[-1]
    etaZ <- as.numeric(CppMatrix::matrixVectorMultiply(Z, alpha_env))
    env_terms <- build_select_env_terms(
      Z, fitZ, noncs_env_var = w_noncs_var,
      min_eta_env_var = min_etaW_var,
      noncs_max_abs_cor = noncs_max_abs_cor
    )
    if (is.null(env_terms$ZCS)) stop_select_env_no_cs()
    ZCS <- env_terms$ZCS
    ZCS_refit <- env_terms$ZCS_refit

    ssX <- ocat_suffstat_block(
      X, y_int = y_int, eta = eta, Znui = cbind(ZCS_refit, WCS_refit),
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
        noncs_main <- build_noncs_refit_term(
          X = X, fit = fitX, CSdt = CSdt, cs_indices = cs_indices,
          XCS = XCS, noncs_var = x_noncs_var,
          min_eta_var = min_etaX_var,
          noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Z
        )
        XCS_refit <- append_noncs_refit_term(
          XCS_refit, noncs_main, "Main_noncs_res",
          corr_threshold = noncs_max_abs_cor
        )
      }
    }

    XCS_W <- XCS_refit
    if (!is.null(XCS_refit)) {
      W <- get_pairwise_interactions(XCS_W, ZCS,
                                     include_x_squared = if (main_no_cs) FALSE else include_x_squared)
      if (!interaction_design_available(W, iter, min.iter, allow_empty = main_no_cs)) {
        W <- NULL
        fitW <- NULL
        WCS <- NULL
        WCS_refit <- NULL
      } else {
      ssW <- ocat_suffstat_block(
        W, y_int = y_int, eta = eta, Znui = cbind(ZCS_refit, XCS_refit),
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
        XCS = XCS, Z = Z,
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
    } else {
      W <- NULL
      WCS <- NULL
      WCS_refit <- NULL
      fitW <- NULL
    }

    pred_offset <- ocat_predictor_data(
      Xextra = cbind(ZCS_refit, XCS_refit, WCS_refit), n = n
    )
    eta_offset <- NULL
    penalty_names <- refit_penalty_terms(colnames(pred_offset))
    penalty_V <- ocat_penalty_variance(
      clm_link, fitX, fitW, penalty_names, fitZ = fitZ
    )
    fit_final <- ocat_fit_explicit(
      y, pred_offset, alpha_start = ocat_thresholds(fit_final), offset = eta_offset,
      clm_link = clm_link, family = NULL, penalty_V = penalty_V
    )
    theta <- clean_coef(ocat_thresholds(fit_final))

    errX <- sqrt(mean((beta - beta_prev)^2))
    errZ <- sqrt(mean((alpha_env - alpha_env_prev)^2))
    errT <- sqrt(mean((theta - theta_prev)^2))
    err <- max(errX, errZ, errT)
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
    Xextra = cbind(ZCS_refit, XCS_refit, WCS_refit), n = n
  )
  penalty_names <- refit_penalty_terms(colnames(pred_final))
  penalty_V <- ocat_penalty_variance(
    clm_link, fitX, fitW, penalty_names, fitZ = fitZ
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
  EnvIndex <- Identifying_EnvEffect(fitZ, colnames(Z))
  EnvIndex <- safe_add_p(EnvIndex, G)

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
    fitZ = fitZ,
    fitJoint = fit_final,
    main_discoveries = MainIndex,
    interaction_discoveries = IntIndex,
    environment_discoveries = EnvIndex
  )
  if (returnModel) AA$FinalModel <- pred_final
  AA$report <- extract_direction_table(AA, G)
  return(AA)
}
