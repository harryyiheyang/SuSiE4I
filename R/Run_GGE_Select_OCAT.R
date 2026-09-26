Run_GGE_Select_OCAT <- function(X, Z, y, family = NULL, mgcv_model = NULL, Lmain, Lenv, Lint,
    max.iter, min.iter, max.eps, susie_para_main, susie_para_int, susie_para_env, verbose = TRUE, n_threads = 1, x_noncs_var = 0.1,
    w_noncs_var = 0.1, noncs_max_abs_cor = 0.9, include_x_squared = FALSE, suff_block_size = 10000L, L.init = 1, returnModel = FALSE) {
    run_start <- proc.time()[["elapsed"]]
    n <- NROW(y)
    p <- ncol(X)
    suff_block_size <- validate_suff_block_size(suff_block_size)
    if (is.null(dim(Z)))
        Z <- matrix(Z, ncol = 1)
    if (is.null(colnames(Z)))
        colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    eta_clip_range <- c(-50, 50)
    weight_cutoff <- 0.0025
    min_etaX_var <- 1e-07
    min_etaW_var <- 1e-07
    y_info <- ocat_prepare_response(y, family)
    y <- y_info$y
    if (is.null(family)) family <- mgcv::ocat(R = y_info$ncat)
    if (!ocat_is_family(family)) {
        stop("Run_GGE_Select_OCAT requires mgcv::ocat().")
    }
    response_info <- list(
        data = data.frame(y = as.integer(y)), response = "y", n = length(y)
    )
    if (response_info$n != nrow(X))
        stop("Length(y) must equal nrow(X).")
    env_init <- initial_continuous_env(Z, y, Lenv, n_threads = n_threads, block_size = suff_block_size)
    ZCS <- env_init$ZCS
    pred0 <- mgcv_predictor_data(Z = ZCS, n = n)
    Data0 <- cbind(response_info$data, pred0)
    fit_final <- {
        mgcv_fit_explicit(response_info$response, colnames(pred0), Data0, family, mgcv_model = mgcv_model)
    }
    g <- c()
    beta <- rep(0, p)
    alpha <- rep(0, ncol(Z))
    gamma <- 0
    beta_prev <- beta
    alpha_prev <- alpha
    eta0 <- coef(fit_final)[1L]
    if (!is.finite(eta0))
        eta0 <- 0
    etaX <- 0
    etaW <- 0
    etaZ <- 0
    ZCS_refit <- NULL
    XCS_refit <- NULL
    WCS_refit <- NULL
    fitZ <- NULL
    fitX <- NULL
    fitW <- NULL
    fitX_no_cs_streak <- 0L
    for (iter in 1:max.iter) {
        beta_prev <- beta
        alpha_prev <- alpha
        eta_work <- pmin(pmax(as.numeric(fit_final$linear.predictors),
            eta_clip_range[1]), eta_clip_range[2])
        work_residual <- as.numeric(fit_final$residuals)
        work_weight <- fit_final$working.weights
        if (is.null(work_weight)) work_weight <- fit_final$weights
        work_weight <- as.numeric(work_weight)
        if (length(work_residual) != n || length(work_weight) != n) {
            stop("mgcv::ocat did not return one working residual and weight per observation.")
        }
        pseudo_response <- eta_work + work_residual
        bad <- !is.finite(pseudo_response) | !is.finite(work_weight) |
            work_weight <= 0
        if (mean(bad) > 0.9) stop("Too many invalid OCAT working observations.")
        if (any(bad)) {
            pseudo_response[bad] <- 0
            work_weight[bad] <- 0
        }
        work_weight <- robust_weight(work_weight, cutoff = weight_cutoff)
        work <- list(
            pseudo_response = pseudo_response, weights = work_weight,
            phi0 = 1, n_eff = sum(work_weight)^2 / sum(work_weight^2)
        )
        pseudo_response <- work$pseudo_response
        W_diag <- work$weights
        ZI_env <- cbind(Intercept = rep(1, n), XCS_refit, WCS_refit)
        ssZ <- weighted_projected_suffstats(X = Z, y = pseudo_response, ZI = ZI_env, weights = W_diag,
            nuisance_precision = projection_penalty_precision(ZI_env, fitX, fitW, fitZ), n_threads = n_threads,
            block_size = suff_block_size)
        ZtZ_w <- ssZ$XtX
        Zty_w <- ssZ$Xty
        yty4Z <- ssZ$yty
        fitZ <- .fit_susie_stage(structural = list(XtX = ZtZ_w, Xty = Zty_w, yty = yty4Z, n = max(0.95 * n, work$n_eff), L = Lenv),
            susie_para = susie_para_env, stage = "env", iter = iter, min.iter = min.iter, gaussian = FALSE, residual_variance = work$phi0)
        alpha <- coef(fitZ)[-1]
        etaZ <- matrixVectorMultiply(Z, alpha)
        env_terms <- build_select_env_terms(Z, fitZ, noncs_env_var = w_noncs_var, min_eta_env_var = min_etaW_var, noncs_max_abs_cor = noncs_max_abs_cor)
        if (is.null(env_terms$ZCS))
            stop_select_env_no_cs()
        ZCS <- env_terms$ZCS
        ZCS_refit <- env_terms$ZCS_refit
        ZI_main <- cbind(Intercept = 1, ZCS_refit, WCS_refit)
        ssX <- weighted_projected_suffstats(X = X, y = pseudo_response, ZI = ZI_main, weights = W_diag,
            nuisance_precision = projection_penalty_precision(ZI_main, fitX, fitW, fitZ), n_threads = n_threads,
            block_size = suff_block_size)
        XtX <- {
            ssX$XtX
        }
        Xty <- ssX$Xty
        yty4X <- ssX$yty
        fitX <- .fit_susie_stage(structural = list(XtX = XtX, Xty = Xty, yty = yty4X, n = max(0.95 * n, work$n_eff), L = Lmain),
            susie_para = susie_para_main, stage = "main", iter = iter, min.iter = min.iter, gaussian = FALSE, residual_variance = work$phi0)
        beta <- coef.susie(fitX)[-1]
        etaX <- matrixVectorMultiply(X, beta)
        CSdt <- summary(fitX)$vars
        x_component <- build_component_design_from_fit(X, fitX, "Main_CS")
        cs_indices <- x_component$cs_indices
        fitX_no_cs_streak <- if (length(cs_indices)) 0L else fitX_no_cs_streak + 1L
        main_no_cs <- is.null(x_component$design)
        if (main_no_cs) {
            noncs_main <- build_full_noncs_refit_term(X, fitX)
            XCS <- NULL
            if (is.null(noncs_main)) {
                XCS_refit <- NULL
            }
            else {
                XCS_refit <- matrix(noncs_main, ncol = 1)
                colnames(XCS_refit) <- "Main_noncs_res"
            }
        }
        else {
            XCS <- x_component$design
            XCS_refit <- XCS
            {
                noncs_main <- build_noncs_refit_term(X = X, fit = fitX, CSdt = CSdt, cs_indices = cs_indices, XCS = XCS,
                  noncs_var = x_noncs_var, min_eta_var = min_etaX_var, noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Z)
                if (!is.null(noncs_main)) {
                  XCS_refit <- cbind(XCS_refit, Main_noncs_res = noncs_main)
                }
            }
        }
        XCS_W <- XCS_refit
        W <- get_pairwise_interactions(XCS_W, ZCS, include_x_squared = if (main_no_cs)
            FALSE
        else include_x_squared)
        WCS <- NULL
        WCS_refit <- NULL
        if (!interaction_design_available(W, iter, min.iter, allow_empty = main_no_cs)) {
            W <- NULL
            fitW <- NULL
            WCS <- NULL
            WCS_refit <- NULL
        }
        else {
            ZI_int <- cbind(Intercept = 1, ZCS_refit, XCS_refit)
            ssW <- weighted_projected_suffstats(W, pseudo_response, ZI_int, W_diag,
              nuisance_precision = projection_penalty_precision(ZI_int, fitX, fitW, fitZ),
              n_threads = n_threads, block_size = suff_block_size)
            WtW <- ssW$XtX
            Wty <- ssW$Xty
            yty4W <- ssW$yty
            fitW <- .fit_susie_stage(structural = list(XtX = WtW, Xty = Wty, yty = yty4W, n = max(0.95 * n, work$n_eff), L = Lint),
                susie_para = susie_para_int, stage = "int", iter = iter, min.iter = min.iter, gaussian = FALSE, residual_variance = work$phi0)
            CSdt_w <- summary(fitW)$vars
            cs_indices_w <- unique(CSdt_w$cs[CSdt_w$cs > 0])
            w_component <- build_component_design_from_fit(W, fitW, "Int_CS")
            WCS <- w_component$design
            WCS_refit <- WCS
            {
                w_noncs <- build_w_noncs_refit_term(W = W, fitW = fitW, WCS = WCS, etaX = etaX, XCS = XCS, Z = Z, w_noncs_var = w_noncs_var,
                  min_etaW_var = min_etaW_var, noncs_max_abs_cor = noncs_max_abs_cor)
                if (!is.null(w_noncs)) {
                  if (is.null(WCS_refit)) {
                    WCS_refit <- matrix(w_noncs, ncol = 1)
                    colnames(WCS_refit) <- "Int_noncs_res"
                  }
                  else {
                    WCS_refit <- append_noncs_refit_term(WCS_refit, w_noncs, "Int_noncs_res", corr_threshold = noncs_max_abs_cor)
                  }
                }
            }
        }
        gamma <- if (is.null(fitW))
            0
        else coef.susie(fitW)[-1]
        pred <- mgcv_predictor_data(Xextra = cbind(ZCS_refit, XCS_refit, WCS_refit), n = n)
        Data <- cbind(response_info$data, pred)
        penalty_names <- refit_penalty_terms(colnames(pred))
        penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names, fitZ = fitZ)
        fit_final <- {
            mgcv_fit_fixed_ridge(response_info$response, colnames(pred), Data, family, penalty_V, dispersion = work$phi0,
                mgcv_model = mgcv_model)
        }
        coefs <- coef(fit_final)
        idx <- 1
        eta0 <- coefs[1]
        if (!is.finite(eta0))
            eta0 <- 0
        ZCSbeta <- coefs[(idx + 1):(idx + ncol(ZCS_refit))]
        idx <- idx + ncol(ZCS_refit)
        etaZ <- matrixVectorMultiply(ZCS_refit, ZCSbeta)
        if (!is.null(XCS_refit)) {
            XCSbeta <- coefs[(idx + 1):(idx + ncol(XCS_refit))]
            idx <- idx + ncol(XCS_refit)
            etaX <- matrixVectorMultiply(XCS_refit, XCSbeta)
        }
        else {
            etaX <- rep(0, n)
        }
        if (!is.null(WCS_refit)) {
            WCSbeta <- coefs[(idx + 1):(idx + ncol(WCS_refit))]
            etaW <- matrixVectorMultiply(WCS_refit, WCSbeta)
        }
        else {
            etaW <- 0
        }
        errX <- sqrt(mean((beta - beta_prev)^2))
        errZ <- sqrt(mean((alpha - alpha_prev)^2))
        err <- max(errX, errZ)
        g[iter] <- err
        if (verbose)
            cat(sprintf("Iteration %d: err = %.3e\n", iter, err))
        if (fitX_no_cs_streak >= 3L) {
          if (verbose) cat("No main credible set detected in 3 consecutive iterations; stopping.\n")
          break
        }
        if (err < max.eps && iter > min.iter) {
            if (verbose)
                cat("Converged!\n")
            break
        }
    }
    pred_list <- list()
    if (!is.null(ZCS_refit))
        pred_list$ZCS <- ZCS_refit
    pred_list$XCS <- XCS_refit
    if (!is.null(WCS_refit))
        pred_list$WCS <- WCS_refit
    pred <- as.data.frame(do.call(cbind, pred_list))
    Dat <- cbind(response_info$data, pred)
    refit_dispersion <- mgcv_refit_dispersion(fit_final)
    penalty_names <- refit_penalty_terms(colnames(pred))
    penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names, fitZ = fitZ)
    {
        fit_final <- mgcv_fit_fixed_ridge(response_info$response, colnames(pred), Dat, family, penalty_V, dispersion = refit_dispersion,
            mgcv_model = mgcv_model)
    }
    fit_final$n_eff <- work$n_eff
    MainIndex <- Identifying_MainEffect(fitX, colnames(X))
    IntIndex <- Identifying_IntEffect(fitW, colnames(W))
    IntIndex <- filter_noncs_interactions(IntIndex)
    EnvIndex <- Identifying_EnvEffect(fitZ, colnames(Z))
    G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
    MainIndex <- safe_add_p(MainIndex, G)
    IntIndex <- safe_add_p(IntIndex, G)
    EnvIndex <- safe_add_p(EnvIndex, G)
    if (verbose) {
        plot(g, type = "o", col = "black", pch = 16, xlab = "Iteration", ylab = "Max Parameter Change", main = "Convergence Trace (Max |Delta| in alpha and beta)")
        for (i in seq_along(g)) {
            text(x = i, y = g[i], labels = formatC(g[i], format = "e", digits = 1), pos = 3, cex = 0.7, col = "red")
        }
    }
    AA <- list(diagnostics = make_diagnostics(iter, g, run_start), fitX = fitX, fitW = fitW, fitZ = fitZ, fitJoint = fit_final,
        main_discoveries = MainIndex, interaction_discoveries = IntIndex, environment_discoveries = EnvIndex)
    if (returnModel)
        AA$FinalModel <- Dat
    AA$report <- extract_direction_table(AA, G)
    return(AA)
}
