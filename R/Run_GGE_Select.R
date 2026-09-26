Run_GGE_Select <- function(X, Z, y, mgcv_model = NULL, crossprodX = NULL, Lmain, Lenv, Lint, max.iter,
    min.iter, max.eps, susie_para_main, susie_para_int, susie_para_env, verbose = TRUE, n_threads = 1, x_noncs_var = 0.1,
    w_noncs_var = 0.1, noncs_max_abs_cor = 0.9, include_x_squared = FALSE, suff_block_size = 10000L, L.init = 1, returnModel = FALSE) {
    family <- gaussian()
    run_start <- proc.time()[["elapsed"]]
    n <- NROW(y)
    p <- ncol(X)
    if (!is.null(crossprodX) && !identical(dim(as.matrix(crossprodX)), c(p, p))) {
        stop("crossprodX must be a square ncol(X)-by-ncol(X) matrix.")
    }
    suff_block_size <- validate_suff_block_size(suff_block_size)
    if (is.null(dim(Z)))
        Z <- matrix(Z, ncol = 1)
    if (is.null(colnames(Z)))
        colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    min_etaX_var <- 1e-07
    min_etaW_var <- 1e-07
    response_info <- mgcv_prepare_response(y, family)
    if (response_info$n != nrow(X))
        stop("Length(y) must equal nrow(X).")
    XtX <- if (is.null(crossprodX)) {
        blockwise_crossprod(X, n_threads = n_threads,
            block_size = suff_block_size)
    }
    else {
        as.matrix(crossprodX)
    }
    ZtZ <- blockwise_crossprod(Z, n_threads = n_threads,
        block_size = suff_block_size)
    env_init <- initial_continuous_env(Z, y, Lenv, n_threads = n_threads, block_size = suff_block_size)
    ZCS <- env_init$ZCS
    pred0 <- mgcv_predictor_data(Z = ZCS, n = n)
    Data0 <- cbind(response_info$data, pred0)
    fit_final <- {
        gaussian_ridge_refit(y, pred0, numeric(0), dispersion = stats::var(y), n_threads = n_threads, block_size = suff_block_size)
    }
    g <- c()
    beta <- rep(0, p)
    alpha <- rep(0, ncol(Z))
    gamma <- 0
    beta_prev <- beta
    alpha_prev <- alpha
    eta0 <- fit_final$coefficients[1L]
    if (!is.finite(eta0))
        eta0 <- 0
    etaX <- 0
    etaW <- 0
    etaZ <- 0
    ZCS_refit <- NULL
    XCS_refit <- NULL
    WCS_refit <- NULL
    {
        warm <- greedy_lm_warm_start(X = X, y = y, Z = ZCS, L.init = L.init)
        if (length(warm$selected)) {
            warm_coef <- coef(warm$fit)
            warm_coef[is.na(warm_coef)] <- 0
            beta[warm$selected] <- utils::tail(warm_coef, length(warm$selected))
            etaX <- matrixVectorMultiply(X[, warm$selected, drop = FALSE], beta[warm$selected])
        }
    }
    fitX_no_cs_streak <- 0L
    for (iter in 1:max.iter) {
        beta_prev <- beta
        alpha_prev <- alpha
        phi0 <- fit_final$dispersion
        rZ <- y - etaX - etaW - eta0
        Zty <- as.vector(matrixMultiply(matrix(rZ, nrow = 1L), Z))
        yty4Z <- sum(rZ^2)
        fitZ <- .fit_susie_stage(structural = list(XtX = ZtZ, Xty = Zty, yty = yty4Z, n = n, L = Lenv),
            susie_para = susie_para_env, stage = "env", iter = iter, min.iter = min.iter, gaussian = TRUE, residual_variance = phi0)
        alpha <- coef(fitZ)[-1]
        etaZ <- matrixVectorMultiply(Z, alpha)
        env_terms <- build_select_env_terms(Z, fitZ, noncs_env_var = w_noncs_var, min_eta_env_var = min_etaW_var, noncs_max_abs_cor = noncs_max_abs_cor)
        if (is.null(env_terms$ZCS))
            stop_select_env_no_cs()
        ZCS <- env_terms$ZCS
        ZCS_refit <- env_terms$ZCS_refit
        rX <- y - etaZ - etaW - eta0
        Xty <- as.vector(matrixMultiply(matrix(rX, nrow = 1L), X))
        yty4X <- sum(rX^2)
        fitX <- .fit_susie_stage(structural = list(XtX = XtX, Xty = Xty, yty = yty4X, n = n, L = Lmain),
            susie_para = susie_para_main, stage = "main", iter = iter, min.iter = min.iter, gaussian = TRUE, residual_variance = phi0)
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
            rW <- y - etaZ - etaX - eta0
            WtW <- blockwise_crossprod(W, n_threads = n_threads,
                block_size = suff_block_size)
            Wty <- as.vector(matrixMultiply(matrix(rW, nrow = 1L), W))
            yty4W <- sum(rW^2)
            fitW <- .fit_susie_stage(structural = list(XtX = WtW, Xty = Wty, yty = yty4W, n = n, L = Lint),
                susie_para = susie_para_int, stage = "int", iter = iter, min.iter = min.iter, gaussian = TRUE, residual_variance = phi0)
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
            gaussian_ridge_refit(y, pred, penalty_V, phi0, n_threads = n_threads, block_size = suff_block_size)
        }
        coefs <- fit_final$coefficients
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
        {
            phi_err <- abs(fit_final$dispersion - phi0)/max(1, fit_final$dispersion, phi0)
            err <- max(err, phi_err)
        }
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
    refit_dispersion <- fit_final$dispersion
    penalty_names <- refit_penalty_terms(colnames(pred))
    penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names, fitZ = fitZ)
    {
        ridge_final <- gaussian_ridge_refit(y, pred, penalty_V, refit_dispersion, fixed_point = TRUE, n_threads = n_threads,
            block_size = suff_block_size)
        fit_final <- mgcv_fit_fixed_ridge(response_info$response, colnames(pred), Dat, family, penalty_V, dispersion = ridge_final$dispersion,
            mgcv_model = mgcv_model)
        attr(fit_final, "gaussian_ridge_fixed_point") <- list(iterations = ridge_final$fixed_point_iterations, dispersion = ridge_final$dispersion,
            penalty_dispersion = ridge_final$penalty_dispersion)
    }
    fit_final$n_eff <- n
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
