Run_GGE_GLM <- function(X, Z, y, family = binomial(link = "logit"), mgcv_model = NULL, Lmain, Lint, max.iter, min.iter, max.eps,
    susie_para_main, susie_para_int, noint_env = NULL, verbose = TRUE, n_threads = 1, L.init = 1, x_noncs_var = 0.1, w_noncs_var = 0.1,
    noncs_max_abs_cor = 0.9, include_x_squared = FALSE, suff_block_size = 10000L, returnModel = FALSE) {
    run_start <- proc.time()[["elapsed"]]
    n <- NROW(y)
    p <- ncol(X)
    if (is.null(dim(Z)))
        Z <- matrix(Z, ncol = 1)
    if (is.null(colnames(Z)))
        colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    nZ <- ncol(Z)
    suff_block_size <- validate_suff_block_size(suff_block_size)
    validate_mgcv_irls_family(family)
    eta_clip_range <- c(-50, 50)
    weight_cutoff <- 0.0025
    min_etaX_var <- 1e-07
    min_etaW_var <- 1e-07
    response_info <- mgcv_prepare_response(y, family)
    if (response_info$n != nrow(X))
        stop("Length(y) must equal nrow(X).")
    warm <- {
        mgcv_greedy_warm_start(X = X, response_info = response_info, Z = Z, family = family, L.init = L.init, mgcv_model = mgcv_model)
    }
    fit_final <- warm$fit
    g <- c()
    err <- Inf
    beta <- rep(0, p)
    alpha <- rep(0, nZ)
    eta0 <- 0
    etaX <- 0
    etaW <- 0
    etaZ <- 0
    coef_init <- coef(fit_final)
    coef_init[is.na(coef_init)] <- 0
    eta0 <- coef_init[1]
    if (nZ > 0) {
        alpha <- coef_init[2:(1 + nZ)]
        etaZ <- matrixVectorMultiply(Z, alpha)
    }
    if (length(warm$selected) > 0) {
        init_beta <- coef_init[(2 + nZ):length(coef_init)]
        beta[warm$selected] <- init_beta
        etaX <- matrixVectorMultiply(X[, warm$selected, drop = FALSE], init_beta)
    }
    XCS <- NULL
    WCS <- NULL
    XCS_refit <- NULL
    WCS_refit <- NULL
    W <- NULL
    fitX <- NULL
    fitW <- NULL
    fitX_no_cs_streak <- 0L
    for (iter in 1:max.iter) {
        beta_prev <- beta
        alpha_prev <- alpha
        work <- {
            extract_mgcv_working(fit_final, weight_cutoff = weight_cutoff, eta_clip_range = eta_clip_range)
        }
        pseudo_response <- work$pseudo_response
        W_diag <- work$weights
        ZI_main <- cbind(Intercept = 1, Z, WCS_refit)
        ssX <- weighted_projected_suffstats(X = X, y = pseudo_response, ZI = ZI_main, weights = W_diag,
            nuisance_precision = projection_penalty_precision(ZI_main, fitX, fitW), n_threads = n_threads,
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
                main_noncs <- build_noncs_refit_term(X = X, fit = fitX, CSdt = CSdt, cs_indices = cs_indices, XCS = XCS,
                  noncs_var = x_noncs_var, min_eta_var = min_etaX_var, noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Z)
                XCS_refit <- append_noncs_refit_term(XCS_refit, main_noncs, "Main_noncs_res", corr_threshold = noncs_max_abs_cor)
            }
        }
        XCS_W <- XCS_refit
        W <- get_pairwise_interactions(XCS_W, Z = Z, noint_env = noint_env, include_x_squared = if (main_no_cs)
            FALSE
        else include_x_squared)
        WCS <- NULL
        WCS_refit <- NULL
        if (!interaction_design_available(W, iter, min.iter, allow_empty = main_no_cs)) {
            W <- NULL
            fitW <- NULL
        }
        else {
            ZI_int <- cbind(Intercept = 1, Z, XCS_refit)
            ssW <- weighted_projected_suffstats(W, pseudo_response, ZI_int, W_diag,
              nuisance_precision = projection_penalty_precision(ZI_int, fitX, fitW),
              n_threads = n_threads, block_size = suff_block_size)
            WtW <- ssW$XtX
            Wty <- ssW$Xty
            yty4W <- ssW$yty
            fitW <- .fit_susie_stage(structural = list(XtX = WtW, Xty = Wty, yty = yty4W, n = max(0.95 * n, work$n_eff), L = Lint),
                susie_para = susie_para_int, stage = "int", iter = iter, min.iter = min.iter, gaussian = FALSE, residual_variance = work$phi0)
            CSdt_w <- summary(fitW)$vars
            cs_indices_w <- sort(unique(CSdt_w$cs[CSdt_w$cs > 0]))
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
        pred <- if (!is.null(WCS_refit)) {
            mgcv_predictor_data(Z = Z, Xextra = cbind(XCS_refit, WCS_refit), n = n)
        }
        else {
            mgcv_predictor_data(Z = Z, Xextra = XCS_refit, n = n)
        }
        Data <- cbind(response_info$data, pred)
        penalty_names <- refit_penalty_terms(colnames(pred))
        penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
        fit_final <- {
            mgcv_fit_fixed_ridge(response_info$response, colnames(pred), Data, family, penalty_V, dispersion = work$phi0,
                mgcv_model = mgcv_model)
        }
        coefs <- coef(fit_final)
        coefs[is.na(coefs)] <- 0
        eta0 <- coefs[1]
        alpha <- coefs[2:(1 + nZ)]
        idx <- 1 + nZ
        etaZ <- matrixVectorMultiply(Z, alpha)
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
    XCS_final <- XCS_refit
    WCS_final <- WCS_refit
    if (!is.null(WCS_final)) {
        pred <- mgcv_predictor_data(Z = Z, Xextra = cbind(XCS_final, WCS_final), n = n)
    }
    else {
        pred <- mgcv_predictor_data(Z = Z, Xextra = XCS_final, n = n)
    }
    Dat <- cbind(response_info$data, pred)
    refit_dispersion <- mgcv_refit_dispersion(fit_final)
    penalty_names <- refit_penalty_terms(colnames(pred))
    penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
    {
        fit_final <- mgcv_fit_fixed_ridge(response_info$response, colnames(pred), Dat, family, penalty_V, dispersion = refit_dispersion,
            mgcv_model = mgcv_model)
    }
    fit_final$n_eff <- work$n_eff
    G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
    MainIndex <- Identifying_MainEffect(fitX, colnames(X))
    MainIndex <- safe_add_p(MainIndex, G)
    IntIndex <- Identifying_IntEffect(fitW, colnames(W))
    IntIndex <- filter_noncs_interactions(IntIndex)
    IntIndex <- safe_add_p(IntIndex, G)
    if (verbose) {
        plot(g, type = "o", col = "black", pch = 16, xlab = "Iteration", ylab = "Max Parameter Change", main = "Convergence Trace (GLM)")
        for (i in seq_along(g)) {
            text(x = i, y = g[i], labels = formatC(g[i], format = "e", digits = 1), pos = 3, cex = 0.7, col = "red")
        }
    }
    diagnostics <- make_diagnostics(iter, g, run_start)
    AA <- list(diagnostics = diagnostics, fitX = fitX, fitW = fitW, fitJoint = fit_final, main_discoveries = MainIndex, interaction_discoveries = IntIndex)
    if (returnModel)
        AA$FinalModel <- Dat
    AA$report <- extract_direction_table(AA, G)
    return(AA)
}
