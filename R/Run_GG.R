Run_GG <- function(X, y, mgcv_model = NULL, crossprodX = NULL, Lmain, Lint, max.iter, min.iter, max.eps, susie_para_main,
    susie_para_int, verbose = TRUE, n_threads = 1, L.init = 1, x_noncs_var = 0.1, w_noncs_var = 0.1, noncs_max_abs_cor = 0.9,
    include_x_squared = FALSE, suff_block_size = 10000L, returnModel = FALSE) {
    family <- gaussian()
    run_start <- proc.time()[["elapsed"]]
    n <- NROW(y)
    p <- ncol(X)
    if (!is.null(crossprodX) && !identical(dim(as.matrix(crossprodX)), c(p, p))) {
        stop("crossprodX must be a square ncol(X)-by-ncol(X) matrix.")
    }
    suff_block_size <- validate_suff_block_size(suff_block_size)
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
    warm <- {
        warm_lm <- greedy_lm_warm_start(X = X, y = y, Z = NULL, L.init = L.init)
        fit0 <- list(coefficients = coef(warm_lm$fit), linear.predictors = stats::fitted(warm_lm$fit), dispersion = residual_variance_from_lm(warm_lm$fit))
        class(fit0) <- "gaussian_ridge_refit"
        list(fit = fit0, selected = warm_lm$selected)
    }
    fit_final <- warm$fit
    g <- c()
    err <- Inf
    beta <- rep(0, p)
    eta0 <- 0
    etaX <- 0
    etaW <- 0
    coef_init <- fit_final$coefficients
    coef_init[is.na(coef_init)] <- 0
    eta0 <- coef_init[1]
    if (length(warm$selected) > 0) {
        beta[warm$selected] <- coef_init[-1]
        etaX <- matrixVectorMultiply(X[, warm$selected, drop = FALSE], coef_init[-1])
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
        phi0 <- fit_final$dispersion
        rX <- y - etaW - eta0
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
        if (is.null(x_component$design)) {
            noncs_main <- build_full_noncs_refit_term(X, fitX)
            XCS <- NULL
            XCS_refit <- if (is.null(noncs_main))
                NULL
            else matrix(noncs_main, ncol = 1)
            if (!is.null(XCS_refit))
                colnames(XCS_refit) <- "Main_noncs_res"
            W <- NULL
            WCS <- NULL
            WCS_refit <- NULL
            fitW <- NULL
        }
        else {
            XCS <- x_component$design
            XCS_refit <- XCS
            main_noncs <- build_noncs_refit_term(X = X, fit = fitX, CSdt = CSdt, cs_indices = cs_indices, XCS = XCS, noncs_var = x_noncs_var,
                min_eta_var = min_etaX_var, noncs_max_abs_cor = noncs_max_abs_cor)
            XCS_refit <- append_noncs_refit_term(XCS_refit, main_noncs, "Main_noncs_res", corr_threshold = noncs_max_abs_cor)
            XCS_W <- XCS_refit
            W <- get_pairwise_interactions(XCS_W, include_x_squared = include_x_squared)
            WCS <- NULL
            WCS_refit <- NULL
            if (!interaction_design_available(W, iter, min.iter)) {
                W <- NULL
                fitW <- NULL
            }
            else {
                rW <- y - etaX - eta0
                WtW <- blockwise_crossprod(W, n_threads = n_threads,
                    block_size = suff_block_size)
                Wty <- as.vector(matrixMultiply(matrix(rW, nrow = 1L), W))
                yty4W <- sum(rW^2)
                fitW <- .fit_susie_stage(structural = list(XtX = WtW, Xty = Wty, yty = yty4W, n = n,
                  L = Lint), susie_para = susie_para_int, stage = "int", iter = iter, min.iter = min.iter, gaussian = TRUE,
                  residual_variance = phi0)
                CSdt_w <- summary(fitW)$vars
                cs_indices_w <- sort(unique(CSdt_w$cs[CSdt_w$cs > 0]))
                w_component <- build_component_design_from_fit(W, fitW, "Int_CS")
                WCS <- w_component$design
                WCS_refit <- WCS
                w_noncs <- build_w_noncs_refit_term(W = W, fitW = fitW, WCS = WCS, etaX = etaX, XCS = XCS, w_noncs_var = w_noncs_var,
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
            mgcv_predictor_data(Xextra = cbind(XCS_refit, WCS_refit), n = n)
        }
        else {
            mgcv_predictor_data(Xextra = XCS_refit, n = n)
        }
        Data <- cbind(response_info$data, pred)
        penalty_names <- refit_penalty_terms(colnames(pred))
        penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
        fit_final <- {
            gaussian_ridge_refit(y, pred, penalty_V, phi0, n_threads = n_threads, block_size = suff_block_size)
        }
        coefs <- fit_final$coefficients
        coefs[is.na(coefs)] <- 0
        eta0 <- coefs[1]
        idx <- 1
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
        err <- sqrt(mean((beta - beta_prev)^2))
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
    XCS_final <- XCS_refit
    WCS_final <- WCS_refit
    if (!is.null(XCS_final) && !is.null(WCS_final)) {
        pred <- mgcv_predictor_data(Xextra = cbind(XCS_final, WCS_final), n = n)
    }
    else if (!is.null(XCS_final)) {
        pred <- mgcv_predictor_data(Xextra = XCS_final, n = n)
    }
    else if (!is.null(WCS_final)) {
        pred <- mgcv_predictor_data(Xextra = WCS_final, n = n)
    }
    else {
        pred <- mgcv_predictor_data(n = n)
    }
    Dat <- cbind(response_info$data, pred)
    refit_dispersion <- {
        fit_final$dispersion
    }
    penalty_names <- refit_penalty_terms(colnames(pred))
    penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
    {
        ridge_final <- gaussian_ridge_refit(y, pred, penalty_V, refit_dispersion, fixed_point = TRUE, n_threads = n_threads,
            block_size = suff_block_size)
        fit_final <- mgcv_fit_fixed_ridge(response_info$response, colnames(pred), Dat, family, penalty_V, dispersion = ridge_final$dispersion,
            mgcv_model = mgcv_model)
        attr(fit_final, "gaussian_ridge_fixed_point") <- list(iterations = ridge_final$fixed_point_iterations, dispersion = ridge_final$dispersion,
            penalty_dispersion = ridge_final$penalty_dispersion)
    }
    fit_final$n_eff <- n
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
    AA <- list(diagnostics = make_diagnostics(iter, g, run_start), fitX = fitX, fitW = fitW, fitJoint = fit_final, main_discoveries = MainIndex,
        interaction_discoveries = IntIndex)
    if (returnModel)
        AA$FinalModel <- Dat
    AA$report <- extract_direction_table(AA, G)
    return(AA)
}
