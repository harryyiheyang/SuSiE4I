# Run_GGE_GLM with the linear Z replaced by a GAM null model (see SuSiE4I_GAM):
# the null smooths enter every projection as B with precision S_lambda / phi
# (Vp / phi of the GAM); interactions are built from Main_CS only.
Run_GAM <- function(X, null, family = gaussian(), mgcv_model = NULL, Lmain, Lint, max.iter, min.iter, max.eps,
    susie_para_main, susie_para_int, verbose = TRUE, n_threads = 1, x_noncs_var = 0.1, w_noncs_var = 0.1,
    noncs_max_abs_cor = 0.9, include_x_squared = FALSE, suff_block_size = 10000L, int_suggested_coverage = NULL,
    returnModel = FALSE) {
    run_start <- proc.time()[["elapsed"]]
    n <- nrow(X)
    p <- ncol(X)
    Z <- null$Zint
    B <- null$B
    eta_clip_range <- c(-50, 50)
    weight_cutoff <- 0.0025
    min_etaX_var <- 1e-07
    min_etaW_var <- 1e-07
    fit_final <- null$fit
    g <- c()
    err <- Inf
    beta <- rep(0, p)
    XCS <- NULL
    WCS <- NULL
    XCS_refit <- NULL
    WCS_refit <- NULL
    W <- NULL
    fitX <- NULL
    fitW <- NULL
    int_smooth <- NULL
    fml_refit <- null$formula
    fitX_no_cs_streak <- 0L
    for (iter in 1:max.iter) {
        beta_prev <- beta
        work <- {
            extract_mgcv_working(fit_final, weight_cutoff = weight_cutoff, eta_clip_range = eta_clip_range)
        }
        pseudo_response <- work$pseudo_response
        W_diag <- work$weights
        S_null <- gam_null_penalty(null$fit, fit_final) / work$phi0
        Bm <- B
        Sm <- S_null
        if (!is.null(int_smooth)) {
            # s(z, by = Main_CS) terms of the last refit are projected with their REML penalty
            Lp <- stats::predict(fit_final, type = "lpmatrix")
            keep <- !colnames(Lp) %in% names(attr(fit_final, "refit_penalty")$V)
            Bm <- Lp[, keep, drop = FALSE]
            Sm <- (gam_null_penalty(fit_final, fit_final) / work$phi0)[keep, keep, drop = FALSE]
        }
        ZI_main <- cbind(Bm, WCS_refit)
        P_main <- diag(c(rep(0, ncol(Bm)), projection_penalty_precision(WCS_refit, fitX, fitW)), ncol(ZI_main))
        P_main[seq_len(ncol(Bm)), seq_len(ncol(Bm))] <- Sm
        ssX <- weighted_projected_suffstats(X = X, y = pseudo_response, ZI = ZI_main, weights = W_diag,
            nuisance_precision = P_main, n_threads = n_threads,
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
        W <- if (main_no_cs) NULL else get_pairwise_interactions(XCS, Z = Z, include_x_squared = include_x_squared)
        WCS <- NULL
        WCS_refit <- NULL
        int_smooth <- NULL
        if (!interaction_design_available(W, iter, min.iter, allow_empty = main_no_cs)) {
            W <- NULL
            fitW <- NULL
        }
        else {
            ZI_int <- cbind(B, XCS_refit)
            P_int <- diag(c(rep(0, ncol(B)), projection_penalty_precision(XCS_refit, fitX, fitW)), ncol(ZI_int))
            P_int[seq_len(ncol(B)), seq_len(ncol(B))] <- S_null
            ssW <- weighted_projected_suffstats(W, pseudo_response, ZI_int, W_diag,
              nuisance_precision = P_int,
              n_threads = n_threads, block_size = suff_block_size)
            WtW <- ssW$XtX
            Wty <- ssW$Xty
            yty4W <- ssW$yty
            fitW <- .fit_susie_stage(structural = list(XtX = WtW, Xty = Wty, yty = yty4W, n = max(0.95 * n, work$n_eff), L = Lint),
                susie_para = susie_para_int, stage = "int", suggested_coverage = int_suggested_coverage, iter = iter, min.iter = min.iter, gaussian = FALSE, residual_variance = work$phi0)
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
            # An Int CS whose z (z or f(z)) is smooth in the null model is refitted as s(z, by = Xi).
            csw <- susie_cs_list(fitW)
            zw <- sub("^f\\((.*)\\)$", "\\1", sub("\\*Main_CS[0-9]+$", "", colnames(W)))
            for (k in seq_along(csw$index)) {
                v <- csw$vars[[k]][zw[csw$vars[[k]]] %in% names(environment(null$formula)$.s4i_ibases)]
                if (length(v)) {
                    v <- v[which.max(fitW$alpha[csw$index[k], v])]
                    int_smooth[paste0("Int_CS", csw$index[k])] <- paste0("s(", zw[v], "):", sub("^.*\\*", "", colnames(W)[v]))
                }
            }
            if (!is.null(int_smooth)) {
                WCS_refit <- WCS_refit[, !colnames(WCS_refit) %in% names(int_smooth), drop = FALSE]
                if (!ncol(WCS_refit)) WCS_refit <- NULL
            }
        }
        fml_refit <- null$formula
        if (!is.null(int_smooth)) {
            zs <- sub("^s\\((.*)\\):.*$", "\\1", int_smooth)
            xs <- sub("^.*:", "", int_smooth)
            fml_refit <- stats::update(null$formula, stats::as.formula(paste(". ~ . +", paste(unique(sprintf(
                "s(%s, by = %s, bs = \"s4iAM\", xt = list(base = .s4i_ibases[[\"%s\"]]))", zs, xs, zs)), collapse = " + "))))
        }
        pred <- if (!is.null(WCS_refit)) {
            mgcv_predictor_data(Xextra = cbind(XCS_refit, WCS_refit), n = n)
        }
        else {
            mgcv_predictor_data(Xextra = XCS_refit, n = n)
        }
        Data <- cbind(null$data, pred)
        penalty_names <- refit_penalty_terms(colnames(pred))
        penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
        fit_final <- {
            mgcv_fit_fixed_ridge(null$response, colnames(pred), Data, family, penalty_V, dispersion = work$phi0,
                mgcv_model = mgcv_model, formula = fml_refit)
        }
        coefs <- coef(fit_final)
        coefs[is.na(coefs)] <- 0
        if (!is.null(XCS_refit)) {
            etaX <- matrixVectorMultiply(XCS_refit, coefs[colnames(XCS_refit)])
        }
        else {
            etaX <- rep(0, n)
        }
        if (!is.null(WCS_refit)) {
            etaW <- matrixVectorMultiply(WCS_refit, coefs[colnames(WCS_refit)])
        }
        else {
            etaW <- 0
        }
        err <- sqrt(mean((beta - beta_prev)^2))
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
        pred <- mgcv_predictor_data(Xextra = cbind(XCS_final, WCS_final), n = n)
    }
    else {
        pred <- mgcv_predictor_data(Xextra = XCS_final, n = n)
    }
    Dat <- cbind(null$data, pred)
    refit_dispersion <- mgcv_refit_dispersion(fit_final)
    penalty_names <- refit_penalty_terms(colnames(pred))
    penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
    {
        fit_final <- mgcv_fit_fixed_ridge(null$response, colnames(pred), Dat, family, penalty_V, dispersion = refit_dispersion,
            mgcv_model = mgcv_model, formula = fml_refit)
    }
    fit_final$n_eff <- work$n_eff
    G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
    if (!is.null(int_smooth)) {
        St <- summary(fit_final)$s.table[int_smooth, , drop = FALSE]
        St[, 1:3] <- NA
        rownames(St) <- names(int_smooth)
        G <- rbind(G, St)
    }
    MainIndex <- Identifying_MainEffect(fitX, colnames(X))
    MainIndex <- safe_add_p(MainIndex, G)
    IntIndex <- Identifying_IntEffect(fitW, colnames(W))
    IntIndex <- filter_noncs_interactions(IntIndex)
    IntIndex <- safe_add_p(IntIndex, G)
    if (verbose) {
        plot(g, type = "o", col = "black", pch = 16, xlab = "Iteration", ylab = "Max Parameter Change", main = "Convergence Trace (GAM)")
        for (i in seq_along(g)) {
            text(x = i, y = g[i], labels = formatC(g[i], format = "e", digits = 1), pos = 3, cex = 0.7, col = "red")
        }
    }
    diagnostics <- make_diagnostics(iter, g, run_start)
    AA <- list(diagnostics = diagnostics, fitNull = null$fit, fitX = fitX, fitW = fitW, fitJoint = fit_final, main_discoveries = MainIndex, interaction_discoveries = IntIndex)
    if (returnModel)
        AA$FinalModel <- Dat
    AA$report <- extract_direction_table(AA, G)
    return(AA)
}
