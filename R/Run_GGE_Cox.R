Run_GGE_Cox <- function(X, Z, y, status,
                        Lmain, Lint, max.iter, min.iter, max.eps,
                        susie_para_main, susie_para_int,
                        noint_env = NULL,
                        verbose = TRUE, n_threads = 1,
                        L.init = 1,
                        x_noncs_var = 0.1,
                        w_noncs_var = 0.1,
                        noncs_max_abs_cor = 0.9,
                        include_x_squared = FALSE,
                        suff_block_size = NULL,
                        returnModel = FALSE) {

run_start <- proc.time()[["elapsed"]]
n <- length(y)
p <- ncol(X)
if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
nZ <- ncol(Z)
suff_block_size <- validate_suff_block_size(suff_block_size)
ridge <- 1e-8
eta_clip_range <- c(-50, 50)
min_etaX_var <- 1e-7
min_etaW_var <- 1e-7

surv_y <- survival::Surv(y, status)
warm <- greedy_cox_warm_start(X = X, y = y, status = status, Z = Z,
                              L.init = L.init)
fit_final <- warm$fit

g <- c()
beta <- rep(0, p)
alpha <- rep(0, nZ)
coef_init <- coef(fit_final)
coef_init[is.na(coef_init)] <- 0
if (length(coef_init) >= nZ) {
alpha <- coef_init[seq_len(nZ)]
}
if (length(warm$selected) > 0) {
init_beta <- coef_init[(nZ + 1):length(coef_init)]
beta[warm$selected] <- init_beta
}
XCS <- NULL
WCS <- NULL
WCS_refit <- NULL
W <- NULL
fitX <- NULL
fitW <- NULL
Data <- NULL

fitX_no_cs_streak <- 0L
for (iter in 1:max.iter) {
beta_prev <- beta
alpha_prev <- alpha

eta <- fit_final$linear.predictors
eta <- pmin(pmax(eta, eta_clip_range[1]), eta_clip_range[2])

ssX <- cox_suffstat_block(X, eta, cbind(Z, WCS_refit), y, status,
                          nuisance_precision = projection_penalty_precision(cbind(Z, WCS_refit), fitX, fitW),
                          n_threads = n_threads, ridge = ridge,
                           block_size = suff_block_size)
fitX <- .fit_susie_stage(
structural = list(XtX = ssX$XtX, Xty = ssX$Xty, yty = n - 1, n = n, L = Lmain),
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
} else {
XCS <- x_component$design
XCS_refit <- XCS
{
main_noncs <- build_noncs_refit_term(
X = X, fit = fitX, CSdt = CSdt, cs_indices = cs_indices,
XCS = XCS, noncs_var = x_noncs_var,
min_eta_var = min_etaX_var,
noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Z
)
XCS_refit <- append_noncs_refit_term(XCS_refit, main_noncs, "Main_noncs_res",
                                     corr_threshold = noncs_max_abs_cor)
}

Data <- as.data.frame(cbind(Z, XCS_refit, WCS_refit))
penalty_names <- refit_penalty_terms(names(Data))
penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
fit_final <- cox_fit_fixed_ridge(y, status, Data, penalty_V)
eta <- fit_final$linear.predictors
eta <- pmin(pmax(eta, eta_clip_range[1]), eta_clip_range[2])

XCS_W <- XCS_refit
W <- get_pairwise_interactions(XCS_W, Z = Z, noint_env = noint_env,
                               include_x_squared = if (is.null(XCS)) FALSE else include_x_squared)
WCS <- NULL
WCS_refit <- NULL
if (!interaction_design_available(W, iter, min.iter, allow_empty = main_no_cs)) {
W <- NULL
fitW <- NULL
} else {
ssW <- cox_suffstat_block(W, eta, cbind(Z, XCS_refit), y, status,
                          nuisance_precision = projection_penalty_precision(cbind(Z, XCS_refit), fitX, fitW),
                          n_threads = n_threads, ridge = ridge,
                           block_size = suff_block_size)
fitW <- .fit_susie_stage(
structural = list(XtX = ssW$XtX, Xty = ssW$Xty, yty = n - 1, n = n, L = Lint),
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
W = W, fitW = fitW, WCS = WCS, etaX = matrixVectorMultiply(X, beta),
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
}

if (!is.null(WCS_refit)) {
Data <- as.data.frame(cbind(Z, XCS_refit, WCS_refit))
} else {
Data <- as.data.frame(cbind(Z, XCS_refit))
}
penalty_names <- refit_penalty_terms(names(Data))
penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
fit_final <- cox_fit_fixed_ridge(y, status, Data, penalty_V)
alpha <- coef(fit_final)[seq_len(nZ)]
alpha[is.na(alpha)] <- 0

errX <- sqrt(mean((beta - beta_prev)^2))
errZ <- if (length(alpha)) sqrt(mean((alpha - alpha_prev)^2)) else 0
err <- max(errX, errZ)
g[iter] <- err
if (verbose) cat(sprintf("Iteration %d: err = %.3e, events = %d\n", iter, err, ssX$n_eff))
if (fitX_no_cs_streak >= 3L) {
  if (verbose) cat("No main credible set detected in 3 consecutive iterations; stopping.\n")
  break
}
if (err < max.eps && iter > min.iter) {
if (verbose) cat("Converged!\n")
break
}
}

XCS_final <- XCS_refit
WCS_final <- WCS_refit
if (!is.null(WCS_final)) {
Data <- as.data.frame(cbind(Z, XCS_final, WCS_final))
} else {
Data <- as.data.frame(cbind(Z, XCS_final))
}
penalty_names <- refit_penalty_terms(names(Data))
penalty_V <- refit_penalty_variance(fitX, fitW, penalty_names)
fit_final <- cox_fit_fixed_ridge(y, status, Data, penalty_V)
fit_final$n_eff <- ssX$n_eff

G <- cox_coef_table(fit_final)
MainIndex <- Identifying_MainEffect(fitX, colnames(X))
MainIndex <- safe_add_p(MainIndex, G)
IntIndex <- Identifying_IntEffect(fitW, colnames(W))
IntIndex <- filter_noncs_interactions(IntIndex)
IntIndex <- safe_add_p(IntIndex, G)

if (verbose) {
plot(g, type = "o", col = "black", pch = 16,
     xlab = "Iteration", ylab = "Max Parameter Change",
     main = "Convergence Trace (Cox PH, Breslow)")
for (i in seq_along(g)) {
text(x = i, y = g[i], labels = formatC(g[i], format = "e", digits = 1),
     pos = 3, cex = 0.7, col = "red")
}
}

AA <- list(diagnostics = make_diagnostics(iter, g, run_start),
           fitX = fitX,
           fitW = fitW,
           fitJoint = fit_final,
           main_discoveries = MainIndex,
           interaction_discoveries = IntIndex)
if (returnModel) AA$FinalModel <- Data
AA$report <- extract_direction_table(AA, G)
return(AA)
}
