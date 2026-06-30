# FABLE coverage simulation, block loadings, k=10, n=200, p=1000.
# rho solved directly so average entrywise coverage = 1-alpha (FABLE eq. 16),
# instead of the default average-b approximation. k forced to 10 via the C++.

setwd("/Users/peterdunson/Desktop/FABLE")

library(FABLE)

out_dir   <- "test_fits"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N         <- 200
B         <- 5
base_seed <- 42
K         <- 10
P         <- 1000
MC        <- 1000
alpha     <- 0.05

make_block_lambda <- function(P, K, val) {
  Lambda <- matrix(0, nrow = P, ncol = K)
  block  <- P / K
  for (h in seq_len(K)) {
    rows <- ((h - 1) * block + 1):(h * block)
    Lambda[rows, h] <- val
  }
  Lambda
}

sigma_unequal <- rep(c(0.1, 0.1, 0.1, 0.5, 0.5), length.out = P)
sigma_equal   <- rep(1, P)

param_sets <- list(
  list(Lambda = make_block_lambda(P, K, 1), Sigma = diag(sigma_equal)),
  list(Lambda = make_block_lambda(P, K, 1), Sigma = diag(sigma_unequal)),
  list(Lambda = make_block_lambda(P, K, 5), Sigma = diag(sigma_equal)),
  list(Lambda = make_block_lambda(P, K, 5), Sigma = diag(sigma_unequal))
)
param_sets <- lapply(param_sets, function(p) {
  p$Omega <- tcrossprod(p$Lambda) + p$Sigma
  p
})

generate_data <- function(params, N, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  L <- chol(params$Omega)
  Z <- matrix(rnorm(N * nrow(params$Lambda)), nrow = N)
  Z %*% L
}

# Solve rho so mean entrywise asymptotic coverage = 1-alpha (FABLE eq. 16).
# Plug-ins from FABLEHyperParameters: ||mu_u||^2 = G[u,u], mu_u'mu_v = G[u,v],
# sig_u = SigmaSqEstimate[u].
solve_rho_avg <- function(G, sigsq, alpha = 0.05, rho_max = 50) {
  z  <- qnorm(1 - alpha / 2)
  d  <- diag(G)
  
  ut <- which(upper.tri(G, diag = TRUE), arr.ind = TRUE)
  u  <- ut[, 1]; v <- ut[, 2]
  is_diag <- (u == v)
  
  du <- d[u]; dv <- d[v]
  su <- sigsq[u]; sv <- sigsq[v]
  guv <- G[cbind(u, v)]
  
  A_off  <- sv * du + su * dv
  S2_off <- sv * du + su * dv + du * dv + guv^2
  c_diag <- 2 * su^2
  A_diag <- 4 * su * du
  S2_diag <- 2 * (du + su)^2
  
  S <- sqrt(ifelse(is_diag, S2_diag, S2_off))
  
  mean_cov <- function(rho) {
    l0sq <- ifelse(is_diag, c_diag + (rho^2) * A_diag, (rho^2) * A_off)
    mean(2 * pnorm(z * sqrt(l0sq) / S) - 1)
  }
  
  f <- function(rho) mean_cov(rho) - (1 - alpha)
  if (f(1) >= 0) return(1)
  if (f(rho_max) < 0) return(rho_max)
  uniroot(f, interval = c(1, rho_max), tol = 1e-6)$root
}

fit_fable_fixedk <- function(dat, k_fixed = K, MC = 1000, alpha = 0.05,
                             gamma0 = 1, delta0sq = 1) {
  
  Y <- as.matrix(dat)
  P <- ncol(Y)
  
  t <- system.time({
    svdY   <- svd(Y)
    U_Y    <- svdY$u
    V_Y    <- svdY$v
    svalsY <- svdY$d
    
    hyp <- FABLEHyperParameters(Y, U_Y, V_Y, svalsY, kEst = k_fixed)
    
    rho_star     <- solve_rho_avg(hyp$G, hyp$SigmaSqEstimate, alpha = alpha)
    varInflation <- rho_star^2
    
    samp <- CPPFABLESampler(Y, gamma0, delta0sq, MC,
                            U_Y, V_Y, svalsY,
                            kEst = k_fixed,
                            varInflation = varInflation)
    
    pp <- CPPCCFABLEPostProcessing(samp, alpha)
  })
  time_sec <- unname(t["elapsed"])
  
  LLt_hat   <- samp$G
  Sigma_hat <- as.numeric(samp$SigmaSqEstimatePostMean)
  Omega_hat <- LLt_hat + diag(Sigma_hat)
  
  list(
    LLt_hat   = LLt_hat,
    Sigma_hat = Sigma_hat,
    Omega_hat = Omega_hat,
    Omega_lo  = pp$LowerQuantileMatrix,
    Omega_hi  = pp$UpperQuantileMatrix,
    time_sec  = time_sec,
    var_infl  = varInflation
  )
}

run_simulation_fable <- function(param_sets, N, B = 100, base_seed = 42,
                                 k_fixed = K, MC = 1000, alpha = 0.05) {
  
  results <- vector("list", length(param_sets) * B)
  idx <- 1L
  
  for (p_idx in seq_along(param_sets)) {
    params     <- param_sets[[p_idx]]
    Omega_true <- params$Omega
    Pp         <- nrow(Omega_true)
    ut         <- upper.tri(Omega_true, diag = TRUE)
    offdiag    <- upper.tri(Omega_true, diag = FALSE)
    dg         <- diag(TRUE, Pp)
    LLt_true   <- tcrossprod(params$Lambda)
    
    for (b in seq_len(B)) {
      dat    <- generate_data(params, N, seed = base_seed + p_idx * B + b)
      result <- fit_fable_fixedk(dat, k_fixed = k_fixed, MC = MC, alpha = alpha)
      
      diff_Omega  <- result$Omega_hat - Omega_true
      diff_LLt    <- result$LLt_hat   - LLt_true
      diff_sigma2 <- result$Sigma_hat - diag(params$Sigma)
      
      LLt_hat <- result$LLt_hat
      denom   <- sqrt(sum(LLt_hat^2) * sum(LLt_true^2))
      rv_LLt  <- if (denom > 0) sum(LLt_hat * LLt_true) / denom else 0
      
      covered <- (result$Omega_lo <= Omega_true) & (Omega_true <= result$Omega_hi)
      width   <- result$Omega_hi - result$Omega_lo
      
      results[[idx]] <- data.frame(
        param_set           = p_idx, b = b,
        mse_Omega           = mean(diff_Omega^2),
        mse_LLt             = mean(diff_LLt^2),
        mse_sigma2          = mean(diff_sigma2^2),
        rv_LLt              = rv_LLt,
        cover_Omega         = mean(covered[ut]),
        cover_Omega_diag    = mean(covered[dg]),
        cover_Omega_offdiag = if (any(offdiag)) mean(covered[offdiag]) else NA,
        ci_width_Omega      = mean(width[ut]),
        time_sec            = result$time_sec,
        var_infl            = result$var_infl
      )
      idx <- idx + 1L
    }
  }
  
  do.call(rbind, results)
}

results_fable <- run_simulation_fable(
  param_sets = param_sets,
  N          = N,
  B          = B,
  base_seed  = base_seed,
  k_fixed    = K,
  MC         = MC,
  alpha      = alpha
)

fname <- file.path(out_dir, sprintf("sim_fable_solve_%d_%d.rds", K, P))
saveRDS(results_fable, fname)

set_label <- c("lambda=1, Sigma=I", "lambda=1, Sigma uneq",
               "lambda=5, Sigma=I", "lambda=5, Sigma uneq")

for (ps in sort(unique(results_fable$param_set))) {
  r <- results_fable[results_fable$param_set == ps, ]
  print(data.frame(
    set       = set_label[ps],
    cover     = round(mean(r$cover_Omega), 3),
    ci_width  = round(mean(r$ci_width_Omega), 3),
    rho2      = round(mean(r$var_infl), 2),
    mse_Omega = round(mean(r$mse_Omega), 4)
  ), row.names = FALSE)
}