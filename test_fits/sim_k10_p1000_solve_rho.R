# FABLE coverage simulation, block loadings, k=10, n=200, p=1000.
# Adds a direct solve for rho^2: pick rho so the AVERAGE entrywise asymptotic
# coverage equals 1-alpha exactly (FABLE eq. 16), instead of the default
# average-b approximation. Set RHO_MODE to compare.

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

# "default"  = average-b approximation, varInflation = (sum(CC)/(p(p+1)/2))^2
# "solve"    = root-find rho so mean entrywise coverage = 1-alpha (eq. 16)
RHO_MODE  <- "solve"

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

# ── Direct rho solve (FABLE Section 3.3 / eq. 16) ─────────────────────────────
# Plug-in entrywise coverage q_uv(rho) = 2*Phi(z * l0_uv(rho)/S0_uv) - 1, with
#   off-diag: l0^2 = rho^2 (sig_v ||mu_u||^2 + sig_u ||mu_v||^2)
#             S0^2 = sig_v ||mu_u||^2 + sig_u ||mu_v||^2 + ||mu_u||^2 ||mu_v||^2 + (mu_u'mu_v)^2
#   diag:     l0^2 = 2 sig_u^2 + 4 rho^2 sig_u ||mu_u||^2
#             S0^2 = 2 (||mu_u||^2 + sig_u)^2
# Plug-ins from FABLEHyperParameters: ||mu_u||^2 = G[u,u], mu_u'mu_v = G[u,v],
# sig_u = SigmaSqEstimate[u].
# Solve mean_uv q_uv(rho) = 1-alpha for rho >= 1 by uniroot.
# Vectorised over the upper triangle (incl diagonal) for speed at p=1000.

solve_rho_avg <- function(G, sigsq, alpha = 0.05, rho_max = 50) {
  P  <- nrow(G)
  z  <- qnorm(1 - alpha / 2)
  d  <- diag(G)                       # ||mu_u||^2 for each u
  
  ut <- which(upper.tri(G, diag = TRUE), arr.ind = TRUE)
  u  <- ut[, 1]; v <- ut[, 2]
  is_diag <- (u == v)
  
  du <- d[u]; dv <- d[v]
  su <- sigsq[u]; sv <- sigsq[v]
  guv <- G[cbind(u, v)]
  
  # rho-independent pieces
  # off-diagonal
  A_off <- sv * du + su * dv                                  # multiplies rho^2 in l0^2
  S2_off <- sv * du + su * dv + du * dv + guv^2
  # diagonal
  c_diag <- 2 * su^2                                          # rho-independent part of l0^2
  A_diag <- 4 * su * du                                       # multiplies rho^2 in l0^2
  S2_diag <- 2 * (du + su)^2
  
  S2 <- ifelse(is_diag, S2_diag, S2_off)
  S  <- sqrt(S2)
  
  # mean coverage as a function of rho
  mean_cov <- function(rho) {
    l0sq <- ifelse(is_diag,
                   c_diag + (rho^2) * A_diag,
                   (rho^2) * A_off)
    q <- 2 * pnorm(z * sqrt(l0sq) / S) - 1
    mean(q)
  }
  
  # mean_cov is increasing in rho; find root of mean_cov(rho) - (1-alpha)
  target <- 1 - alpha
  f <- function(rho) mean_cov(rho) - target
  
  # at rho = 1 coverage may already exceed target (then rho = 1 is the floor)
  if (f(1) >= 0) return(1)
  if (f(rho_max) < 0) return(rho_max)   # cap if even rho_max undershoots
  uniroot(f, interval = c(1, rho_max), tol = 1e-6)$root
}

# ── Fixed-rank FABLE fit ──────────────────────────────────────────────────────

fit_fable_fixedk <- function(dat, k_fixed = K, MC = 1000, alpha = 0.05,
                             gamma0 = 1, delta0sq = 1, rho_mode = "solve") {
  
  Y <- as.matrix(dat)
  P <- ncol(Y)
  
  t <- system.time({
    svdY   <- svd(Y)
    U_Y    <- svdY$u
    V_Y    <- svdY$v
    svalsY <- svdY$d
    
    hyp <- FABLEHyperParameters(Y, U_Y, V_Y, svalsY, kEst = k_fixed)
    
    if (rho_mode == "default") {
      CCMatrix     <- cov_correct_matrix(hyp$SigmaSqEstimate, hyp$G)
      varInflation <- (sum(CCMatrix) / (P * (P + 1) / 2))^2
    } else {
      rho_star     <- solve_rho_avg(hyp$G, hyp$SigmaSqEstimate, alpha = alpha)
      varInflation <- rho_star^2
    }
    
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
    est_rank  = samp$EstimatedRank,
    var_infl  = varInflation,
    ess_mean    = NA_real_,
    ess_per_sec = NA_real_
  )
}

run_simulation_fable <- function(param_sets, N, B = 100, base_seed = 42,
                                 k_fixed = K, MC = 1000, alpha = 0.05,
                                 rho_mode = "solve") {
  
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
      result <- fit_fable_fixedk(dat, k_fixed = k_fixed, MC = MC,
                                 alpha = alpha, rho_mode = rho_mode)
      
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
        aabias_Omega        = mean(abs(diff_Omega)),
        mabias_Omega        = max(abs(diff_Omega)),
        mse_LLt             = mean(diff_LLt^2),
        aabias_LLt          = mean(abs(diff_LLt)),
        mabias_LLt          = max(abs(diff_LLt)),
        mse_sigma2          = mean(diff_sigma2^2),
        aabias_sigma2       = mean(abs(diff_sigma2)),
        mabias_sigma2       = max(abs(diff_sigma2)),
        rv_LLt              = rv_LLt,
        cover_Omega         = mean(covered[ut]),
        cover_Omega_diag    = mean(covered[dg]),
        cover_Omega_offdiag = if (any(offdiag)) mean(covered[offdiag]) else NA,
        ci_width_Omega      = mean(width[ut]),
        time_sec            = result$time_sec,
        ess_mean            = result$ess_mean,
        ess_per_sec         = result$ess_per_sec,
        est_rank            = result$est_rank,
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
  alpha      = alpha,
  rho_mode   = RHO_MODE
)

results_fable$method <- paste0("fable_", RHO_MODE)

fname <- file.path(out_dir, sprintf("sim_fable_%s_%d_%d.rds", RHO_MODE, K, P))
saveRDS(results_fable, fname)

set_label <- c("lambda=1, Sigma=I", "lambda=1, Sigma uneq",
               "lambda=5, Sigma=I", "lambda=5, Sigma uneq")

cat(sprintf("\nrho_mode = %s\n", RHO_MODE))
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
