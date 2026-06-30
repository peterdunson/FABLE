## fable_sim_k10_p1000.R
## Same structure and seeds as the k10_p200 FABLE run, but P=1000 with n=200.
## Purpose: test whether the ~99.9% overcoverage at n=p=200 moves toward the
## nominal 0.95 as sqrt(n)/p shrinks.
##   n=p=200  : sqrt(n)/p = 0.071
##   n=200, p=1000 : sqrt(n)/p = 0.014   (5x closer to the large-p regime)
##
## FIXED RANK k=10 (bypasses FABLE's internal JIC, which overestimates k here).

setwd("/Users/peterdunson/Desktop/FABLE")

library(FABLE)

out_dir   <- "test_fits"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N         <- 200
B         <- 10
base_seed <- 42
K         <- 10        # fixed number of factors
P         <- 1000      # <-- only change from the p200 run
MC        <- 1000
alpha     <- 0.05

# ── Same DGP, extended to P=1000 ──────────────────────────────────────────────
# Block size becomes 1000/10 = 100. sigma_unequal recycles the same
# c(0.1,0.1,0.1,0.5,0.5) pattern. Identical to the p200 setup otherwise.

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

# ── Fixed-rank FABLE fit (identical to the p200 script) ───────────────────────

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
    
    CCMatrix     <- cov_correct_matrix(hyp$SigmaSqEstimate, hyp$G)
    varInflation <- (sum(CCMatrix) / (P * (P + 1) / 2))^2
    
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
    var_infl  = varInflation,        # track rho^2 to see if it shrinks vs p200
    ess_mean    = NA_real_,
    ess_per_sec = NA_real_
  )
}

# ── Simulation loop ───────────────────────────────────────────────────────────

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
    
    cat(sprintf("\n  Param set %d/%d\n", p_idx, length(param_sets)))
    flush.console()
    
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
      
      cat(sprintf("    rep %3d/%d  [%.1fs, k=%d, rho2=%.2f]\n",
                  b, B, result$time_sec, result$est_rank, result$var_infl))
      flush.console()
    }
  }
  
  do.call(rbind, results)
}

# ── Run and save ──────────────────────────────────────────────────────────────

cat(sprintf("\n── FABLE (fixed k=%d) | K=%d, P=%d, n=%d ──\n", K, K, P, N))
cat(sprintf("   sqrt(n)/p = %.4f\n", sqrt(N) / P))
flush.console()

results_fable <- run_simulation_fable(
  param_sets = param_sets,
  N          = N,
  B          = B,
  base_seed  = base_seed,
  k_fixed    = K,
  MC         = MC,
  alpha      = alpha
)

results_fable$method <- "fable"

fname <- file.path(out_dir, sprintf("sim_fable_%d_%d.rds", K, P))
saveRDS(results_fable, fname)
cat(sprintf("\nSaved: %s\n", fname))