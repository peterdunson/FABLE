## fable_analysis_k10_p1000.R
## Reads the FABLE simulation output and reports the same point-estimate,
## coverage, and efficiency summaries as the factorverse analysis script.

setwd("/Users/peterdunson/Desktop/FABLE")

library(ggplot2)
library(tidyr)
library(dplyr)

K <- 10
P <- 1000

fname <- file.path("test_fits", sprintf("sim_fable_%d_%d.rds", K, P))

results <- readRDS(fname)
results$method <- factor("fable")

set_label <- c("lambda=1, Sigma=I", "lambda=1, Sigma uneq",
               "lambda=5, Sigma=I", "lambda=5, Sigma uneq")

# ═══════════════════════════════════════════════════════════════════════════
#  PART 1: Point-estimate recovery (MSE on Omega and LambdaLambda', RV)
# ═══════════════════════════════════════════════════════════════════════════

cat("\n========== POINT-ESTIMATE RECOVERY ==========\n")
for (ps in sort(unique(results$param_set))) {
  r <- results[results$param_set == ps, ]
  if (nrow(r) == 0) next
  cat(sprintf("\n== Set %d [%s] ==\n", ps, set_label[ps]))
  cat(sprintf("  mse_Omega : mean %.4f / med %.4f / sd %.4f\n",
              mean(r$mse_Omega), median(r$mse_Omega), sd(r$mse_Omega)))
  cat(sprintf("  mse_LLt   : mean %.4f / med %.4f\n",
              mean(r$mse_LLt), median(r$mse_LLt)))
  cat(sprintf("  rv_LLt    : mean %.4f  (1 = perfect subspace recovery)\n",
              mean(r$rv_LLt)))
  cat(sprintf("  mse_sigma2: mean %.4f\n", mean(r$mse_sigma2)))
  if ("est_rank" %in% names(r))
    cat(sprintf("  k_hat     : %s  (true k = %d)\n",
                paste(sort(unique(r$est_rank)), collapse = ","), K))
}

# ═══════════════════════════════════════════════════════════════════════════
#  PART 2: Coverage of Omega (target 0.95)
# ═══════════════════════════════════════════════════════════════════════════

cat("\n========== COVERAGE OF OMEGA (target 0.95) ==========\n")
cover_tbl <- results |>
  group_by(param_set) |>
  summarise(
    cover_all     = mean(cover_Omega),
    cover_diag    = mean(cover_Omega_diag),
    cover_offdiag = mean(cover_Omega_offdiag, na.rm = TRUE),
    ci_width      = mean(ci_width_Omega),
    .groups = "drop"
  ) |>
  arrange(param_set) |>
  mutate(set = set_label[param_set])

print(as.data.frame(cover_tbl), digits = 3)

cover_long <- cover_tbl |>
  pivot_longer(c(cover_all, cover_diag, cover_offdiag),
               names_to = "target", values_to = "coverage") |>
  mutate(target = recode(target,
                         cover_all     = "All entries",
                         cover_diag    = "Diagonal (variances)",
                         cover_offdiag = "Off-diag (covariances)"))

p_cover <- ggplot(cover_long, aes(x = factor(param_set), y = coverage, fill = target)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Parameter set", y = "Empirical coverage", fill = NULL,
       title = "FABLE coverage of Omega 95% CIs (K=10, P=1000; dashed = 0.95)") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_cover)

# ═══════════════════════════════════════════════════════════════════════════
#  PART 3: MSE boxplots (Omega, LambdaLambda', Sigma)
# ═══════════════════════════════════════════════════════════════════════════

p_mse <- results |>
  pivot_longer(c(mse_Omega, mse_LLt, mse_sigma2),
               names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                         mse_Omega  = "MSE (Omega)",
                         mse_LLt    = "MSE (Lambda Lambda')",
                         mse_sigma2 = "MSE (Sigma)")) |>
  ggplot(aes(x = factor(param_set), y = value)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.4, fill = "steelblue") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Parameter set", y = NULL,
       title = "FABLE recovery MSE (K=10, P=1000, sparse blocks)") +
  theme_minimal()
print(p_mse)

# ═══════════════════════════════════════════════════════════════════════════
#  PART 4: Computational efficiency
# ═══════════════════════════════════════════════════════════════════════════

cat("\n========== COMPUTATIONAL EFFICIENCY ==========\n")
eff_tbl <- results |>
  group_by(param_set) |>
  summarise(
    time_mean = mean(time_sec, na.rm = TRUE),
    time_sd   = sd(time_sec, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(param_set) |>
  mutate(set = set_label[param_set])

print(as.data.frame(eff_tbl), digits = 4)
cat("\n(FABLE produces independent draws, so ESS / ESS-per-sec are not defined.)\n")