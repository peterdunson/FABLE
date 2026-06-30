# Data generator. Block loadings, k=10, n=200, p=1000, 4 scenarios.

N         <- 200
base_seed <- 42
K         <- 10
P         <- 1000

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

# Example: one dataset from scenario 3 (lambda=5, Sigma=I)
# dat <- generate_data(param_sets[[3]], N, seed = base_seed)