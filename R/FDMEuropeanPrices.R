
# Explicit finite-difference scheme for pricing European put options

# TODO: Update to include call/put

#' European put price using FDM
#'
#' @param sigma Volatility
#' @param r Risk-free interest rate
#' @param K Strike price
#' @param Smax Truncation of asset price
#' @param N Number of spatial intervals
#' @param M Number of time-steps
#'
#' @return Vector of prices e.g. V[1,] give the prices at time 0 for each S
#' @export
explicitFDMEuropeanPut <- function(sigma, r, T, K, Smax, N, M) {

  dt <- T/M
  dS <- Smax/N

  A <- numeric(N + 1)
  B <- numeric(N + 1)
  C <- numeric(N + 1)
  V <- matrix(0, nrow = M + 1, ncol = N + 1)

  # Compute coefficients
  for (n in 0:N) {
    A[n + 1] <- 0.5 * n^2 * sigma^2 * dt - 0.5 * n * r * dt
    B[n + 1] <- 1 - n^2 * sigma^2 * dt - r * dt
    C[n + 1] <- 0.5 * n^2 * sigma^2 * dt + 0.5 * n * r * dt
  }

  # Terminal conditions
  for (n in 0:N) {
    V[M + 1, n + 1] <- max(K - n * dS, 0)
  }

  # Time stepping
  for (m in M:1) {
    V[m, 1] <- B[1] * V[m + 1, 1] + C[1] * V[m + 1, 2]
    for (n in 1:(N - 1)) {
      V[m, n + 1] <- A[n + 1] * V[m + 1, n] + B[n + 1] * V[m + 1, n + 1]
              + C[n + 1] * V[m + 1, n + 2]
    }
    V[m, N + 1] <- A[N + 1] * V[m + 1, N] + B[N + 1] * V[m + 1, N + 1]
  }

  # Price matrix
  return(V)
}




