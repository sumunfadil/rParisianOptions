
#' Explicit finite-difference scheme for pricing
#' European put options
#'
#' @param sigma Volatility
#' @param r Risk-free interest rate
#' @param K Strike price
#' @param Smax Truncation of asset price
#' @param N Number of spatial intervals
#' @param M Number of time-steps
#'
#' @return Vector of prices V. V[m,] give the prices
#' at time-step m for each S
#' @export
explicitFDMEuropeanPut <- function(sigma, r, T, K, Smax, N, M) {

  delta_t <- T / M
  delta_S <- Smax / N

  # Initialize vectors
  A <- numeric(N + 1)
  B <- numeric(N + 1)
  C <- numeric(N + 1)
  V <- matrix(0, nrow = M + 1, ncol = N + 1)

  # Compute coefficients
  for (n in 0:N) {
    A[n + 1] <- 0.5 * n^2 * sigma^2 * delta_t - 0.5 * n * r * delta_t
    B[n + 1] <- 1 - n^2 * sigma^2 * delta_t - r * delta_t
    C[n + 1] <- 0.5 * n^2 * sigma^2 * delta_t + 0.5 * n * r * delta_t
  }

  # Terminal conditions
  for (n in 0:N) {
    V[M + 1, n + 1] <- max(K - n * delta_S, 0)
  }

  # Time stepping
  for (m in M:1) {
    V[m, 1] <- B[1] * V[m + 1, 1] + C[1] * V[m + 1, 2]
    for (n in 1:(N - 1)) {
      V[m, n + 1] <- A[n + 1] * V[m + 1, n] + B[n + 1] * V[m + 1, n + 1] + C[n + 1] * V[m + 1, n + 2]
    }
    V[m, N + 1] <- A[N + 1] * V[m + 1, N] + B[N + 1] * V[m + 1, N + 1]
  }
  return(V)
}


#' Implicit finite-difference scheme for pricing
#' European put options
#'
#' @param sigma Volatility
#' @param r Risk-free interest rate
#' @param K Strike price
#' @param Smax Truncation of asset price
#' @param N Number of spatial intervals
#' @param M Number of time-steps
#'
#' @return Vector of prices V. V[m,] give the prices
#' at time-step m for each S
#' @export
implicitFDMEuropeanPut <- function(sigma, r, T, K, Smax, N, M) {

  delta_t <- T/M
  delta_S <- Smax/N

  # Initialize vectors
  a <- numeric(N-1)
  b <- numeric(N)
  c <- numeric(N-1)
  KMatrix <- matrix(0, nrow = N, ncol = N)
  V <- matrix(0, nrow = M + 1, ncol = N + 1)

  # Compute coefficients
  for (n in 0:(N-1)) {

    if (n>0) {
      a[n] <- -0.5 * n^2 * sigma^2 * delta_t + 0.5 * n * r * delta_t
    }

    b[n+1] <- 1 + n^2 * sigma^2 * delta_t + r * delta_t

    if (n<(N-1)) {
      c[n+1] <- -0.5 * n^2 * sigma^2 * delta_t - 0.5 * n * r * delta_t
    }
  }

  KMatrix <- tridiagonalMatrix(a, b, c)

  # Terminal conditions
  for (n in 0:N) {
    V[M + 1, n + 1] <- max(K - n * delta_S, 0)
  }

  # VNm = 0 (by default zero entries)

  # Matrix multiplication
  for (m in (M+1):1) {
    V[m-1,1:N] = thomasAlgorithmSolver(KMatrix, V[m,1:N])
  }

  return(V)

}



