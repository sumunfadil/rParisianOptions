

#' Plot of grid points
#'
#' @return Grid points
#' @export
plotGridPoints <- function() {

  S <- seq(from = 0, to = 1, length.out = 9)
  t <- seq(from = 0, to = 1, length.out = 5)
  plot(NA, xlim = range(S), ylim = range(t), xlab = "S", ylab = "t")
  points(rep(S, each = length(t)), rep(t, times = length(S)), pch = 15)

}

#' Plotting function for FDM solution
#'
#' @param S Asset price
#' @param V Matrix of FDM values
#' @param BS Black-Scholes prices
#' @param m Time-step
#'
#' @return A plot for the solution of FDM at grid points
#' compared with Black-Scholes price
#' @export
#'
#' @examples
#' plotFDMV(S, V, BSPutPrice)
plotFDMV <- function(S, V, BS, m=1) {

  par(mfrow = c(1, 1))
  par(mar = c(4, 4.5, 4.5, 1) + 0.1)

  # V[m,]: price at time-step m (note: indexing starts at 1)
  plot(S, V[1,], type = "p", pch = 15, cex=1,
       xlab = expression(S ~ "," ~ S[n]),
       ylab = expression(V(S,0) ~ "," ~ V[n]^0))
  lines(S, BS, lwd=2, col = "blue")
  xlim <- range(S)
  ylim <- range(V[m,])
  segments(xlim[1], 0, xlim[2], 0, lty = 2)
  segments(0, ylim[1], 0, ylim[2], lty = 2)

  legend("topright",
         legend = c(expression(V[n]^0), expression(V(S,0))),
         col = c("black", "blue"),
         pch = c(15, NA),
         lwd = c(NA, 2),
         cex = 0.9,
         inset = c(0.02, 0.02),
         box.lwd = 2,
         bty = "n")

  #lines(BS,V[14,], col = "red")
  #lines(BS,V[17,], col = "purple")
  par(mar = c(5, 4, 4, 2) + 0.1)

}

#' Error plots for explicit FDM scheme at time 0
#'
#' @param N Number of spatial grid points
#' @param M Number of time steps
#' @param error If TRUE, produces error plot
#'
#' @return A plot of the difference between the explicit FDM
#' solution and the Black-Scholes price, against the asset price at
#' time 0, for any N and M
#' @export
plotsExplicitFDM <- function(N,M,error=TRUE) {

  sigma <- 0.4
  r <- 0.05
  T <- 1
  K <- 0.25
  Smax <- 1
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

  S <- seq(from = 0, to=1, by = delta_S)
  V0 <- V[1,]

  if (error){

    par(mfrow = c(1, 1))
    par(mar = c(4, 4.5, 4.5, 1) + 0.1)

    BSPut <- sapply(S, BlackScholesPutPrice, T=T, K=K, r=r, sigma=sigma, t=0, q=0)
    PutPriceErrors <- V0 - BSPut
    plot(S, PutPriceErrors, type = "p", pch = 15,
         xlab = expression(S[n]), ylab = expression(V[n]^0 - V(S[n],0)), cex = 0.9,
         main = paste("N =", N, ", M =", M))

    par(mar = c(5, 4, 4, 2) + 0.1)

  } else {
    plot(S, V0, type = "o", col = "black", pch = 15, xlab = "S", ylab = "V(S,0)")
    lines(S, V0, col = "blue")
  }

}


#' Error plots for implicit FDM scheme at time 0
#'
#' @param N Number of spatial grid points
#' @param M Number of time steps
#' @param error If TRUE, produces error plot
#'
#' @return A plot of the difference between the implicit FDM
#' solution and the Black-Scholes price, against the asset price at
#' time 0, for any N and M
#' @export
plotsImplicitFDM <- function(N,M,error=TRUE) {

  sigma <- 0.4
  r <- 0.05
  T <- 1
  K <- 0.25
  Smax <- 1
  delta_t <- T / M
  delta_S <- Smax / N

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

  S <- seq(from = 0, to=1, by = delta_S)
  V0 <- V[1,]

  if (error){

    par(mfrow = c(1, 1))
    par(mar = c(4, 4.5, 4.5, 1) + 0.1)

    BSPut <- sapply(S, BlackScholesPutPrice, T=T, K=K, r=r, sigma=sigma, t=0, q=0)
    PutPriceErrors <- V0 - BSPut
    plot(S, PutPriceErrors, type = "p", pch = 15,
         xlab = expression(S[n]), ylab = expression(V[n]^0 - V(S[n],0)), cex = 0.9,
         main = paste("N =", N, ", M =", M))

    par(mar = c(5, 4, 4, 2) + 0.1)

  } else {
    plot(S, V0, type = "o", col = "black", pch = 15, xlab = "S", ylab = "V(S,0)")
    lines(S, V0, col = "blue")
  }

}



