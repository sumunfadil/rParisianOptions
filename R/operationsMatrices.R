

#' Making a tridiagonal matrix
#'
#' @param a Sub-diagonal
#' @param b Main diagonal
#' @param c Super-diagonal
#'
#' @return A tridiagonal matrix
#' @export
tridiagonalMatrix <- function(a, b, c) {

  N <- length(b)
  mat <- matrix(0, nrow = N, ncol = N)

  # Main diagonal
  diag(mat) <- b

  # Sub-diagonal
  if (N > 1) {
    mat[cbind(2:N, 1:(N-1))] <- a
  }

  # Super-diagonal
  if (N > 1) {
    mat[cbind(1:(N-1), 2:N)] <- c
  }

  return(mat)
}


#' Thomas algorithm
#'
#' @param A Square matrix
#' @param b A vector
#'
#' @return The solution x of the equation
#' Ax = b
#' @export
thomasAlgorithmSolver <- function(A, b) {

  N <- length(b)

  # Extract the diagonals
  a <- numeric(N - 1)
  b_diag <- numeric(N)
  c <- numeric(N - 1)

  for (i in 1:(N - 1)) {
    a[i] <- A[i + 1, i]
    c[i] <- A[i, i + 1]
  }
  for (i in 1:N) {
    b_diag[i] <- A[i, i]
  }

  # Forward elimination
  c_prime <- numeric(N - 1)
  d_prime <- numeric(N)

  c_prime[1] <- c[1] / b_diag[1]
  d_prime[1] <- b[1] / b_diag[1]

  for (i in 2:(N - 1)) {
    c_prime[i] <- c[i] / (b_diag[i] - a[i - 1] * c_prime[i - 1])
  }

  for (i in 2:N) {
    d_prime[i] <- (b[i] - a[i - 1] * d_prime[i - 1]) / (b_diag[i] - a[i - 1] * c_prime[i - 1])
  }

  # Backward substitution
  x <- numeric(N)
  x[N] <- d_prime[N]

  for (i in (N - 1):1) {
    x[i] <- d_prime[i] - c_prime[i] * x[i + 1]
  }

  return(x)
}





