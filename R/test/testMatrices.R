
###############################
### Tridiagonal matrix ########
###############################

a <- c(1, 2, 3)
b <- c(4, 5, 6, 7)
c <- c(8, 9, 10)
rParisianOptions::tridiagonalMatrix(a, b, c)

###############################
### Thomas algorithm ##########
###############################

A <- matrix(c(
  4, 1, 0, 0,
  1, 4, 1, 0,
  0, 1, 4, 1,
  0, 0, 1, 3
), nrow = 4, byrow = TRUE)

b <- c(15, 15, 15, 10)

xThomas <- rParisianOptions::thomasAlgorithmSolver(A, b)
x <- solve(A) %*% b
sum(abs(xThomas - x))

# TODO: Check robustness of tridiagonalMatrix()
