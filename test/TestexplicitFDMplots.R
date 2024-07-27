
# Testing
# Initialize parameters
sigma <- 0.4
r <- 0.05
T <- 1
K <- 0.25
Smax <- 1
N <- 16
M <- 32

delta_S <- Smax/N
S <- seq(from = 0, to=Smax, by = delta_S)

V <- explicitFDMEuropeanPut(sigma, r, T, K, Smax, N, M)
BSPutPrice <- sapply(S, BlackScholesPutPrice, T=T, K=K, r=r, sigma=sigma, t=0, q=0)
plotFDMV(S, V, BSPutPrice)
