

# Initialize parameters ---------------------------------------------------


sigma <- 0.4
r <- 0.05
T <- 1
K <- 0.25
Smax <- 1
N <- 16
M <- 32

delta_S <- Smax/N
S <- seq(from = 0, to=Smax, by = delta_S)


# Explicit scheme ---------------------------------------------------------


# Plot of explicit FDM solutions vs Black-Scholes prices
V <- rParisianOptions::explicitFDMEuropeanPut(sigma, r, T, K, Smax, N, M)
BSPutPrice <- sapply(S, rParisianOptions::BlackScholesPutPrice, T=T, K=K, r=r, sigma=sigma, t=0, q=0)
rParisianOptions::plotFDMV(S, V, BSPutPrice)

# Error plots of explicit FDM
rParisianOptions::plotFDMError(16,16)
rParisianOptions::plotFDMError(16,32)
rParisianOptions::plotFDMError(32,32)
rParisianOptions::plotFDMError(512,65536)

#rParisianOptions::plotsExplicitFDM(32,48)
#rParisianOptions::plotsExplicitFDM(48,48)

# TODO: Update to include call/put in explicitFDMEuropeanPut()



# Implicit scheme ---------------------------------------------------------


# Plot of implicit FDM solutions vs Black-Scholes prices
V <- rParisianOptions::implicitFDMEuropeanPut(sigma, r, T, K, Smax, N, M)
rParisianOptions::plotFDMV(S, V, BSPutPrice)

# Error plots of implicit FDM
rParisianOptions::plotFDMError(16,16,scheme=rParisianOptions::implicitFDMEuropeanPut)
rParisianOptions::plotFDMError(16,32,scheme=rParisianOptions::implicitFDMEuropeanPut)
rParisianOptions::plotFDMError(32,32,scheme=rParisianOptions::implicitFDMEuropeanPut)
rParisianOptions::plotFDMError(512,65536,scheme=rParisianOptions::implicitFDMEuropeanPut)

#rParisianOptions::plotsImplicitFDM(32,256,scheme=rParisianOptions::implicitFDMEuropeanPut)
#rParisianOptions::plotsImplicitFDM(256,16384,scheme=rParisianOptions::implicitFDMEuropeanPut)


# Error plots for implicit FDM (4 plots)
P1 <- rParisianOptions::plotFDMError(16,16,scheme=rParisianOptions::implicitFDMEuropeanPut,error=FALSE,return=TRUE)
P2 <- rParisianOptions::plotFDMError(16,32,scheme=rParisianOptions::implicitFDMEuropeanPut,error=FALSE,return=TRUE)
P3 <- rParisianOptions::plotFDMError(32,32,scheme=rParisianOptions::implicitFDMEuropeanPut,error=FALSE,return=TRUE)
P4 <- rParisianOptions::plotFDMError(512,65536,scheme=rParisianOptions::implicitFDMEuropeanPut,error=FALSE,return=TRUE)
rParisianOptions::plotErrorMultiple(P1,P2,P3,P4)


# Crank-Nicolson scheme ---------------------------------------------------


# Plot of theta FDM solutions vs Black-Scholes prices
V <- rParisianOptions::thetaFDMEuropeanPut(sigma, r, T, K, Smax, N, M)
rParisianOptions::plotFDMV(S, V, BSPutPrice)

# Error plots of theta FDM
rParisianOptions::plotFDMError(16,16,scheme=rParisianOptions::thetaFDMEuropeanPut, theta=0.5)
rParisianOptions::plotFDMError(16,32,scheme=rParisianOptions::thetaFDMEuropeanPut, theta=0.5)
rParisianOptions::plotFDMError(32,32,scheme=rParisianOptions::thetaFDMEuropeanPut, theta=0.5)
rParisianOptions::plotFDMError(512,65536,scheme=rParisianOptions::thetaFDMEuropeanPut, theta=0.5)

#rParisianOptions::plotFDMError(32,256,scheme=rParisianOptions::thetaFDMEuropeanPut)
#rParisianOptions::plotFDMError(256,16384,scheme=rParisianOptions::thetaFDMEuropeanPut)


# Error plots for implicit FDM (4 plots)
P4 <- rParisianOptions::plotFDMError(16,16,scheme=rParisianOptions::thetaFDMEuropeanPut,error=FALSE,return=TRUE, theta=0.5)
P3 <- rParisianOptions::plotFDMError(16,32,scheme=rParisianOptions::thetaFDMEuropeanPut,error=FALSE,return=TRUE, theta=0.5)
P2 <- rParisianOptions::plotFDMError(32,32,scheme=rParisianOptions::thetaFDMEuropeanPut,error=FALSE,return=TRUE, theta=0.5)
P1 <- rParisianOptions::plotFDMError(512,65536,scheme=rParisianOptions::thetaFDMEuropeanPut,error=FALSE,return=TRUE, theta=0.5)
rParisianOptions::plotErrorMultiple(P1,P2,P3,P4)











