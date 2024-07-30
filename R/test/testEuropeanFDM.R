
###############################
### Initialize parameters #####
###############################

sigma <- 0.4
r <- 0.05
T <- 1
K <- 0.25
Smax <- 1
N <- 16
M <- 32

delta_S <- Smax/N
S <- seq(from = 0, to=Smax, by = delta_S)

###############################
### Explicit scheme ###########
###############################

# Plot of explicit FDM solutions vs Black-Scholes prices
V <- rParisianOptions::explicitFDMEuropeanPut(sigma, r, T, K, Smax, N, M)
BSPutPrice <- sapply(S, rParisianOptions::BlackScholesPutPrice, T=T, K=K, r=r, sigma=sigma, t=0, q=0)
rParisianOptions::plotFDMV(S, V, BSPutPrice)

# Error plots of explicit FDM
rParisianOptions::plotsExplicitFDM(16,16, error=TRUE)
rParisianOptions::plotsExplicitFDM(16,32, error=TRUE)
rParisianOptions::plotsExplicitFDM(32,32, error=TRUE)
rParisianOptions::plotsExplicitFDM(512,65536, error=TRUE)

#rParisianOptions::plotsExplicitFDM(32,48, error=TRUE)
#rParisianOptions::plotsExplicitFDM(48,48, error=TRUE)

# TODO: Update to include call/put in explicitFDMEuropeanPut()

###############################
### Implicit scheme ###########
###############################

# Plot of implicit FDM solutions vs Black-Scholes prices
V <- rParisianOptions::implicitFDMEuropeanPut(sigma, r, T, K, Smax, N, M)
rParisianOptions::plotFDMV(S, V, BSPutPrice)

# Error plots of implicit FDM
rParisianOptions::plotsImplicitFDM(16,16, error=TRUE)
rParisianOptions::plotsImplicitFDM(16,32, error=TRUE)
rParisianOptions::plotsImplicitFDM(32,32, error=TRUE)
rParisianOptions::plotsImplicitFDM(512,65536, error=TRUE)

#rParisianOptions::plotsImplicitFDM(32,256, error = TRUE)
#rParisianOptions::plotsImplicitFDM(256,16384, error = TRUE)




###############################
### Crank-Nicolson scheme #####
###############################











