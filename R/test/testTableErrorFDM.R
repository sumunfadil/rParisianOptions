

# Explicit FDM Error tables for different M and N -------------------------

M <- 2^(4:16)
N <- 2^(4:9)
rParisianOptions::errorTable(N,M)


# Implicit FDM Error tables for different M and N -------------------------

M <- 2^(4:9)
N <- 2^(4:9)
rParisianOptions::errorTable(N,M,scheme=rParisianOptions::implicitFDMEuropeanPutError)


# Crank-Nicolson FDM Error tables for different M and N -------------------

M <- 2^(4:9)
N <- 2^(4:9)
rParisianOptions::errorTable(N,M,scheme=rParisianOptions::thetaFDMEuropeanPutError,theta=0.5)










