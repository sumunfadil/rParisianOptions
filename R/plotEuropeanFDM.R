
# Plot of grid points -----------------------------------------------------


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

# Table of errors for different N and M -----------------------------------


#' Table of errors for different N and M
#'
#' @param N Sequence of spatial grid points
#' @param M Sequence of time steps
#' @param schem FDM scheme to be applied
#' @param view View dataframe
#'
#' @return Outputs table of error for different
#' N and M in latex
#' @export
errorTable <- function(N,M,scheme=explicitFDMEuropeanPutError,view=TRUE, theta=NULL) {

  schemeWrapper <- function(N, M) {

    if (identical(scheme, thetaFDMEuropeanPutError) && !is.null(theta)) {
      return(scheme(N, M, theta))
    } else {
      return(scheme(N, M))
    }
  }

  tableExplicit <- outer(N, M, FUN = Vectorize(schemeWrapper))
  df <- as.data.frame(t(tableExplicit))
  rownames(df) <- M
  colnames(df) <- N
  df <- format(df, scientific = TRUE, digits = 5)

  if (view) {

    View(df)

  } else {

    #library(xtable)
    latex_table <- xtable(df)
    print(latex_table, type = "latex", include.rownames = TRUE,
          include.colnames = TRUE)

    #file_path <- "dataframe_latex_table.tex"
    #capture.output(print(latex_table, type = "latex", include.rownames = TRUE,
    #                     include.colnames = TRUE), file = file_path)

  }

}


# Plotting function for FDM solution --------------------------------------


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



# Error plots for FDM scheme at time 0 ------------------------------------

#' Error plots for FDM scheme at time 0
#'
#' @param N Number of spatial grid points
#' @param M Number of time steps
#' @param error If TRUE, produces error plot
#'
#' @return A plot of the difference between the a FDM solution
#' and the Black-Scholes price, against the asset price at
#' time 0, for any N and M
#' @export
plotFDMError <- function(N,M,scheme=explicitFDMEuropeanPut,error=TRUE,return=FALSE,plotV=FALSE, theta=NULL) {

  sigma <- 0.4
  r <- 0.05
  T <- 1
  K <- 0.25
  Smax <- 1

  delta_S <- Smax / N

  if (identical(scheme, thetaFDMEuropeanPut) && !is.null(theta)) {
    V <- scheme(sigma, r, T, K, Smax, N, M, theta)
  } else {
    V <- scheme(sigma, r, T, K, Smax, N, M)
  }

  S <- seq(from = 0, to=1, by = delta_S)
  V0 <- V[1,]
  BSPut <- sapply(S, BlackScholesPutPrice, T=T, K=K, r=r, sigma=sigma, t=0, q=0)
  PutPriceErrors <- V0 - BSPut

  if (error){

    par(mfrow = c(1, 1))
    par(mar = c(4, 4.5, 4.5, 1) + 0.1)

    plot(S, PutPriceErrors, type = "p", pch = 15,
         xlab = expression(S[n]), ylab = expression(V[n]^0 - V(S[n],0)), cex = 0.9,
         main = paste("N =", N, ", M =", M))

    par(mar = c(5, 4, 4, 2) + 0.1)

  }

  if (plotV) {
    plot(S, V0, type = "o", col = "black", pch = 15, xlab = "S", ylab = "V(S,0)")
    lines(S, V0, col = "blue")
  }

  if (return) {
    return(list(S = S,Error = PutPriceErrors))
  }

}


#' Error plots for implicit FDM scheme at time 0 for
#' different N and M
#'
#' @param P1 Error plot 1 list
#' @param P2 Error plot 2 list
#' @param P2 Error plot 3 list
#' @param P4 Error plot 4 list
#'
#' @return A plot of the difference between the any FDM
#' solution and the Black-Scholes price, against the asset price at
#' time 0, for 4 pairs of (N,M)
#' @export
plotErrorMultiple <- function(P1,P2,P3,P4) {

  par(mfrow = c(1, 1))
  par(mar = c(4, 4.5, 4.5, 1) + 0.1)

  plot(P1$S, P1$Error, type = "p", pch = 0,
       main = "Error plot for FDM scheme",
       xlab = expression(S[n]), ylab = expression(V[n]^0 - V(S[n],0)))
  points(P2$S, P2$Error, type = "p", pch = 3, col = "blue")
  points(P3$S, P3$Error, type = "p", pch = 9, col = "purple")
  points(P4$S, P4$Error, type = "l", col = "red", lwd = 2)

  legend(x = 0.55, y = -0.00175,
         legend = c("N=16, M=16", "N=16, M=32", "N=32, M=32", "N=512, M=65536"),
         col = c("black", "blue", "purple", "red"),
         pch = c(0, 3, 9, NA),
         lty = c(NA, NA, NA, 1),
         lwd = c(NA, NA, NA, 2),
         cex = 0.9,
         bty = "n",
         y.intersp = 0.5
  )

  par(mar = c(5, 4, 4, 2) + 0.1)

}



