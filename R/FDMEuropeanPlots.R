
# Plotting function for solution of FDM compared with BS


#' Plotting function for FDM solution
#'
#' @param S Asset price
#' @param V Matrix of FDM values
#' @param BS Black-Scholes prices
#' @param m Time-step
#'
#' @return NA
#' @export
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
