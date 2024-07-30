

#' Black-Scholes price of European put option
#'
#' @param t Current time (default t=0)
#' @param T Maturity date
#' @param S Current asset price
#' @param K Strike price
#' @param r Risk-free interest rate
#' @param q Dividend yield
#' @param sigma Volatility test
#'
#' @return Price of European put option
#' @export
BlackScholesPutPrice <- function(t=0,T,S,K,r,q=0,sigma) {

  if(T <= t) {
    stop("T (maturity) must be greater than t (current time)")
  }

  dplus <- (log(S/K) + (r-q+0.5*sigma*sigma)*(T-t))/(sigma*sqrt(T-t))
  dminus <- dplus - sigma*base::sqrt(T-t)
  price <- K*exp(-r*(T-t))*stats::pnorm(-dminus) - S*base::exp(-q*(T-t))*stats::pnorm(-dplus)
  return(price)
}


#' Black-Scholes price of European call option
#'
#' @param t Current time (default t=0)
#' @param T Maturity date
#' @param S Current asset price
#' @param K Strike price
#' @param r Risk-free interest rate
#' @param q Dividend yield
#' @param sigma Volatility
#'
#' @return Price of European call option
#' @export
BlackScholesCallPrice <- function(t=0,T,S,K,r,q=0,sigma) {

  if(T <= t) {
    stop("T (maturity) must be greater than t (current time)")
  }

  dplus <- (log(S/K) + (r-q+0.5*sigma*sigma)*(T-t))/(sigma*sqrt(T-t))
  dminus <- dplus - sigma*base::sqrt(T-t)
  price <- S*base::exp(-q*(T-t))*stats::pnorm(dplus) - K*base::exp(-r*(T-t))*stats::pnorm(dminus)
  return(price)

}

# TODO: Sensitivities



























