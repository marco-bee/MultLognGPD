#' Simulation of a multivariate lognormal-GPD mixture, where the GPD has a
#' Gumbel dependence structure
#'
#' This function simulates a multivariate lognormal-GPD mixture.
#' @param n positive integer: number of observations sampled.
#' @param p real, 0<p<1: prior probability
#' @param d positive integer: dimension of the distribution.
#' @param mu (dx1) vector: log-mean of the multivariate lognormal distribution.
#' @param Psi positive definite symmetric (dxd) matrix: covariance matrix of
#' the underlying multivariate normal distribution.
#' @param gdppars (2xd)-vector of positive real numbers: matrix whose i-th 
#' column contains the scale and shape parameters of the i-th
#' marginal generalized Pareto distribution.
#' @return ysim (n x d) matrix: n random vectors from the
#' d-variate lognormal - generalized Pareto mixture.
#' @import stats
#' @export
#' @examples
#' gpdpars <- cbind(c(1,0.5),c(1.5,0.25))
#' ysim <- rMultLognGPD(100,.9,2,c(0,0),diag(2),gpdpars,2)
#'
#' @importFrom Rdpack reprompt


rMultLognGPD <- function(n,p,d,mu,Psi,gpdpars,gammap)
{
  gum.cop <- copula::gumbelCopula(gammap,dim=d)
  n1 <- sum(rbinom(n,1,p))
  n2 <- n - n1
  Y1 <- compositions::rlnorm.rplus(n1,mu,Psi)
  u <- copula::rCopula(n2,gum.cop)
  Y2 <-  matrix(0,n2,d)
  for (i in 1:d)
  {
    Y2[,i] <- evd::qgpd(u[,i],0,gpdpars[1,i],gpdpars[2,i])
  }
  Y <- rbind(Y1,Y2)
  return(Y)
}
