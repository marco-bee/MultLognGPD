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
#' @param xivec d-vector of positive real numbers: shape parameters of the 
#' marginal generalized Pareto distributions.
#' @param betavec d-vector of positive real numbers: scale parameters of the 
#' marginal generalized Pareto distributions.
#' @return ysim (n x d) matrix: n random vectors from the
#' d-variate lognormal - generalized Pareto mixture.
#' @import stats
#' @export
#' @examples
#' ysim <- rMultLognGPD(100,.9,2,c(0,0),diag(2),c(.25,.5),c(1,1.5),2)
#'
#' @importFrom Rdpack reprompt


rMultLognGPD <- function(n,p,d,mu,Psi,xivec,betavec,gammap)
{
  gum.cop <- copula::gumbelCopula(gammap,dim=d)
  n1 <- sum(rbinom(n,1,p))
  n2 <- n - n1
  Y1 <- compositions::rlnorm.rplus(n1,mu,Psi)
  u <- copula::rCopula(n2,gum.cop)
  Y2 <-  matrix(0,n2,d)
  for (i in 1:d)
  {
    Y2[,i] <- evd::qgpd(u[,i],0,betavec[i],xivec[i])
  }
  Y <- rbind(Y1,Y2)
  return(Y)
}
