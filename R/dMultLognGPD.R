#' Density of a multivariate lognormal-GPD mixture, where the GPD has a
#' Gumbel dependence structure
#'
#' This function evaluaates the density of a multivariate lognormal-GPD
#' mixture.
#' @param x (dx1) positive vector: points where the function is evaluated.
#' @param p real, 0<p<1: prior probability
#' @param mu (dx1) vector: log-mean of the multivariate lognormal distribution.
#' @param Psi positive definite symmetric (dxd) matrix: covariance matrix of
#' the underlying multivariate normal distribution.
#' @param gdppars (2xd)-vector of positive real numbers: matrix whose i-th 
#' column contains the scale and shape parameters of the i-th
#' marginal generalized Pareto distribution.
#' @param gammap real > 1: parameter of the Gumbel copula.
#' @return dmix n-vector: n values of the
#' d-variate lognormal - generalized Pareto mixture density.
#' @import stats
#' @export
#' @examples
#' y <- dMultLognGPD(c(1,2),.5,c(0,0),diag(2),rbind(c(1,1.5),c(.25,5)),2)
#'
#' @importFrom Rdpack reprompt

dMultLognGPD <- function(x,p,mu,Psi,gpdpars,gammap)
{
  n <- nrow(x)
  d <- ncol(x)
  u <- matrix(0,n,d)
  f2gpd <- matrix(0,n,d)
  d1 <- compositions::dlnorm.rplus(x,mu,Psi)
  gum.cop <- copula::gumbelCopula(gammap,dim=d)
  for (i in 1:d)
  {
    u[,i] <- evd::pgpd(x[,i],0,gpdpars[1,i],gpdpars[2,i])
    f2gpd[,i] <- evd::dgpd(x[,i],0,gpdpars[1,i],gpdpars[2,i])
  }
  u <- pmin(pmax(u,1e-10),1-1e-10)
  d2 <- copula::dCopula(u,gum.cop) * apply(f2gpd,1,prod)
  dmix <- p * d1 + (1-p) * d2
  return(dmix)
}
