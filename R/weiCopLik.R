#' Weighted Gumbel copula log-likelihood
#'
#' This function evaluates the Gumbel copula log-likelihood
#' function computed with weighted observations.
#' @param x scalar > 1: values of the Gumbel copula parameter. 
#' @param y data matrix (nxd): observed data.
#' @param p_weights numerical vector (nx1) with elements in (0,1): weights
#' of the observations (in the EM algorithm, posterior probabilities).
#' @param gdp_pars matrix (2xd) with non-negative elements: scale (1st
#' row) and shape (2nd row) parameters of the marginal GPDs.
#' @return llik real: numerical value of the weighted log-likelihood
#' function
#' @export
#' @examples
#' ysim <- rMultLognGPD(100,.9,2,c(0,0),diag(2),c(.25,.5),c(1,1.5),2)
#' x0 <- c(.7,.2,1.3,.8,1.7)
#' res <- EMlogngpdmix(x0, y, 1000)
#' llik <- weiGpdLik(c(res$beta,res$xi),y,res$post)
#'
#' @importFrom Rdpack reprompt

weiCopLik <- function(x,y,p_weights,gpdpars)
{
  # x = theta 
  
  d <- ncol(y)
  n <- nrow(y)
  u <- matrix(0,n,d)
  f2gpd <- matrix(0,n,d)
  for (i in 1:d)
  {
    u[,i] <- evd::pgpd(y[,i],0,gpdpars[1,i],gpdpars[2,i])
    f2gpd[,i] <- evd::dgpd(y[,i],0,gpdpars[1,i],gpdpars[2,i])
  }
  u <- pmin(pmax(u,1e-10),1-1e-10)
  type.cop <- copula::gumbelCopula(x,dim=d)
  f2 <- copula::dCopula(u,type.cop) * apply(f2gpd,1,prod)
  f2 <- pmax(f2, .Machine$double.xmin)
  llik = sum(p_weights * log(f2))
  return(llik)
}
