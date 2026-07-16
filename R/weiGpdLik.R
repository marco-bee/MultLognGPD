#' Weighted GPD log-likelihood
#'
#' This function evaluates the zero-mean generalized Pareto log-likelihood
#' function computed with weighted observations.
#' @param x numerical vector (2dx1): values of the parameters \eqn{\beta}
#' and \eqn{\xi} for each marginal: (\eqn{\beta_1},\eqn{\xi_1},\eqn{\beta_2},\eqn{\xi_2},...). 
#' @param y data matrix (nxd): observed data.
#' @param p_weights numerical vector (nx1) with elements in (0,1): weights
#' of the observations (in the EM algorithm, posterior probabilities).
#' @return llik real: numerical value of the weighted GPD log-likelihood
#' function
#' @export
#' @examples
#' ysim <- rMultLognGPD(100,.9,2,c(0,0),diag(2),c(.25,.5),c(1,1.5),2)
#' x0 <- c(.7,.2,1.3,.8,1.7)
#' res <- EMlogngpdmix(x0, y, 1000)
#' llik <- weiGpdLik(c(res$beta,res$xi),y,res$post)
#'
#' @importFrom Rdpack reprompt

weiGpdLik <- function(x,type.cop,y,p_weights)
{
  # x[i] = log(xi_i), x[i+1]= log(beta_i)
  
  n <- nrow(y)
  d <-  ncol(y)
  u <- matrix(0,n,d)
  f2mat <-  matrix(0,n,d)
  for (i in 1:d)
  {
    u[,i] <- evd::pgpd(y[,i],0,exp(x[-1+2*i]),exp(x[2*i]))
    f2mat[,i] <- evd::dgpd(y[,i],0,exp(x[-1+2*i]),exp(x[2*i]))
  }
  u <- pmin(pmax(u,1e-10),1-1e-10)
  f2 <- copula::dCopula(u,type.cop) * apply(f2mat,1,prod)
  f2 <- pmax(f2, .Machine$double.xmin)
  llik = sum(p_weights * log(f2))
  return(llik)
}
