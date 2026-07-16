#' Bootstrap standard errors for the MLEs of a lognormal-GPD mixture
#'
#' This function draws a bootstrap sample and uses it to estimate the parameters of a lognormal-Pareto mixture distribution. Since this is typically called by LPfitEM, see the help of LPfitEM for examples.
#' @param x list: sequence of integers 1,...,K, where K is the mumber of datasets. Set x = 1 in case
#' of a single dataset.
#' @param x0 numerical vector (5x1): initial values of the parameters p,
#' \eqn{\mu}, \eqn{\sigma}, \eqn{\xi}, \eqn{\beta}. 
#' @param y numerical vector: observed sample.
#' @param maxiter non-negative integer: maximum number of iterations of the EM algorithm.
#' @return Estimated parameters obtained from a bootstrap sample.
#' @details At each bootstrap replication, the mixture is estimated via the EM algorithm.
#' @export

EMMultBoot = function(x,y,eps,maxiter)
{
  samSiz <- nrow(y)
  indici = sample(samSiz, samSiz, replace = TRUE)
  yboot = y[indici,]
  temp <- EMMultp(yboot,eps,maxiter)
  res <- temp$pars
  results <- list(res=res)
}
