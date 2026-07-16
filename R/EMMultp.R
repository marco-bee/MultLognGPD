#' Multivariate estimation
#'
#' This function estimates a multivariate lognormal - generalized Pareto mixture
#' by means of the EM algorithm. Optionally, bootstrap standard errors are
#' computed via parallel computing.
#' @param ymix data matrix (nxd): observed data.
#' @param eps positive real: tolerance in the stopping criterion of the EM algorithm.
#' @param maxiter positive integer: maximum number of iterations of the EM algorithm.
#' @param nboot positive integer: number of bootstrap replications for the
#' computation of the standard errors (defaults to 0).
#' @return A list with the following elements is returned:
#' "pars" = estimated values of the parameters,
#' "loglik" = maximimzed log-likelihood,
#' "niter" = number of iterations,
#' "post_p" = posterior probabilities of all observations,
#' bootEst = matrix of parameter estimates at each bootstrap replications (only if nboot > 0).
#' bootStd = bootstrap standard errors of each parameter (only if nboot > 0).
#' @export
#' @examples
#' y <- rMultLognGPD(100,.9,2,c(0,0),diag(2),c(.25,.5),c(1,1.5),2)
#' x0 <- c(.7,.2,1.3,.8,1.7)
#' res <- EMlogngpdmix(x0, y, 1000)
#'
#' @importFrom Rdpack reprompt

EMMultp <- function(ymix,eps,maxiter,nboot = 0)
{
  nit = 1
  n = nrow(ymix)
  d = ncol(ymix)
  mu <- colMeans(log(ymix))
  Psi <- var(log(ymix))
  Psi.low <- Psi[lower.tri(Psi,diag = TRUE)]
  gpdparst <- matrix(0,2,d)
  loglik <- rep(0,maxiter)
  post_p <- rep(0,n)		# open matrix for posterior probabilities
  u <- matrix(0,n,d)		
  f2gpd <- matrix(0,n,d)		
  f1 <- rep(0,n)
  dev1 <- matrix(0,n,d)
  change = 100
  
  # initialization
  
  p <-  .5
  for (i in 1:d)
  {
    fitted_GPD <- try(evd::fpot(ymix[,i], threshold = quantile(ymix[,i],.7)))
    if (inherits(fitted_GPD, "try-error")) {
      gpdparst[,i] <- c(1,.5)
    }
    else
    {
      gpdparst[,i] <- as.vector(fitted_GPD$estimate) # beta, xi
    }
  }
  uu <- copula::pobs(ymix)                # pseudo-observations
  fit.tau <- copula::fitCopula(copula::gumbelCopula(dim=d), uu, method="itau")
  gammap <- as.double(coefficients(fit.tau))
  
  parold <- c(p,mu,Psi.low,matrixcalc::vec(gpdparst),gammap)
  
  while ((change > eps || change < 0) && nit <= maxiter)
  {
    print(nit)
    f1 <- compositions::dlnorm.rplus(ymix,mu,Psi)   
    gum.cop <- copula::gumbelCopula(gammap,dim=d)
    for (i in 1:d)
    {
      u[,i] <- evd::pgpd(ymix[,i],0,gpdparst[1,i],gpdparst[2,i])
      f2gpd[,i] <- evd::dgpd(ymix[,i],0,gpdparst[1,i],gpdparst[2,i])
    }
    u <- pmin(pmax(u,1e-10),1-1e-10)
    f2 <- copula::dCopula(u,gum.cop) * apply(f2gpd,1,prod)
    f <- p * f1 + (1-p) * f2       
    post_p <- p * f1 / f           # E-step
    
    p <- mean(post_p)              # M step: prior probabilities
    for (i in 1:d)
    {
      mu[i] <- post_p %*% log(ymix[,i]) / (n * p) # M step: lognormal parameters
    }
    mu <- as.double(mu)
    for (i in 1:d)
    {
      dev1[,i]=(log(ymix[,i])-mu[i]) * sqrt(post_p)
    }
    Psi <- t(dev1) %*% dev1 /(n*p)
    Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
    Psi.low <- Psi[lower.tri(Psi,diag=TRUE)]
    
    # CM-step 1: marginal GPD parameters 
    
    gpdpars = optim(log(matrixcalc::vec(gpdparst)),weiGpdLik,gr=NULL,
                    gum.cop,ymix,1-post_p,control=list(fnscale=-1))
    gpdparst <- matrix(exp(gpdpars$par),2,d)
    
    # CM-step 2: copula parameter(s) 
    
    coppars = optimize(weiCopLik,c(1,10),
                       ymix,1-post_p,gpdparst,maximum=TRUE)
    gammap <- coppars$maximum
    
    pars <- c(p,mu,Psi.low,matrixcalc::vec(gpdparst),gammap)
    loglik[nit] <- sum(log(dMultLognGPD(ymix,p,mu,Psi,gpdparst,gammap)))        # evaluate log-likelihood function
    
    if (nboot == 0)
    {
      if (nit>2)
      {
        alpha <- (loglik[nit] - loglik[nit-1]) / (loglik[nit-1] - loglik[nit-2]) # Aitken's acceleration (Cui et al. 2026, p. 4)
        ell_inf <- loglik[nit] + (1/(1-alpha)) * (loglik[nit] - loglik[nit-1])
        parold = pars
        change <- ell_inf - loglik[nit]
        if ((change < eps && change > 0) || nit==maxiter)
        {
          out <- list(pars = pars, loglik = loglik, niter = nit, post_p = post_p)
          return(out)
        }
      }
      nit <- nit + 1
    }
    else
    {
      nreps.list <- sapply(1:nboot, list)
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        n.cores <- 2L
      } else {
        n.cores <- parallel::detectCores()
      }
      clust <- parallel::makeCluster(n.cores)
      BootMat = matrix(0,nboot,1+d+d*(d+1)/2+2*d+1)
      temp <- parallel::parLapply(clust,nreps.list, EMMultBoot,ymix,eps,maxiter)
      parallel::stopCluster(cl=clust)
      for (i in 1:nboot)
      {
        BootMat[i,] = as.vector(unlist(temp[[i]]))
      }
      stddev = apply(BootMat,2,sd,na.rm=TRUE)
      out <- list(pars = pars, loglik = loglik, niter = nit, post_p = post_p,bootEst=BootMat,bootStd=stddev)
      return(out)
      # if ((change < eps && change > 0) || nit==maxiter)
      # {
      #   out <- list(pars = pars, loglik = loglik, niter = nit, post_p = post_p)
      #   return(out)
      # }
      nit <-  nit + 1
    }
    # nit <- nit + 1
  }
}
