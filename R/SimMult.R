SimMult <- function(n,p,d,mu,Psi,xivec,betavec,gammap)
{
  # simulate n d-variate observations from a d-variate lognormal-GPD mixture,
  # where the GPD has a Gumbel dependence structure
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
