library(truncnorm)
library(MCMCpack)
library(metaSEM)
library(abind)

# Model Specification 
# p = # of categories
# W_i = X_i B + e_i, e_i ~ N(0,Sigma)
# Sigma = matrix( 1, g', g, Phi + gg' ) by row
# B ~ N(b0,B.inv)
# Phi.inv ~ Wishart(kappa,C)
# C \set (kappa - p + 1)(Delta - B.inv) = (kappa - p + 1)(I - tauI) for prior specification

n = 100
p = 4     # Number of options
k = 5     # Covariates for each category

#  the covariance between Z_i1,...,Z_i(p-1)
kappa = p + 2
tau = 1/100
C = (kappa - p + 1)*(1 - tau)*diag(1,p-2) # This is weakly informative
g = matrix(rnorm(p-2,mean=0,sd=sqrt(tau)),ncol=1)
Phi = riwish(v=kappa, S=C)
Sigma = rbind(c(1,t(g)),cbind(g,Phi+g%*%t(g)))

# beta = matrix(rnorm(k*(p-1)),ncol=1)
beta = matrix(rnorm(k),ncol=1)
  
X.array = array(NA,c(p-1,k,n))
w.array = array(NA,c(p-1,1,n))
y = rep(NA,n)

for(i in 1:n)
{
  X.array[,,i] = matrix(1,nrow=p-1,ncol=1)%x%matrix(rnorm(k),nrow=1) # Same covariates for each category
  w.array[,,i] = X.array[,,i]%*%beta + matrix(rmvnorm(1,mean=rep(0,p-1),sigma=Sigma),ncol=1)
  y[i] = max(w.array[,,i]>0)*which.max(w.array[,,i]) + 0
}



