bigT = 1000
contain = rep(NA,bigT)

for(t in 1:bigT)
{
  n = 100
  p = 3     # Number of options
  k = 1    # Covariates 
  my.grid = matrix(runif(2*n,0,1),ncol=2)
  dist.mat = as.matrix(dist(my.grid))
  lambda = .1
  H = exp(-dist.mat/(lambda)) + diag(1e-8,n)
  TT = riwish(v=p+2, S=diag(1,p-1))
  TT = TT/TT[1,1]
  TT = diag(1,p-1)
  
  beta = matrix(rnorm(k),ncol=1)
  X.nostack = matrix(rnorm(n*k),n,k); X = matrix(X.nostack[sort(rep(1:n,p-1)),],ncol=k)
  w.vec = matrix(rmvnorm(1,mean=X%*%beta,sigma=H%x%TT),ncol=1)
  locs = sort(rep(1:n,p-1)) # This tells me how to group together the w's
  y = apply(array(w.vec,c(p-1,1,n)),3,which.max)*apply(array(w.vec,c(p-1,1,n)),3,function(c)max(c)>0)
  
  I = 1000
  
  # Priors
  kappa = (p-1) + 2                                    # Priors on TT
  C = diag(1,p-1)
  b0 = matrix(rep(0,k),ncol=1)                         # Priors on Beta
  B.inv = diag(1/10,k)     
  
  # Initialize
  w.fit.vec = w.vec                                     #Initialize at truth for now
  beta.store = matrix(NA,I,k); beta.store[1,] = beta
  TT.store = matrix(NA,I,3); TT.store[1,] = c(1,1,0)    # (Variance, Variance, Covariance)
  lambda.store = rep(NA,I); lambda.store[1] = .1
  step = .1 # For Lambda proposal
  lambda.UB = 1
  
  beta.curr = matrix(beta.store[1,],ncol=1)
  TT.curr = diag(1,p-1)
  lambda.curr = lambda.store[1]
  H.curr = exp(-dist.mat/lambda.curr) + diag(1e-8,n)
  H.curr.inv = solve(H.curr)
  
  for(i in 2:I)
  {
    S.inv = H.curr.inv%x%solve(TT.curr)
    
    # Update Latent W's
    for(j in 1:(n*(p-1)))
    {
      loc.ind = which(locs==locs[j])                   # Tells me which elements in w.fit.vec are for this location
      curr.cat = which(w.fit.vec[loc.ind,]==w.fit.vec[j])
      w.other = matrix(w.fit.vec[setdiff(loc.ind,j),],ncol=1)
      
      w.star = w.fit.vec[-j,]
      X.k = matrix(X[-j,],ncol=k)
      x.j = matrix(X[j,],k,1)
      tau2.w = 1/S.inv[j,j]
      mu.w = t(x.j)%*%beta.curr + (S.inv[j,-j]*(-tau2.w))%*%(w.star-X.k%*%beta.curr)
      
      if(y[locs[j]]==curr.cat){w.fit.vec[j,] = rtruncnorm(1,a=max(0,w.other),b=Inf,mean=mu.w,sd=sqrt(tau2.w))}
      if(y[locs[j]]!=curr.cat){w.fit.vec[j,] = rtruncnorm(1,a=-Inf,b=max(0,w.other),mean=mu.w,sd=sqrt(tau2.w))}
    }
    
    # Update Beta
    V = solve(t(X)%*%S.inv%*%X + B.inv)
    E = V%*%(t(X)%*%S.inv%*%w.fit.vec + B.inv%*%b0)
    beta.curr = matrix(rmvnorm(1,mean=E,sigma=V),ncol=1)
    beta.store[i,] = beta.curr/sqrt(TT.curr[1,1])
    
    # Update TT
    w.mat = matrix(w.fit.vec,n,p-1,byrow=TRUE)
    mu.mat = matrix(X%*%beta.curr,n,p-1,byrow=TRUE)
    C.star = t(w.mat-mu.mat)%*%H.curr.inv%*%(w.mat-mu.mat) + C
    kappa.star = kappa + n
    TT.curr = riwish(v=kappa.star,S=C.star)
    TT.store[i,] = c(TT.curr[1,1],TT.curr[2,2],TT.curr[1,2])
    
    # Update lengthscale
    lambda.curr = lambda.store[i-1]
    lambda.prop = rtruncnorm(1,mean=lambda.curr,sd=step,a=0,b=lambda.UB)
    H.prop = exp((-dist.mat)/lambda.prop) + diag(1e-8,n)
    H.prop.inv = solve(H.prop)
    
    numer = ( -.5*(p-1)*determinant(H.prop,logarithm=TRUE)$modulus-.5*t(w.fit.vec-X%*%beta.curr)%*%(H.prop.inv%x%solve(TT.curr))%*%(w.fit.vec-X%*%beta.curr)
              - log(pnorm(lambda.UB,mean=lambda.prop,sd=step)-pnorm(0,mean=lambda.prop,sd=step)) )
    denom = ( -.5*(p-1)*determinant(H.curr,logarithm=TRUE)$modulus-.5*t(w.fit.vec-X%*%beta.curr)%*%(H.curr.inv%x%solve(TT.curr))%*%(w.fit.vec-X%*%beta.curr)
              - log(pnorm(lambda.UB,mean=lambda.curr,sd=step)-pnorm(0,mean=lambda.curr,sd=step))  )
    
    if(log(runif(1))<numer-denom)
    {
      lambda.store[i] = lambda.prop
      H.curr = H.prop
      H.curr.inv = H.prop.inv
    }else
    {
      lambda.store[i] = lambda.curr
    }

  }
    LB = quantile(beta.store,c(.025))
    UB = quantile(beta.store,c(.975))
    contain[t] = (beta>LB)&(beta<UB)
    print(t)
    print(mean(contain[1:t]))
}










