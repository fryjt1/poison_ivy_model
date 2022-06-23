Augment.X = function(X,p)
{
  n = nrow(X)
  X.new = diag(1,p-1)%x%X
  X.new = X.new[c(apply(matrix(1:n,ncol=1),1,function(x)seq(x,n*(p-1)+x-1,by=n))),]
  return(X.new)
}
Sample.W = function(y,locs,w.vec.curr,S.inv,XB)
{
  w.draw = w.vec.curr
  n = length(unique(locs))
  p = length(locs)/n + 1
  
  for(j in 1:(n*(p-1)))
  {
    loc.ind = which(locs==locs[j])  # Tells me which elements in w.vec.curr are for this location
    curr.cat = which(w.draw[loc.ind,]==w.draw[j])
    w.other = matrix(w.draw[setdiff(loc.ind,j),],ncol=1)
    
    w.star = w.draw[-j,]
    tau2.w = 1/S.inv[j,j]
    mu.w = XB[j,] + tcrossprod(matrix(S.inv[j,-j],nrow=1)*(-tau2.w),(w.star-XB[-j,]))
    
    if(y[locs[j]]==curr.cat){w.draw[j,] = rtruncnorm(1,a=max(0,w.other),b=Inf,mean=mu.w,sd=sqrt(tau2.w))}
    if(y[locs[j]]!=curr.cat){w.draw[j,] = rtruncnorm(1,a=-Inf,b=max(0,w.other),mean=mu.w,sd=sqrt(tau2.w))}
  }
  return(w.draw)
}

Sample.Beta = function(X,S.inv,w.vec,B.inv,b0)
{
  XtS.inv = t(X)%*%S.inv
  V = solve(XtS.inv%*%X + B.inv)
  V = symmpart(V)
  E = V%*%(XtS.inv%*%w.vec + B.inv%*%b0)
  draw = matrix(rmvnorm(1,mean=E,sigma=V),ncol=1)
  return(draw)
}

# Currently, TT is fixed!!!!!!!!!
Sample.TT = function(w.vec,XB,H.inv,C,kappa)
{
  n = ncol(H.inv)
  p = nrow(w.vec)/n + 1
  w.mat = matrix(w.vec,n,p-1,byrow=TRUE)
  mu.mat = matrix(XB,n,p-1,byrow=TRUE)
  matrix.SSE = t(w.mat-mu.mat)%*%H.inv%*%(w.mat-mu.mat)
  C.star = matrix.SSE + C
  kappa.star = kappa + n
  draw = riwish(v=kappa.star,S=C.star)
  draw = matrix(c(2,1,1,1,2,1,1,1,2),3,3,byrow=TRUE)
  return(list(draw=draw,matrix.SSE=matrix.SSE))
}

Sample.Lambda = function(lambda.curr,lambda.UB,step,H.curr,H.inv.curr,TT.curr,w.vec,XB,matrix.SSE)
{
  n = nrow(H.curr)
  p = nrow(TT.curr) + 1
  
  w.mat = matrix(w.vec,n,p-1,byrow=TRUE)
  mu.mat = matrix(XB,n,p-1,byrow=TRUE)
  
  lambda.prop = rtruncnorm(1,mean=lambda.curr,sd=step,a=0,b=lambda.UB)
  H.prop = exp((-dist.mat)/lambda.prop) + diag(1e-8,n)
  H.prop.inv = solve(H.prop)
  
  numer = ( -.5*(p-1)*determinant(H.prop,logarithm=TRUE)$modulus-.5*tr(solve(TT.curr)%*%t(w.mat-mu.mat)%*%H.prop.inv%*%(w.mat-mu.mat))
            - log(pnorm(lambda.UB,mean=lambda.prop,sd=step)-pnorm(0,mean=lambda.prop,sd=step)) )
  denom = ( -.5*(p-1)*determinant(H.curr,logarithm=TRUE)$modulus-.5*tr(solve(TT.curr)%*%matrix.SSE)
            - log(pnorm(lambda.UB,mean=lambda.curr,sd=step)-pnorm(0,mean=lambda.curr,sd=step))  )
  # matrix.SSE has H.curr.inv in it --> see TT update
  
  if(log(runif(1))<numer-denom)
  {
    lambda.draw = lambda.prop
    H.draw = H.prop
    H.inv.draw = H.prop.inv
  }else
  {
    lambda.draw = lambda.curr
    H.draw = H.curr
    H.inv.draw = H.inv.curr
  }
  
  return(list(lambda=lambda.draw,H=H.draw,H.inv=H.inv.draw))
}

Impute.Cont = function(H.inv,TT.inv,w.vec,X.nostack.miss,B,miss.ind,col.nums,m0,C0)
{
  n = nrow(H.inv)
  p = nrow(TT.inv) + 1
  B.mat = matrix(B,ncol=p-1)
  w.mat = matrix(w.vec,n,p-1,byrow=TRUE) - X.nostack.miss[,-col.nums]%*%B.mat[-col.nums,]
  x = matrix(X.nostack.miss[,col.nums],ncol=1)
  Bm = matrix(B.mat[col.nums,],ncol=1)
  x1 = matrix(x[-miss.ind,],ncol=1)
  
  V = solve(H.inv[miss.ind,miss.ind]*c(t(Bm)%*%TT.inv%*%Bm) + solve(C0))
  E = V%*%( (H.inv[miss.ind,-miss.ind]%*%(w.mat[-miss.ind,]-x1%*%t(Bm)) 
      + H.inv[miss.ind,miss.ind]%*%w.mat[miss.ind,])%*%TT.inv%*%Bm + solve(C0)%*%m0)
  
  draw = rmvnorm(1,mean=E,sigma=V)
  return(draw)
}

#H.inv = H.inv.curr; TT.inv = solve(TT); w.vec = w.fit.vec;X.nostack.miss = X.nostack.miss
#B = beta; miss.inds = miss.ind; col.nums = col.nums

Impute.Cat = function(H.inv,TT.inv,w.vec,X.nostack.miss,B,miss.inds,col.nums) 
# Needs to do one missing observation at a time
# col.nums tells me the columns of X.nostack that indicate categories
{
  n = nrow(H.inv)
  p = nrow(TT.inv) + 1
  X = X.nostack.miss
  B.mat = matrix(B,ncol=p-1)
  B.curl = B.mat[col.nums,]
  cat.ind = rep(NA,length(miss.inds))
  
  for(i in 1:length(miss.inds))
  {
    curr.ind = miss.inds[i]
    W.curl =  matrix(w.vec,n,p-1,byrow=TRUE) - X[,-col.nums]%*%B.mat[-col.nums,]
    X.curl = X[,col.nums]
    
    x2 = matrix(X.curl[curr.ind,],ncol=1)
    X1 = X.curl[-curr.ind,]
    z2 = matrix(W.curl[curr.ind,],ncol=1)
    W1 = W.curl[-curr.ind,]
    w2 = matrix(W.curl[curr.ind,],ncol=1)
    
    S.21X1 = H.inv[-curr.ind,curr.ind]%*%X1
    S.21W1 = H.inv[-curr.ind,curr.ind]%*%W1
    
    log.probs = rep(NA,length(col.nums)+1)  # one for each cat in X, one for all 0's
    x.options = rbind(rep(0,length(col.nums)),diag(1,length(col.nums)))
    T.invB = TT.inv%*%t(B.curl)
    BtT.invB = B.curl%*%T.invB
    
    for(j in 1:nrow(x.options))
    {
      x2 = matrix(x.options[j,],ncol=1)

      log.probs[j] = -.5*tr( x2%*%S.21X1%*%BtT.invB + t(x2%*%S.21X1)%*%BtT.invB + x2%*%H.inv[curr.ind,curr.ind]%*%t(x2)%*%BtT.invB
             -2*(x2%*%S.21W1%*%T.invB + x2%*%H.inv[curr.ind,curr.ind]%*%t(w2)%*%T.invB) 
              )
    }
    my.max = which.max(log.probs)
    log.denom = log.probs[my.max] + log(1+sum(exp(log.probs[-my.max]-log.probs[my.max])))
    probs = exp(log.probs-log.denom)
    cat.ind[i] = sample(1:nrow(x.options),size=1,prob=probs)
    x.draw = x.options[cat.ind[i],]
    X[curr.ind,col.nums] = x.draw
  }
  return(list(X=X,cat.ind=cat.ind))
}


# ###################### SCRATCH WORK
# 
# W = matrix(w.vec,n,p-1,byrow=TRUE)
# log.probs2 = rep(NA,length(col.nums)+1) 
# log.probs3 = log.probs2
# for(k in 1:nrow(x.options))
# {
#   X.temp = X.nostack.miss
#   X.temp[curr.ind,col.nums] = x.options[k,]
#   X.curl.temp = X.curl
#   X.curl.temp[curr.ind,] = x.options[k,]
#   
#   log.probs2[k] = -.5*tr( TT.inv%*%t(W - X.temp%*%B.mat)%*%H.inv%*%(W - X.temp%*%B.mat) )
#   log.probs3[k] = -.5*( tr(TT.inv%*%t(B.curl)%*%t(X.curl.temp)%*%H.inv%*%X.curl.temp%*%B.curl 
#                           - 2*(t(X.curl.temp)%*%H.inv%*%W.curl%*%TT.inv%*%t(B.curl)) )
#                           )
# }
# 
# my.max = which.max(log.probs)
# log.denom = log.probs[my.max] + log(1+sum(exp(log.probs[-my.max]-log.probs[my.max])))
# my.max = which.max(log.probs2)
# log.denom2 = log.probs2[my.max] + log(1+sum(exp(log.probs2[-my.max]-log.probs2[my.max])))
# my.max = which.max(log.probs2)
# log.denom3 = log.probs3[my.max] + log(1+sum(exp(log.probs3[-my.max]-log.probs3[my.max])))
# 
# probs = exp(log.probs-log.denom)
# probs2 = exp(log.probs2-log.denom2)
# probs3 = exp(log.probs3-log.denom3)
