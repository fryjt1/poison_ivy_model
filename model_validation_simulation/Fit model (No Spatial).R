B = 10000

w.fit = array(w.vec,c(p-1,1,n))
Sigma.rep.inv = solve(bdiagRep(Sigma,n))
X.big.mat = apply(X.array,2,function(c)c)
beta.store = matrix(0,B,k); beta.store[1,] = beta
Sigma = TT

for(b in 2:B)
{
  beta.curr = matrix(beta.store[b-1,],ncol=1)
  
  for(j in 1:(p-1))
  {
    S.j.j = Sigma[j,j]
    S.nj.nj = Sigma[-j,-j]
    S.j.nj = matrix(Sigma[j,-j],nrow=1)
    V = S.j.j - S.j.nj%*%solve(S.nj.nj)%*%t(S.j.nj)

    for(i in 1:n)
    {
      E = X.nostack[i,]%*%beta.curr +  S.j.nj%*%solve(S.nj.nj)%*%(w.fit[-j,,i]-X.nostack[c(2,3),]%*%beta.curr)

      # if(y[i]==0){w.fit[j,,i] = rtruncnorm(1,a=-Inf,b=0,mean=E,sd=sqrt(V))}
      # if(y[i]==j){w.fit[j,,i] = rtruncnorm(1,a=max(0,w.fit[-j,,i]),b=Inf,mean=E,sd=sqrt(V))}
      # if((y[i]!=j)&&(y[i]!=0)){w.fit[j,,i] = rtruncnorm(1,a=-Inf,b=w.fit[y[i],,i],mean=E,sd=sqrt(V))}
      if(y[i]==j){w.fit[j,,i] = rtruncnorm(1,a=max(0,w.fit[-j,,i]),b=Inf,mean=E,sd=sqrt(V))}
      if(y[i]!=j){w.fit[j,,i] = rtruncnorm(1,a=-Inf,b=max(0,w.fit[-j,,i]),mean=E,sd=sqrt(V))}
    }
  }
  
  V = solve( t(X.big.mat)%*%Sigma.rep.inv%*%X.big.mat + diag(1/10000,k))
  E = V%*%t(X.big.mat)%*%Sigma.rep.inv%*%apply(w.fit,2,function(c)c)
  beta.store[b,] = rmvnorm(1,mean=E,sigma=V)
  
  if(b%%100==0){print(b)}
}

par(mfrow=c(3,3))
for(i in 1:length(beta))
{
  plot(beta.store[,i],type="l")
  abline(h=beta[i],col="chartreuse")
}