library(truncnorm)
library(MCMCpack)
library(metaSEM)
library(abind)

n = 2000
p = 4     # Number of options
k = 50    # Covariates 
my.grid = matrix(runif(2*n,0,1),ncol=2)
dist.mat = as.matrix(dist(my.grid))
lambda = .1
H = exp(-dist.mat/(lambda)) + diag(1e-8,n)
TT = riwish(v=p+2, S=diag(1,p-1))
TT = TT/TT[1,1]
TT = diag(1,p-1)


X.nostack = matrix(rnorm(n*k),n,k)
beta = matrix(rnorm(k*(p-1)),ncol=1)
w.vec = matrix(rmvnorm(1,mean=Augment.X(X.nostack,p)%*%beta,sigma=H%x%TT),ncol=1)
locs = sort(rep(1:n,p-1)) # This tells me how to group together the w's
y = apply(array(w.vec,c(p-1,1,n)),3,which.max)*apply(array(w.vec,c(p-1,1,n)),3,function(c)max(c)>0)
table(y)

cols = rainbow(p)
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
text(my.grid[,1],my.grid[,2],labels=y,col=cols[y+1])

system.time({ my.fit = Unordered.Probit.GP(X=X.nostack,y=y,my.grid=my.grid,step=.1,lambda.UB=sqrt(3),b0.scalar=0,
                             B.inv.scalar=.1,maxIt=2) })

par(mfrow=c(3,2))
for(i in 1:length(beta))
{
  plot(my.fit$B.fit[,i],type="l")
  abline(h=beta[i],col="chartreuse")
}
