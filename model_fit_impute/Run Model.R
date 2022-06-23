library(truncnorm)
library(MCMCpack)
library(mvtnorm)
library(Matrix)
library(psych)
library(abind)

X.nostack.miss = design.mat

b0.scalar = 0
B.inv.scalar = .1
lambda.UB = sqrt(2)
step = .01

I = 250000
X = Augment.X(X.nostack.miss,p)
y = PI.type
b0 = matrix(rep(b0.scalar,k*(p-1)),ncol=1)
B.inv = diag(B.inv.scalar,k*(p-1))

kappa = (p-1) + 2  # Priors on TT
C = diag(1,p-1)

# mixed.ind = which(colnames(design.mat)=="LC_describe.newMixed Forest")
# mixed.ind = c(mixed.ind,mixed.ind+50,mixed.ind+100)
# temp = c(beta.temp)
# beta.initial = append(temp,NA,after=mixed.ind[1]-1)
# beta.initial[which(is.na(beta.initial))] = beta.initial[which(is.na(beta.initial))-2]
# beta.initial = append(beta.initial,NA,after=mixed.ind[2]-1)
# beta.initial[which(is.na(beta.initial))] = beta.initial[which(is.na(beta.initial))-2]
# beta.initial = append(beta.initial,NA,after=mixed.ind[3]-1)
# beta.initial[which(is.na(beta.initial))] = beta.initial[which(is.na(beta.initial))-2]

# Initialize

locs = sort(rep(1:n,p-1))
w.fit.vec = matrix(w.initial,ncol=1)
w.fit.store = matrix(NA,I,length(w.fit.vec)); w.fit.store[1,] = w.fit.vec
beta.store = matrix(NA,I,k*(p-1)); beta.store[1,] = beta.initial
lambda.store = rep(NA,I); lambda.store[1] = lambda.initial
TT.store = matrix(NA,nrow=I,ncol=(p-1)^2); TT.store[1,] = c(TT.initial)
water.store = matrix(0,I,length(water.miss.ind))
bedrock.store = matrix(0,I,length(bedrock.miss.ind))
soil.store = matrix(0,I,length(soil.miss.ind))
beta.curr = matrix(beta.store[1,],ncol=1)
lambda.curr = lambda.store[1]
XB = X%*%beta.curr
TT.curr = TT.initial
lambda.curr = lambda.store[1]
H.curr = exp(-dist.mat/lambda.curr) + diag(1e-8,n)
H.inv.curr = solve(H.curr)

system.time({
for(i in 2:I)
{
  S.inv = H.inv.curr%x%solve(TT.curr)

  # Update Latent W's

  w.fit.vec = Sample.W(y=y,locs=locs,w.vec.curr=w.fit.vec,S.inv=S.inv,XB=XB)
  w.fit.store[i,] = w.fit.vec

  # Update Beta

  beta.curr = Sample.Beta(X=X,S.inv=S.inv,w.vec=w.fit.vec,B.inv=B.inv,b0=b0)
  beta.store[i,] = beta.curr/sqrt(TT.curr[1,1])
  XB = X%*%beta.curr

  # Update TT

  draw.list = Sample.TT(w.vec=w.fit.vec,XB=XB,H.inv=H.inv.curr,C=C,kappa=kappa)
  TT.curr = draw.list$draw
  matrix.SSE = draw.list$matrix.SSE
  TT.store[i,] = c(TT.curr)
  
  # Update lengthscale

  draw = Sample.Lambda(lambda.curr=lambda.curr,lambda.UB=lambda.UB,step=step,H.curr=H.curr,
                H.inv.curr=H.inv.curr,TT.curr=TT.curr,w.vec=w.fit.vec,XB=XB,matrix.SSE=matrix.SSE)
  lambda.curr = draw$lambda
  lambda.store[i] = lambda.curr
  H.curr = draw$H
  H.inv.curr = draw$H.inv
  
  # Impute missing values for water storage

  col.nums = which(colnames(X.nostack.miss)=="WaterStorage")
  draw = Impute.Cont(H.inv=H.inv.curr,TT.inv=solve(TT.curr),w.vec=w.fit.vec,
                     X.nostack.miss=X.nostack.miss,B=beta.curr,miss.ind=water.miss.ind,
                     col.nums=col.nums,m0=m0.water.scale, C0=C0.water.scale)
  water.store[i,] = draw
  X.nostack.miss[water.miss.ind,col.nums] = draw
  X = Augment.X(X.nostack.miss,p)

  # Impute missing values for bedrock depth

  col.nums = which(colnames(X.nostack.miss)=="BedrockDepth")
  draw = Impute.Cont(H.inv=H.inv.curr,TT.inv=solve(TT.curr),w.vec=w.fit.vec,
                     X.nostack.miss=X.nostack.miss,B=beta.curr,miss.ind=bedrock.miss.ind,
                     col.nums=col.nums,m0=m0.bedrock.scale, C0=C0.bedrock.scale)
  bedrock.store[i,] = draw
  X.nostack.miss[bedrock.miss.ind,col.nums] = draw
  X = Augment.X(X.nostack.miss,p)

  # Impute missing values for soil type
  col.nums = which(colnames(X.nostack.miss)%in%c("soiltexture.newFine Sandy Loam","soiltexture.newLoam",
                                                 "soiltexture.newRock Outcrop","soiltexture.newSandy Loam"))
  draw = Impute.Cat(H.inv=H.inv.curr,TT.inv=solve(TT.curr),w.vec=w.fit.vec,X.nostack.miss=X.nostack.miss,
                    B=beta.curr,miss.inds=soil.miss.ind,col.nums=col.nums)
  X.nostack.miss = draw$X
  X = Augment.X(draw$X,p)
  soil.store[i,] = draw$cat.ind
  
  print(i)
  
  if(i%%25000==0)
  {
    par(mfrow=c(3,3))
    for(k in 1:9)
    {
      plot(beta.store[20:i,k],type="l")
    }
  }
  
  if(i%%25000==0)
  {
    file = paste("09_26_17_cat_run_",i,sep="")
    path = paste("~/Dropbox/Leman Research/Poison Ivy Spatial/",file,".Rdata",sep="")
    save.image(file = path)
  }
}
  
})

# save.image(file = "~/Dropbox/Leman Research/Poison Ivy Spatial/09_19_17_cat_run_complete.Rdata")

