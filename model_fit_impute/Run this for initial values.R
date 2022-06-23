load("~/Dropbox/Leman Research/Poison Ivy Spatial/04_22_17_cat_run_complete.Rdata")
w.initial = w.fit.store[nrow(w.fit.store),]
beta.initial = beta.store[nrow(beta.store),]
lambda.initial = lambda.store[length(lambda.store)]
TT.initial = TT.curr
rm(list= ls()[!(ls() %in% c("w.initial","beta.initial","lambda.initial","TT.initial"))])

beta.initial = beta.initial[-c(5,55,105)] # Removing FlowAcc
length(beta.initial)
