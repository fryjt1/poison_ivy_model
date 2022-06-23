library(glmnet)
# Fit a model for Bedrock Depth missing data

missing.cols = c("BedrockDepth","WaterStorage","soiltexture.new")
MM = model.matrix(~.-1,data=model.data[,!colnames(model.data)%in%missing.cols])[,-1]

miss.ind = which(is.na(model.data$BedrockDepth))
bedrock.miss.ind = miss.ind
MM.BedrockDepth = MM[-miss.ind,]
y = model.data$BedrockDepth[-miss.ind]
cv = cv.glmnet(x=MM.BedrockDepth,y=y,alpha=1,nfolds=5)
plot(cv)
Bedrock.lasso = glmnet(x=MM.BedrockDepth,y=y,standardize=TRUE,lambda=cv$lambda.1se)
MM.new = MM.BedrockDepth[,-which(Bedrock.lasso$beta==0)]
Bedrock.fit = lm(y~.,data=as.data.frame(MM.new))
newX = MM[miss.ind,-which(Bedrock.lasso$beta==0)]
temp = predict(Bedrock.fit,as.data.frame(newX),se.fit=TRUE)
m0.bedrock = temp$fit
C0.bedrock = 3*diag(temp$se.fit)

# Fit a model for Water Storage missing data

missing.cols = c("BedrockDepth","WaterStorage","soiltexture.new","Soil.Series")
MM = model.matrix(~.-1,data=model.data[,!colnames(model.data)%in%missing.cols])[,-1]

miss.ind = which(is.na(model.data$WaterStorage))
water.miss.ind = miss.ind
MM.WaterStorage = MM[-miss.ind,]
y = model.data$WaterStorage[-miss.ind]
cv = cv.glmnet(x=MM.WaterStorage,y=y,alpha=1,nfolds=5)
plot(cv)
Water.lasso = glmnet(x=MM.WaterStorage,y=y,standardize=TRUE,lambda=cv$lambda.1se)
MM.new = MM.WaterStorage[,-which(Water.lasso$beta==0)]
WaterStorage.fit = lm(y~.,data=as.data.frame(MM.new))
newX = MM[miss.ind,-which(Water.lasso$beta==0)]
temp = predict(WaterStorage.fit,as.data.frame(newX),se.fit=TRUE)

m0.water = temp$fit
C0.water = 3*diag(temp$se.fit)

# Find proportions of each soiltexture
miss.ind = which(is.na(model.data$soiltexture.new))
soil.miss.ind = miss.ind
soil.prop = table(model.data$soiltexture.new[-miss.ind])/sum(table(model.data$soiltexture.new[-miss.ind]))
