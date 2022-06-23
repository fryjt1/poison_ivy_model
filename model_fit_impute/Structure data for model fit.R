library(mvtnorm)
library(MCMCpack)
library(plgp)

# These are columns I'm choosing to remove after some variable selection
remove.cols = c("Soil.SeriesMarbleyard","Soil.SeriesOriskany",
                "Soil.SeriesShelocta","soiltexture.newSilt Loam",
                "soiltexture.newSilty Clay Loam")

# List the columns with missing data and find row index for missing data
missing.cols = c("BedrockDepth","WaterStorage","soiltexture.new")
missing.Bedrock = which(is.na(model.data$BedrockDepth))
missing.Water = which(is.na(model.data$WaterStorage))
missing.soiltexture = which(is.na(model.data$soiltexture.new))

# List off the categories (in the correct order) for the missing categorical variables
soiltexture.list = names(table(IvyData$soiltexture.new))
soiltexture.list = soiltexture.list[-c(6,7)]

# Copy the model data and populate the missing values with some initial values
model.data.c = model.data
model.data.c$BedrockDepth[missing.Bedrock] = m0.bedrock # Initializing missing values
model.data.c$WaterStorage[missing.Water] = m0.water # Initializing missing values
model.data.c$soiltexture.new[missing.soiltexture] = "Rock Outcrop"

# Convert the current model data into a model matrix
y = model.data$PI.type
design.mat = model.matrix(~.,data=model.data.c)[,-c(1:2)]

# Use latitude, longitude, elevation (in unit cube) to create a distance matrix
my.grid = design.mat[,1:3]
my.grid = cbind(design.mat[,1:3])
mins = apply(my.grid,2,function(c) min(c))
ranges = apply(my.grid,2,function(c) (max(c)-min(c)))
my.coord = scale(my.grid,mins,ranges)
dist.mat = as.matrix(dist(my.grid))

# Remove columns for non-selected variables
remove.ind = !colnames(design.mat)%in%remove.cols
design.mat = design.mat[,remove.ind]

# Find the column means and sds (not for categorical)
X.means = colMeans(design.mat)
X.sds = apply(design.mat,2,sd)
X.ranges = apply(design.mat,2,function(c) max(c)-min(c))
X.means[X.ranges==1] = 0
X.sds[X.ranges==1] = 1

# Scale continuous columns, but not categorical
design.mat = scale(design.mat,X.means,X.sds)
n = nrow(design.mat)
design.mat = cbind(rep(1,n),design.mat); colnames(design.mat)[1] = "Intercept"
k = ncol(design.mat)
p = max(y) + 1

m0.bedrock.scale = (m0.bedrock-X.means[which(names(X.means)=="BedrockDepth")])/
                    X.sds[which(names(X.sds)=="BedrockDepth")]
C0.bedrock.scale = C0.bedrock/X.sds[which(names(X.sds)=="BedrockDepth")]^2

m0.water.scale = (m0.water-X.means[which(names(X.means)=="WaterStorage")])/
  X.sds[which(names(X.sds)=="WaterStorage")]
C0.water.scale = C0.water/X.sds[which(names(X.sds)=="WaterStorage")]^2


