################################################################################################################################
# Import and Merge the 2 Datasets
################################################################################################################################

IvyValues_PRISM = read.csv("~/Dropbox/Leman Research/Poison Ivy Spatial/IvyValues_PRISM.csv",sep=",",strip.white=TRUE)
IvyValues_all = read.csv("~/Dropbox/Leman Research/Poison Ivy Spatial/IvyValues_all_Revised.csv",sep=",",strip.white=TRUE)
IvyValues_PRISM = IvyValues_PRISM[ ,!(names(IvyValues_PRISM)) %in% names(IvyValues_all)[-2]]
IvyData = merge(IvyValues_all,IvyValues_PRISM,by="id")

################################################################################################################################
# Subset the Temperature,Precipitation, and Solar Radiation Columns
################################################################################################################################

solRad.cols = IvyData[,colnames(IvyData)%in%c("SolRad_JAN","SolRad_FEB","SolRad_MAR","SolRad_APR","SolRad_MAY","SolRad_JUN",
                                                      "SolRad_JUL","SolRad_AUG","SolRad_SEP","SolRad_OCT","SolRad_NOV","SolRad_DEC")]
tmean.cols = IvyData[,colnames(IvyData)%in%c("tmean_JAN","tmean_FEB","tmean_MAR","tmean_APR","tmean_MAY","tmean_JUN",
                                                      "tmean_JUL","tmean_AUG","tmean_SEP","tmean_OCT","tmean_NOV","tmean_DEC")]
tmin.cols = IvyData[,colnames(IvyData)%in%c("tmin_JAN","tmin_FEB","tmin_MAR","tmin_APR","tmin_MAY","tmin_JUN",
                                            "tmin_JUL","tmin_AUG","tmin_SEP","tmin_OCT","tmin_NOV","tmin_DEC")]
tmax.cols = IvyData[,colnames(IvyData)%in%c("tmax_JAN","tmax_FEB","tmax_MAR","tmax_APR","tmax_MAY","tmax_JUN",
                                                     "tmax_JUL","tmax_AUG","tmax_SEP","tmax_OCT","tmax_NOV","tmax_DEC")]
tmin.cols = IvyData[,colnames(IvyData)%in%c("tmin_JAN","tmin_FEB","tmin_MAR","tmin_APR","tmin_MAY","tmin_JUN",
                                                     "tmin_JUL","tmin_AUG","tmin_SEP","tmin_OCT","tmin_NOV","tmin_DEC")]
ppt.cols = IvyData[,colnames(IvyData)%in%c("ppt_JAN","ppt_FEB","ppt_MAR","ppt_APR","ppt_MAY","ppt_JUN",
                                                    "ppt_JUL","ppt_AUG","ppt_SEP","ppt_OCT","ppt_NOV","ppt_DEC")]

################################################################################################################################
# Create some transformed variables
################################################################################################################################

IvyData$PI.ind = 1
IvyData$PI.ind[which(IvyData$plant_type=="A")] = 0
IvyData$PI.type = 0
IvyData$PI.type[which(IvyData$plant_type=="C")] = 1
IvyData$PI.type[which(IvyData$plant_type=="S")] = 2
IvyData$PI.type[which(IvyData$plant_type=="V")] = 3

# IvyData$PI.ind = factor(IvyData$PI.ind)
IvyData$Sin.Aspect = sin(IvyData$Aspect*pi/180)
IvyData$Cos.Aspect = cos(IvyData$Aspect*pi/180)
IvyData$tmean.growing = rowMeans(tmean.cols[,-c(1,2,11,12)])
IvyData$ppt.growing =  rowMeans(ppt.cols[,-c(1,2,11,12)])
IvyData$solRads.growing =  rowMeans(solRad.cols[,-c(1,2,11,12)])
# Rename some columns with ridiculously long names
colnames(IvyData)[which(colnames(IvyData)=="Available.Water.Storage.0.50.cm...Weighted.Average")] = "WaterStorage"
colnames(IvyData)[which(colnames(IvyData)=="Drainage.Class...Dominant.Condition")] = "DrainageClass"
colnames(IvyData)[which(colnames(IvyData)=="Bedrock.Depth...Minimum..in.")] = "BedrockDepth"

################################################################################################################################
# Code new groups for Land Cover
################################################################################################################################

IvyData$LC_describe.new = NA
IvyData$LC_describe.new[IvyData$LandCover%in%c(11,41)] = "Deciduous Forest"
IvyData$LC_describe.new[IvyData$LandCover%in%c(21,22,23)] = "Developed"
IvyData$LC_describe.new[IvyData$LandCover%in%c(42)] = "Evergreen Forest"
IvyData$LC_describe.new[IvyData$LandCover%in%c(43)] = "Mixed Forest"
IvyData$LC_describe.new[IvyData$LandCover%in%c(81,82)] = "Planted/Cultivated"

################################################################################################################################
# Code new groups for Soil Texture
################################################################################################################################

IvyData$Surface.Soil.Texture = toupper(IvyData$Surface.Soil.Texture)
IvyData$soiltexture.new = NA
IvyData$soiltexture.new[IvyData$Surface.Soil.Texture%in%c("VERY GRAVELLY SILT LOAM","VERY CHANNERY SILT LOAM","SILT LOAM",
                                                          "GRAVELLY SILT LOAM","CHANNERY SILT LOAM")] = "Silt Loam"
IvyData$soiltexture.new[IvyData$Surface.Soil.Texture%in%c("CHANNERY SILTY CLAY LOAM","SILTY CLAY LOAM")] = "Silty Clay Loam"
IvyData$soiltexture.new[IvyData$Surface.Soil.Texture%in%c("LOAM","GRAVELLY LOAM","COBBLY LOAM")] = "Loam"
IvyData$soiltexture.new[IvyData$Surface.Soil.Texture%in%c("DISTURBED")] = "Disturbed"
IvyData$soiltexture.new[IvyData$Surface.Soil.Texture%in%c("ROCK OUTCROP","ROCK OUTCROP COMPLEX")] = "Rock Outcrop"
IvyData$soiltexture.new[IvyData$Surface.Soil.Texture%in%c("SANDY LOAM","CHANNERY SANDY LOAM","COBBLY SANDY LOAM",
                                                          "VERY COBBLY SANDY LOAM","VERY GRAVELLY SANDY LOAM")] = "Sandy Loam"
IvyData$soiltexture.new[IvyData$Surface.Soil.Texture%in%c("CHANNERY FINE SANDY LOAM","COBBLY FINE SANDY LOAM",
                                                          "FINE SANDY LOAM")] = "Fine Sandy Loam"

################################################################################################################################
# Remove columns that aren't needed
################################################################################################################################

keep.cols = c("PI.type","latitude","longitude","Elevation","Curvature","SlopeDeg","Soil.Series","BedrockDepth",
              "WaterStorage","Sin.Aspect","Cos.Aspect","tmean.growing","ppt.growing","solRads.growing","LC_describe.new","soiltexture.new")

model.data = IvyData[-which.max(IvyData$FlowAcc),keep.cols]


# keep.cols = c("PI.type","latitude","longitude","Elevation","FlowAcc","Curvature","SlopeDeg","Soil.Series","BedrockDepth",
#               "WaterStorage","Sin.Aspect","Cos.Aspect","tmean.growing","ppt.growing","solRads.growing","LC_describe.new","soiltexture.new")
# model.data = IvyData[,keep.cols]
# model.data = model.data[-which.max(FlowAcc),]

attach(model.data)
# plot(NA,NA,ylim=c(min(latitude),max(latitude)),xlim=c(min(longitude),max(longitude)),xlab="Longitude",ylab="Latitude")
# points(longitude[PI.type==0],latitude[PI.type==0],cex=.65,col="black")
# points(longitude[PI.type==1],latitude[PI.type==1],cex=.65,col="red")
# points(longitude[PI.type==2],latitude[PI.type==2],cex=.65,col="blue")
# points(longitude[PI.type==3],latitude[PI.type==3],cex=.65,col="green")
# head(model.data)
# 
# rank(FlowAcc)
# r = rev(rank(FlowAcc))
# cols = heat.colors(length(FlowAcc))[r]
# plot(NA,NA,xlim=range(longitude),ylim=range(latitude),main="Flow Accumulation",xlab="Longitude",
#      ylab="Latitude")
# points(longitude,latitude,col=cols,pch=19,cex=.25)
# ind = which(FlowAcc%in%sort(FlowAcc)[(length(FlowAcc)-11):(length(FlowAcc))])
# text(longitude[ind],latitude[ind],labels=FlowAcc[ind],cex=.75)

