library(raster)
library(sp)
library(rasterVis)
library(hydroGOF)

source('R/functions.R')

load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")

## Build a brick containing all data
alldata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
names(alldata) <- c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")

##remove outliers
#hist(alldata)
alldata$band1[alldata$band1 < 0 | alldata$band1 > 750] <- NA
alldata$band2[alldata$band2 < 0 | alldata$band2 > 1000] <- NA
alldata$band3[alldata$band3 < 0 | alldata$band3 > 1250] <- NA
alldata$band4[alldata$band4 < 0 | alldata$band4 > 4000] <- NA
alldata$band5[alldata$band5 < 0 | alldata$band5 > 3500] <- NA
alldata$band7[alldata$band7 < 0 | alldata$band7 > 2500] <- NA
alldata$VCF[alldata$VCF < 0 | alldata$VCF > 100] <- NA
par(mfrow=c(1,1))

## plot the relationship between the Landsat bands and the VCF tree cover
#pairs(alldata)
par(mfrow=c(1,1))

## Extract all data to a data.frame
df <- as.data.frame(getValues(alldata))

##define linear model
linmod <- lm(VCF~band1+band2+band3+band4+band5+band7, df)
summary(linmod)

##predict values
predictVCF <- predict(alldata[[c(1:6)]], linmod)
names(predictVCF) <- "Predicted_VCF"
predictVCF[predictVCF < 0 | predictVCF > 100] <- NA 

##plot and compare both tree cover maps
VCFrasters <- brick(predictVCF, alldata$VCF)
#levelplot(VCFrasters, col.regions=rev(terrain.colors(255)))

##difference map
diff <- calc(VCFrasters, difference)
#plot(diff, main='Difference between predicted and observed tree cover')

##compute RMSE
df_VCF <- as.data.frame(getValues(VCFrasters))
rmse <- rmse(df_VCF$Predicted_VCF, df_VCF$VCF, na.rm=T)

##calculate RMSE for each different land cover class
load('data/trainingPoly.rda')
classes <- rasterize(trainingPoly, predictVCF, field='Class', na.rm=T)
diffsq <- calc(VCFrasters, sq_diff)
avgdiffsq <- zonal(diffsq, classes, fun='mean', na.rm=T)
df_rmse_classes <- data.frame(avgdiffsq)
df_rmse_classes$rmse <- sqrt(df_rmse_classes$mean)
df_rmse_classes$class <- c('cropland', 'forest', 'wetland')
