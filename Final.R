## GEOG 418 Final Project 
## Andie Wood V00844811

# #Install packages as needed
# install.packages("sf")
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("spdep")
# install.packages("GISTools")
# install.packages("raster")
# install.packages("maptools")
# install.packages("rgdal")
# install.packages("spatstat")
# install.packages("sp")
# install.packages("tmap")
# install.packages("gstat")
# install.packages("gtable")

#Libraries
library(sf)
library(plyr)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(spatstat)
library(sp)
library(spatstat)
library(tmap)
library(gstat)
library(gtable)
library(spgwr)
library(grid)
library(gridExtra)
library(bcmaps)
library(bcmapsdata)

#########################
## Load and Clean Data ##
#########################

#Set working directory
setwd("C:/Users/Andie/Desktop/GEOG418/FinalProject")

get_layer("bc_bound", class = "sp")
get_layer("bc_cities", class = "sp")

bcmap <- tm_shape(bc_bound) + 
  tm_fill()+ 
  tm_layout(title = "British Columbia, Canada", 
            legend.position = c("left", "bottom"))+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"))+
  tm_compass(position = c("right", "top"))

png("Figures/BCMap.png")
bcmap
dev.off()

##Reading in particulate matter dataset
pm25 <- read.csv("PM25.csv") #Read in PM2.5 data
#Select only columns 1 and 2
pm25 <- pm25[,1:2]
#Change the column names, ignore N/A values 
colnames(pm25) <- c("POSTALCODE", "PM25")
pm25 <- na.omit(pm25)

##Reading in postal code shapefile
postalcodes <- shapefile("BC_Postal_Codes") #Read in related postal code data

##Reading in dissemination tract and income data
income <- read.csv("Income.csv") #Read in census income data  
colnames(income) <- c("DAUID", "Income") #Select only ID and Income columns
census.tracts <- shapefile("BC_DA.shp") #Read in dissemination tract shapefile
income.tracts <- merge(census.tracts,income, by = "DAUID") #Merge income and dissemination data
nrow(income.tracts) #Determine the number of columns in the dataframe
income.tracts <- income.tracts[!is.na(income.tracts$Income),]

map_MedIncome <- tm_shape(income.tracts) +
  tm_polygons(col = "Income", 
              title = "Median Income", 
              style = "jenks", 
              palette = "RdPu", 
              n = 6, 
              border.alpha = 0.5) +
  tm_layout(title = "Median Incomes by Census Tract in GVRD - 2015", 
            legend.position = c("left", "bottom"))+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"))+
  tm_compass(position = c("right", "top"))

png("Figures/MedianIncomeMap.png")  
map_MedIncome
dev.off()

#Select postal codes that fall within dissemination tracts)
postalcodes <- intersect(postalcodes,income.tracts)

png("Figures/PostalCode")
plot(postalcodes) #See what the data looks like spatially
dev.off()
head(postalcodes) #See what the data looks like in tabular form

#Join PM2.5 data with postal code data
pm25.spatial <- merge(postalcodes,pm25,by = "POSTALCODE")

#Aggregate the PM2.5 values in each DA in order to have a single value per DA. Here we aggregate based on the max.
pm25.aggregate <- aggregate((as.numeric(pm25.spatial$PM25)/10)~pm25.spatial$DAUID,FUN=max)

#Re-join aggregated data to the income.tracts layer.
colnames(pm25.aggregate) <- c("DAUID", "PM25AGG") #Select only ID and Income columns
income.pm25 <- merge(income.tracts,pm25.aggregate, by = "DAUID") #Merge income and dissemination data

#Re-join aggregated data to the pm25.spatial points layer.
pm25.points.aggregate <- merge(pm25.spatial, pm25.aggregate, by = "DAUID")

##Create a subsample of the datapoints provided in the PM2.5 dataset
## Sample size = I was suppose to be 410 but my computer would crash waiting for things
## to run, so I chose to half it. 

sampleSize=205
spSample <- pm25.points.aggregate[sample(1:length(pm25.points.aggregate),sampleSize),]
spSample2 <- pm25.points.aggregate[sample(1:length(pm25.points.aggregate),100),]


#Create a grid called grd to use in your interpolation

grd <- as.data.frame(spsample(spSample, "regular", n=5000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(spSample)

###############################################
## 1. Spatial Autocorrelation on Income Data ##
###############################################

## Create a neighbourhood list, Queens method = true
gvrd.nb <- poly2nb(income.tracts)
## Create neighbourhood net
gvrd.net <- nb2lines(gvrd.nb, coords = coordinates(income.tracts))

## Create lagged weights matrix from neighbourhood list
gvrd.lw <- nb2listw(gvrd.nb, zero.policy = TRUE, style = "W")
print.listw(gvrd.lw, zero.policy = TRUE)

## Global Moran's I 
mi <- moran.test(income.tracts$Income, gvrd.lw, zero.policy = TRUE)
mi

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(gvrd.lw)

mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- ((mI - eI) / sqrt(var))
z 

# Use z value and 95% significance to determine spatial pattern 
if( z > 1.96 | z < -1.96){
  # z is significant
  if(mI > eI){
    spat_pattern = "Clustered"
  }else{
    spat_pattern = "Dispersed"
  }
}else{
  # z is not significant
  spat_pattern = "Random"
}
spat_pattern

## Local Moran's I 
lisa.test <- localmoran(income.tracts$Income, gvrd.lw)

income.tracts$Ii <- lisa.test[,1]
income.tracts$E.Ii<- lisa.test[,2]
income.tracts$Var.Ii<- lisa.test[,3]
income.tracts$Z.Ii<- lisa.test[,4]
income.tracts$P<- lisa.test[,5]

map_LISA <- tm_shape(income.tracts) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fixed",
              palette = "-RdYlBu",
              midpoint = NA, 
              n = 3,
              breaks = c(-Inf, -0.01, 0.01, Inf), 
              labels = c("Negative", "Random", "Positive"),
              border.alpha = 0.5)+
  tm_layout(title = "Local Moran's I for Income in GVRD Census Tracts", 
            legend.position = c("left", "bottom"))+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"))+
  tm_compass(position = c("right", "top"))
png("Figures/LocalMorans.png")
map_LISA
dev.off()

## Plot Local Moran's I against lagged values
png("Figures/LMPlot.png")
moran.plot(income.tracts$Income, gvrd.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Income", 
           ylab="Spatially Lagged Income", quiet=NULL)
dev.off()

##############################
## 2. Spatial Interpolation ## 
##############################

## IDW (Lab 4)
P.idw <- gstat::idw(PM25AGG ~ 1, spSample, newdata=grd, idp=4)
r       <- raster(P.idw)
r.m     <- mask(r, income.tracts)

#map IDW
mapIDW <- tm_shape(r.m) + 
  tm_raster(n=10,palette = "RdPu",alpha = 0.8,
            title="Predicted PM2.5\n(in ug/m3)") + 
  tm_shape(spSample) + 
  tm_dots(size=0.05) +
  tm_legend(legend.outside=TRUE)+
  tm_scale_bar(position = "right") +
  tm_layout(main.title = "Inverse Distance Weighted Spatial Interpolation", main.title.size = 0.8)

png("Figures/IDW.png")
mapIDW
dev.off()

# Leave-one-out validation routine
IDW.out <- vector(length = length(spSample))
for (i in 1:length(spSample)) {
  IDW.out[i] <- gstat::idw(PM25AGG ~ 1, spSample[-i,], spSample[i,], idp=4)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))

png("Figures/LeaveOneOut4.png")
plot(IDW.out ~ spSample$PM25AGG, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ spSample$PM25AGG), col="red", lw=2,lty=2)
abline(0,1)
par(OP)

dev.off()

#RMSE
sqrt( sum((IDW.out - spSample$PM25AGG)^2) / length(spSample))

# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(PM25AGG~1, spSample, newdata=grd, idp=3)
n   <- length(spSample)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(PM25AGG~1, spSample[-i,], newdata=grd, idp=3)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}


# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Southern California
r <- raster(img.sig, layer="v")
r.m <- mask(r, income.tracts)

# Plot the map
tm_shape(r.m) + tm_raster(n=7,title="95% confidence interval \n(in ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.1) +
  tm_legend(legend.outside=TRUE)

### Trend Surface Analysis

# First Order 
f.1 <- as.formula(PM25AGG ~ X + Y) 

# Add X and Y to P - predict attribute based on x and y coordinates
spSample$X <- coordinates(spSample)[,1]
spSample$Y <- coordinates(spSample)[,2]

# Run the regression model (lm is a linear regression model)
lm.1 <- lm( f.1, data=spSample)

# Use the regression model output to interpolate the surface
dat.1st <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.1, newdata=grd))) 

# Clip the interpolated raster to Study Area
r   <- raster(dat.1st)
r.m <- mask(r, income.tracts)

# Plot the map
firstOrder <- tm_shape(r.m) + 
  tm_raster(n=10, palette="RdPu", 
            title="Predicted PM2.5 \n(ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.05) +
  tm_layout(main.title = "First Order: Linear Trend ",
            main.title.size = 1) +
  tm_legend(legend.outside=TRUE)

png("Figures/FirstOrder.png")
firstOrder
dev.off()

## Second Order
f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))
# Add X and Y to P - predict attribute based on x and y coordinates
spSample$X <- coordinates(spSample)[,1]
spSample$Y <- coordinates(spSample)[,2]

# Run the regression model (lm is a linear regression model)
lm.2 <- lm( f.2, data=spSample)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

# Clip the interpolated raster to Study Area
r   <- raster(dat.2nd)
r.m <- mask(r, income.tracts)

# Plot the map
secondOrder <- tm_shape(r.m) + 
  tm_raster(n=10, palette="RdPu", 
            title="Predicted PM2.5 \n(ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.05) +
  tm_layout(main.title = "Second Order: Quadratic Trend ",
            main.title.size = 1) +
  tm_legend(legend.outside=TRUE)

png("Figures/SecondOrder.png")
secondOrder
dev.off()

### Universal Kriging with 2nd Order Polynomial
#Find best fit by adjusting sill, model, range and nugget
var.smpl <- variogram(f.2, spSample, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=0.75, model="Sph", range=8, nugget=0))

png("Figures/Semivari.png")
plot(var.smpl, dat.fit)
dev.off()

# Define the trend model
f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg <- krige( f.2, spSample, grd, dat.fit)
colnames(dat.krg@data)[1] <- "PM25"
### This krige wont work and I think it's cause I have duplicates but I can't figure out why

r <- raster(dat.krg)
r.m <- mask(r, income.tracts)

# Plot the map
mapKrig <- tm_shape(r.m) + 
  tm_raster(n=10, palette="RdPu",  
            title="Predicted PM2.5 \n(ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.05) +
  tm_legend(legend.outside=TRUE)+
  tm_layout(main.title = "Universal Kriging", main.title.size = 1)

png("Figures/MapKrig.png")
mapKrig
dev.off()

r   <- raster(dat.krg, layer="var1.var")
r.mv <- mask(r, income.tracts)

krigVar <- tm_shape(r.mv) + 
  tm_raster(n=10, palette="RdPu",  
            title="Variance map \n(in squared ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.05) +
  tm_legend(legend.outside=TRUE)+
  tm_layout(main.title = "Universal Kriging - Variance", main.title.size = 1)

png("Figures/KrigVar.png")
krigVar
dev.off()

r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.mc <- mask(r, income.tracts)

krigConf <- tm_shape(r.mc) + 
  tm_raster(n=7, palette ="RdPu",
            title="95% CI map \n(in ug/m3)") +tm_shape(spSample) + tm_dots(size=0.05) +
  tm_legend(legend.outside=TRUE)+
  tm_layout(main.title = "Universal Kriging \n95% Confidence Intervals", main.title.size = 1)

png("Figures/KrigConf.png")
krigConf
dev.off()

###################################
## 3. Linear Regression Analysis ##
###################################

#If you have too many cells, you can reduce the number by aggregating values
step.1 <- aggregate(r.m, fact=1, fun=mean)
plot(step.1)

#Convert the raster dataset to points
step.2 <-  rasterToPoints(step.1,fun=NULL, spatial=FALSE, crs=census.tracts)
step.2 <- as.data.frame(step.2) #convert the point dataset to a spatial dataframe
Coords <- step.2[,c("x", "y")]  #assign coordinates to a new object
crs <- crs(census.tracts) #utilize an existing projection
step.3 <- SpatialPointsDataFrame(coords = Coords, data = step.2, proj4string = crs) #create a spatial points dataframe
step.4 <- aggregate(x=step.3,by=income.tracts, FUN=mean) #aggregate points into census tracts
step.5 <- intersect(step.4,income.tracts)  #get the intersection of step.4 with the income.tracts dataset (this will take a while) 

#Now regression
#Dataset with both PM2.5 and Income are stored in a dataset called pm.income.poly.

pm.income.poly <- step.5
colnames(pm.income.poly@data)[3] <- "PM25"

#Plot income and PM2.5 from the pm.income.poly dataset you created
plot(pm.income.poly$PM25~pm.income.poly$Income)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$PM25),]
pm.income.poly <- pm.income.poly[pm.income.poly$PM25 != 0, ]
#Now plot the data again

png("Figures/LinRegPlot.png")
plot(pm.income.poly$PM25~pm.income.poly$Income, ylab = "PM2.5", xlab = "Income")

#Perform a linear regression on the two variables.
#Pollution depends on income 
lm.model <- lm(pm.income.poly$PM25~pm.income.poly$Income)
#Income depends on pollution 
#lm.model <- lm(pm.income.poly$Income~pm.income.poly$PM25)

#Add the regression model to the plot you created
abline(lm.model)
dev.off()
#Get the summary of the results
summary(lm.model)


#You want to determine if the model residuals are spatially clustered. 
#First obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(lm.model))
#Then add the residuals to your spatialpolygon dataframe
pm.income.poly$residuals <- residuals.lm(lm.model)
#Observe the result to make sure it looks correct
head(pm.income.poly)

mapResid <- tm_shape(pm.income.poly)+
  tm_polygons(col = "residuals",
              style = "jenks", 
              palette = "RdPu",
              n = 6, 
              midpoint = NA)+
  tm_legend(legend.outside=TRUE)+
  tm_layout(main.title = "Regression Residuals, PM2.5")
  
png("Figures/MapResid.png")
mapResid
dev.off()

###################################
## 4. Residuals (Moran's I again)##
###################################

## Create a neighbourhood list, Queens method = true
resid.nb <- poly2nb(pm.income.poly)
## Create neighbourhood net
resid.net <- nb2lines(resid.nb, coords = coordinates(pm.income.poly))

## Create lagged weights matrix from neighbourhood list
resid.lw <- nb2listw(resid.nb, zero.policy = TRUE, style = "W")
print.listw(resid.lw, zero.policy = TRUE)

## Global Moran's I 
mi <- moran.test(pm.income.poly$residuals, resid.lw, zero.policy = TRUE)
mi

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(resid.lw)

mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- ((mI - eI) / sqrt(var))
z 

# Use z value and 95% significance to determine spatial pattern 
if( z > 1.96 | z < -1.96){
  # z is significant
  if(mI > eI){
    spat_pattern = "Clustered"
  }else{
    spat_pattern = "Dispersed"
  }
}else{
  # z is not significant
  spat_pattern = "Random"
}
spat_pattern

## Local Moran's I 
relisa.test <- localmoran(pm.income.poly$residuals, resid.lw)

pm.income.poly$Ii <- relisa.test[,1]
pm.income.poly$E.Ii<- relisa.test[,2]
pm.income.poly$Var.Ii<- relisa.test[,3]
pm.income.poly$Z.Ii<- relisa.test[,4]
pm.income.poly$P<- relisa.test[,5]

remap_LISA <- tm_shape(pm.income.poly) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fixed",
              palette = "-RdYlBu",
              midpoint = NA, 
              n = 3,
              breaks = c(-Inf, -0.01, 0.01, Inf), 
              labels = c("Negative", "Random", "Positive"),
              border.alpha = 0.5)+
  tm_layout(title = "Local Moran's I for PM2.5 Residuals", 
            legend.position = c("left", "bottom"))+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"))+
  tm_compass(position = c("right", "top"))

png("Figures/LocalMoransResid.png")
remap_LISA
dev.off()

## Plot Local Moran's I against lagged values
png("Figures/LMPlotResid.png")
moran.plot(pm.income.poly$residuals, resid.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Residuals", 
           ylab="Spatially Lagged Residuals", quiet=NULL)
dev.off()


#############
## 5. GWR  ##
#############
#Let's say you are continuing with your data from the regression analysis. 

#The first thing you need to do is to add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
pm.income.poly.coords <- sp::coordinates(pm.income.poly)
#Observe the result
head(pm.income.poly.coords)
#Now add the coordinates back to the spatialpolygondataframe
pm.income.poly$X <- pm.income.poly.coords[,1]
pm.income.poly$Y <- pm.income.poly.coords[,2]
head(pm.income.poly)


###Could not get GWR to work from here on out. Always crashed or had error :(

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(pm.income.poly$PM25~pm.income.poly$Income, 
                        data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(pm.income.poly$PM25~pm.income.poly$Income, 
                data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
pm.income.poly$localr <- results$localR2

# #Create choropleth map of r-square values
# local.r.square <- pm.income.poly$localr
# shades <- auto.shading(local.r.square, n=6, cols = brewer.pal(6, 'Oranges'))
# choropleth(income.tracts, local.r.square, shades) #map the data with associated colours
# choro.legend(3864000, 1965000, shades) #add a legend (you might need to change the location)

mapLocalr2 <- tm_shape(pm.income.poly) + 
  tm_polygons(col = "localr",
              title = "Local R2",
              style = "jenks",
              palette = "RdPu",
              midpoint = NA, 
              n = 6,
              border.alpha = 0.35)+
  tm_layout(legend.position = c("left", "bottom"), 
            title = "Local R2 of PM2.5 for GWR")+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"))+
  tm_compass(position = c("right", "top"))

mapLocalr2

#Time for more magic. Let's map the coefficients
pm.income.poly$coeff <- results$pm.income.poly.PM25

# #Create choropleth map of the coefficients
# local.coefficient <- pm.income.poly$coeff
# shades <- auto.shading(local.coefficient, n=6, cols = brewer.pal(6, 'Oranges'))
# choropleth(income.tracts, local.coefficient, shades) #map the data with associated colours
# choro.legend(3864000, 1965000, shades) #add a legend (you might need to change the location)

mapCoeff <- tm_shape(pm.income.poly) + 
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              palette = "RdPu",
              n = 6,
              border.alpha = 0.35)+
  tm_layout(legend.position = c("left", "bottom"), 
            title = "GWR Coefficients for PM2.5")+
  tm_scale_bar(position = c("RIGHT", "BOTTOM"))+
  tm_compass(position = c("right", "top"))

mapCoeff


###############################
## 6. Point Pattern Analysis)##
###############################

#####
###KERNEL DENSITY ESTIMATION
#2D (gaussian) kernel, compare how bandwidth (sigma) selection influences the point density estimates
#since data are projected, sigma is represented in metres
#eps is the width and height of the pixels (1000m X 1000m)
#coerce to a SpatialGridDataFrame for plotting
# 
# kde.100 <- density(kma.ppp, sigma = 100, at = "pixels", eps = c(1000, 1000))
# kde.SG <- as(kde.100, "SpatialGridDataFrame")
# kde.500 <- density(kma.ppp, sigma = 500, at = "pixels", eps = c(1000, 1000))
# kde.SG <- cbind(kde.SG, as(kde.500, "SpatialGridDataFrame"))
# 
#   names(kde.SG) <- c(***NAME YOUR SENSITIVITY TEST COLUMNS***)
# #plot
# x11() #opens a new plot window
# spplot(kde.SG)
# 
# #can see how the bandwidth selection influences the density estimates
# summary(kde.SG)
# 
# #use cross-validation to get the bandwidth that minimizes MSE
# bw.d <- bw.diggle(kma.ppp)
# #plot the "optimal" bandwidth
# plot(bw.d, ylim=c(-10, 10), main= ***CHOOSE AN APPROPRIATE TITLE***)
# 
# #density using the cross-validation bandwidth
# kde.bwo <- density(kma.ppp, sigma = bw.d, at = "pixels", eps = c(1000, 1000))
# plot(kde.bwo)
# 


