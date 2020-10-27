### Jennifer Wong-Ala
### Purpose: The purpose of this script is to compare and contrast the velocities particles experience for the entire simulation period

##########################################################
### clear workspace 

rm(list=ls())

##########################################################
### load libraries

library(marmap)
library(raster)
library(rgdal)

##########################################################
### load data 
setwd("/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/model_output_chapter1&2/used_high_res_ROMS/2016")
load('LPT_data_pp_06.13.2020_JWA_6hour_2016_high.RData')


# simulation used is 250 m ROMS ebb tide
# LPT_06.13.2020_JWA_6hour_2016_high.RData

head(data_use)
dim(data_use)

range(data_use$ID) # 1 750 

# data_use$ID[which(data_use$depth<=-700)] # figure out which particle is the really deep one

# plot(data_use$del.t, data_use$depth)

# subset out two particles 
## one that primarily stays in the surface and the second particle is the one that goes all the way to 600 m depths. 

par_shal<-subset(data_use, ID==50)
dim(par_shal) # 3360   16
head(par_shal)

par_deep2<-subset(data_use, ID==268)
dim(par_deep2) # 3432   16
head(par_deep2)

par_deep<-par_deep2[1:nrow(par_shal),]
dim(par_deep)

date<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(par_deep))
# date1<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(par_shal))

sum_deep<-sum(par_deep$w_vel, na.rm=T)*100 # units = m/s (-0.4186448)

sum_shal<-sum(par_shal$w_vel, na.rm=T)*100 # units = m/s (-0.07111202)

avg_shal<-mean(par_shal$w_vel, na.rm=T)*100

avg_deep<-mean(par_deep$w_vel, na.rm=T)*100

sum_shal2<-paste('shallow w_vel sum = ', round(sum_shal, digits=2), 'cm/s' ,sep= ' ')
sum_deep2<-paste('deep w_vel sum = ', round(sum_deep, digits=2), 'cm/s',sep= ' ')

avg_shal2<-paste('shallow w_vel avg = ', round(avg_shal, digits=3),'cm/s' , sep= " ")
avg_deep2<-paste('deep w_vel avg = ', round(avg_deep, digits=3),'cm/s' , sep= " ")


dev.new()
plot(date, par_deep$depth, pch=19, cex=0.5, xlab= "Time", ylab='Depth (m)')
points(date, par_shal$depth, col='red', pch=19, cex=0.5)

text(x=as.POSIXct("2016-04-05 00:00:00"), y=-500, sum_shal2, cex=0.7)
text(x=as.POSIXct("2016-04-05 00:00:00"), y=-530, sum_deep2, cex=0.7)

text(x=as.POSIXct("2016-04-05 00:00:00"), y=-600, avg_shal2, cex=0.7)
text(x=as.POSIXct("2016-04-05 00:00:00"), y=-630, avg_deep2, cex=0.7)


##########################################################
### load data

## convert tiff file to bathy object to be plotted

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/")

str_name<-"250m_ROMS.tiff"

OR_bathy<-raster(str_name) # open using raster fxn 

OR_bathy2<-as.bathy(OR_bathy) # convert to bathy object

## load in the etopo bedrock data as an .xyz file to plot contours
bathy.dat<-read.table('etopo1_bedrock.xyz', sep='')
names(bathy.dat)<-c('lon', "lat", 'depth')

bathy.dat$depth[bathy.dat$depth>0]<-NA # Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]


MR<-readOGR(dsn=path.expand("/Volumes/TOSHIBA EXT/Desktop_backup/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"), layer= 'MPA_MR_COMP_Boundaries_UTM10N')

class(MR) # SpatialPolygonsDataFrame
crs(MR)

extent(MR)

crsmerc=CRS("+proj=longlat +lat_1=43 +lat_2=48 +lat_0=41 +lon_0=-117 +x_0=700000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ") # this transforms the model 

MR_transformed<-spTransform(MR, CRS=crsmerc)


##########################################################
### plot map to select regions 

blues <- c("lightsteelblue4", "lightsteelblue3","lightsteelblue2", "lightsteelblue1")

plot.bathy(OR_bathy2, deep=-1000, shallow=-10, step=50,image=T, land=T, lwd=0.1, bpal=list(c(0, max(OR_bathy2), 'grey'), c(min(OR_bathy2), 0, blues)), xlab=expression(paste("Longitude ("^o,'W)')), ylab=expression(paste("Latitude ("^o,'N)')))

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(50, 100), labcex=0.7,add=T,col='grey', lwd= 0.5)

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons


plot(OR_bathy2, n = 1, lwd = 0.5, add = TRUE) # plot outline of OR coast		

text(-124.15, 42.836294, "Cape Blanco")
text(-123.81, 44.621745, "Newport")

points(par_deep$lon[1],par_deep$lat[1], col='red', pch=19) # released north of Cape Blanco
points(par_shal$lon[1], par_shal$lat[1], col='green', pch=19) # released north of Cape Blanco

deep_depth<-par_deep$layer[1] # both released in the surface layer 
shal_depth<-par_shal$layer[1] 


##########################################################
### subset out X number of particles are advected deeper than 250 m and X number of particles that are kept at depths of 200 m or less

ID_deep<-unique(data_use$ID[which(data_use$depth<=-250)]) # which particles were advected deeper than 250 m 
length(ID_deep) # 31 particles are advected deeper than 250 m 


# ID_shal3<-unique(data_use$ID[which(data_use$depth>=-200)]) # which particles were kept at 200 m or shallower deeper than 250 m 
# length(ID_shal3) # 750 are kept shallower than 200 m 

ID_shal3<-subset(data_use, !(ID %in% ID_deep))
ID_shal2<-unique(ID_shal3$ID) # 719
# randomly subset out 31 shallow particles 

ID_shal<-sample(ID_shal2, size=length(ID_deep), replace=F)
length(ID_shal)


all_shal<-subset(data_use, ID %in% ID_shal)
dim(all_shal)
# head(all_shal)
# tail(all_shal)


all_deep<-subset(data_use, ID %in% ID_deep)
dim(all_deep)
# head(all_deep)
# tail(all_deep)

length(unique(all_deep$ID))

length(unique(all_shal$ID))

##########################################################
### take the sum and average vertical velocity for ALL of the shallow and deep particles 

sum_deep_all<-sum(all_deep$w_vel, na.rm=T)*100 # -849.014 (units = cm/s)

sum_shal_all<-sum(all_shal$w_vel, na.rm=T)*100 # -163.4449 (units = cm/s)

avg_shal_all<-mean(all_shal$w_vel, na.rm=T)*100 # -0.002465232 (units = cm/s)

avg_deep_all<-mean(all_deep$w_vel, na.rm=T)*100 # -0.008314455 (units = cm/s)




