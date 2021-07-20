### Plotting cross shelf transects of data from ROMS 
## Date created: 11/30/2020
## last edited by JWA: 2021-07-15


##########################################################
### clear working directory

rm(list=ls())


##########################################################
### load libraries

library(ncdf4)
library(OceanView)
library(R.matlab)
library(marmap)
library(raster)
library(rgdal)
library(ggplot2)
library(zoo)
library(fields)


##########################################################
### load in data to plot
 
## set working directory
setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing/cross_shelf_transects_data')

# load("vvel_transects_446N_2016_dailyavg_250m.RData")
# load("wvel_TranDaily_44.6N_2016_dailyavg_250m.RData")

## 2 km
plot_dat_low<-vvel_TranDaily_44.6N_2018_2km
contour_dat_low<-wvel_TranDaily_44.6N_2018_2km # w-vel data to plot

## 250 m
plot_dat_hi<-TranDaily_446N_2018_250m # v-velocity data to plot
contour_dat_hi<-wvel_TranDaily_44.6N_2018_250m # w-vel data to plot

## set working directory
setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing')

## load in z_rho data 
load("z_rho_hi.RData")
load("z_rho_low.RData")


##########################################################
### dates to plot
date_low<-seq(as.POSIXct("2018-03-26 02:00:00"), as.POSIXct("2018-05-02 00:00:00"), by='day') # day_2018


date_hi<-seq(as.POSIXct("2018-03-26 02:00:00"), as.POSIXct("2018-05-02 00:00:00"), by='day') # day_2018

## low res 
# seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), by='day') # day_2016 
# seq(as.POSIXct("2017-03-27 02:00:00"), as.POSIXct("2017-05-02 00:00:00"), by='day') # day_2017
# seq(as.POSIXct("2018-03-26 02:00:00"), as.POSIXct("2018-05-02 00:00:00"), by='day') # day_2018

## high res
# seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-02 00:00:00"), by='day') # day_2016 
# seq(as.POSIXct("2017-03-26 00:00:00"), as.POSIXct("2017-04-28 11:00:00"), by='day') # day_2017
# seq(as.POSIXct("2018-03-26 00:00:00"), as.POSIXct("2018-05-01 00:00:00"), by='day') # day_2018



##########################################################
## low res (2 km) transect plot

xc<-1:dim(z_rho_low)[1] # number of dimensions in x 
yc<-1:dim(z_rho_low)[2] # number of dimensions in y

# LOW RES ROMS 
## 2016: i= 28 (04/23/2016)
## 2017: i= 15 (04/10/2017)
## 2018: i= 30 (04/24/2018)

i= 30 # day to plot for data

lat_int<-221 # where do we want to take the transect at? (in terms of grid cell location)

## low res ROMS
# 43.3N = [,149,] 
# 44.6N = [,221,]

low_z<-86

contour_dat2<-contour_dat_low[[i]][1:low_z,] # subset so contour data is in the same dimensions as z_rho so they can be converted using the mapsigma fxn
range(contour_dat2, na.rm=T)

contour_dat2<-(contour_dat2)*100 # convert from m/s to cm/s
range(contour_dat2, na.rm=T) 

## 2 km 
MS<-mapsigma(plot_dat_low[[i]], sigma=z_rho_low[,lat_int,], x=xc, resfac=2) # converts the transect matrix from sigma coordinates to depth values and increases the resolution
MS2<-mapsigma(contour_dat2, sigma=z_rho_low[,lat_int,], x=xc, resfac=2) # converts the transect matrix from sigma coordinates to depth values and increases the resolution

ROMS_title<-'2km'

year<-2018

location<-'44.6N' # 433N, 446N

text_title<-paste(year, ROMS_title, location, sep='_')

# set working directory where to save plots at
setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing')

lon1<-(-126.0909) 
lon2<-(-123.9263) # -126.0909 -123.9263 # low res

# 
# zlim1_range16_433<-(-0.8148693) # range(TranDaily_433N_2016_250m, na.rm=T)
# zlim2_range16_433<-(0.4981138)

# zlim1_range16_446<-(-0.5360371) # range(TranDaily_446N_2016_250m, na.rm=T)
# zlim2_range16_446<-(0.4275637)
# 
# zlim1_range17_433<-(-0.704908) # used range(TranDaily_433N_2017_250m, na.rm=T) for the range
# zlim2_range17_433<-(0.623122)

# zlim1_range17_446<-(-0.4972404) # used range(TranDaily_446N_2017_250m, na.rm=T for the range
# zlim2_range17_446<-(0.6548910)
# 
# zlim1_range18_433<-(-0.6983569)
# zlim2_range18_433<-(0.5496013)

zlim1_range18_446<-(-0.5312972) # range(TranDaily_446N_2018_250m, na.rm=T)
zlim2_range18_446<-(0.5882254)


######### 2 km - Begin plotting data 

quartz(width=8, height=8)
dev.copy(jpeg, paste(text_title,'_', date_low[i], "_",'contour.jpg', sep=''), height=8, width=8, res=300, units='in')

image2D(MS$var, y=MS$depth, x= MS$x, ylim=range(-7.792118e+02,0), NAcol='black', ylab='Depth (m)', xlab= expression(paste("Longitude ("^o,'W)')), main= paste(text_title,'_', date_low[i], sep=''), clab=c("",'m/s'), xaxt='n', xlim=c(40,86),zlim=c(zlim1_range18_446, zlim2_range18_446)) # plots the transect data using image2D fxn ,zlim=c(zlim1_range16_433, zlim2_range17_433) 

### vertical (w) velocity contours
contour(MS2$var, y=MS2$depth, x= MS2$x, lty=2, lwd=0.5, add=T, levels = c(-0.001), col='red', labcex=0.5)
contour(MS2$var, y=MS2$depth, x= MS2$x, lty=1, lwd=0.5, add=T, levels = c(0.001), col='blue', labcex=0.5)

contour(MS2$var, y=MS2$depth, x= MS2$x, lty=2, lwd=0.5, add=T, levels = c(-0.002), col='red', labcex=0.5)
contour(MS2$var, y=MS2$depth, x= MS2$x, lty=1, lwd=0.5, add=T, levels = c(0.002), col='blue', labcex=0.5)


xaxislabs<-seq(-125.1, -123.9531, length.out=length(MS$x))

axis(side=1, at= MS$x, labels= signif(xaxislabs, digits=5), tick=F)

dev.off()


##########################################################
## high res (250 m) transect plot

## 250 m 
xc2<-1:dim(z_rho_hi)[1] # number of dimensions in x 
yc2<-1:dim(z_rho_hi)[2] # number of dimensions in y

# plot single examples 
# HIGH RES ROMS
## 2016: i= 28 (04/23/2016)
## 2017: i= 16 (04/10/2017)
## 2018: i= 30 (04/24/2018)

i= 30 # day to plot for data

lat_int<-1158 # where do we want to take the transect at? (in terms of grid cell location)

## high res ROMS 
# 43.3N = [,580,] 
# 44.6N = [,1158,]

ROMS_title2<-'250m'

year<-2018

location<-'44.6N' # 433N, 446N

text_title2<-paste(year, ROMS_title2, location, sep='_')

# set working directory where to save plots at
setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing')

lon1<-(-125.1) 
lon2<-(-123.9531) # -125.1 to -123.9531 # high res

## zlim ranges for transects (v-velcoity)
# 
# zlim1_range16_433<-(-0.8148693) # range(TranDaily_433N_2016_250m, na.rm=T)
# zlim2_range16_433<-(0.4981138)

# zlim1_range16_446<-(-0.5360371) # range(TranDaily_446N_2016_250m, na.rm=T)
# zlim2_range16_446<-(0.4275637)
# 
# zlim1_range17_433<-(-0.704908) # used range(TranDaily_433N_2017_250m, na.rm=T) for the range
# zlim2_range17_433<-(0.623122)
# # 
# zlim1_range17_446<-(-0.4972404) # used range(TranDaily_446N_2017_250m, na.rm=T for the range
# zlim2_range17_446<-(0.6548910)

# zlim1_range18_433<-(-0.6983569)
# zlim2_range18_433<-(0.5496013)

zlim1_range18_446<-(-0.5312972) # range(TranDaily_446N_2018_250m, na.rm=T)
zlim2_range18_446<-(0.5882254)


hi_z<-371
contour_dat2<-contour_dat_hi[[i]][1:hi_z,] # subset so contour data is in the same dimensions as z_rho so they can be converted using the mapsigma fxn
# range(contour_dat2, na.rm=T)

contour_dat2<-(contour_dat2)*100 # convert from m/s to cm/s
# range(contour_dat2, na.rm=T) 

## 250 m

MS<-mapsigma(plot_dat_hi[[i]], sigma=z_rho_hi[,lat_int,], x=xc2, resfac=2) # converts the transect matrix from sigma coordinates to depth values and increases the resolution
MS2<-mapsigma(contour_dat2, sigma=z_rho_hi[,lat_int,], x=xc2, resfac=2) # converts the transect matrix from sigma coordinates to depth values and increases the resolution


######### 250 m - Begin plotting data 

quartz(width=8, height=8)
dev.copy(jpeg, paste(text_title2,'_', date_hi[i], "_",'contour.jpg', sep=''), height=8, width=8, res=300, units='in')

image2D(MS$var, y=MS$depth, x= MS$x, ylim=range(-7.792118e+02,0), NAcol='black', ylab='Depth (m)', xlab= expression(paste("Longitude ("^o,'W)')), main= paste(text_title2,'_', date_hi[i], sep=''), clab=c("",'m/s'), xaxt='n',zlim=c(zlim1_range18_446, zlim2_range18_446) ) # plots the transect data using image2D fxn 

### vertical (w) velocity contours
# contour(MS2$var, y=MS2$depth, x= MS2$x, lty=2, lwd=0.5, add=T, levels = c(-0.05), col='red', labcex=0.5)
# contour(MS2$var, y=MS2$depth, x= MS2$x, lty=1, lwd=0.5, add=T, levels = c(0.05), col='blue', labcex=0.5)

contour(MS2$var, y=MS2$depth, x= MS2$x, lty=2, lwd=0.5, add=T, levels = c(-0.05), col='red', labcex=0.5)
contour(MS2$var, y=MS2$depth, x= MS2$x, lty=1, lwd=0.5, add=T, levels = c(0.05), col='blue', labcex=0.5)

contour(MS2$var, y=MS2$depth, x= MS2$x, lty=2, lwd=0.5, add=T, levels = c(-0.1), col='red', labcex=0.5)
contour(MS2$var, y=MS2$depth, x= MS2$x, lty=1, lwd=0.5, add=T, levels = c(0.1), col='blue', labcex=0.5)


xaxislabs<-seq(-125.1, -123.9531, length.out=length(MS$x))

axis(side=1, at= MS$x, labels= signif(xaxislabs, digits=5), tick=F)

dev.off()
