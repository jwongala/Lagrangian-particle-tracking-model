### Plotting cross shelf transects of data from ROMS 
## Date: 11/30/2020


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
### source functions

setwd('/Volumes/TOSHIBA EXT/Modularized Model_server')

source('makeGrid.R')

##########################################################
### input variables

latd<-371 # low res= 86; high res= 371
lond<-1446 # low res= 355; high res= 1446
depthd<-40 # both 40
depth_exd<-41
lat_exd<-372 # low res= 87; high res= 372
lon_exd<-1447 # low res= 356; high res= 1447
grid_name<-'MarRes_grid_3s.nc' # low res= 'OR_coast_AK_subdomain_grid.nc'; high res= MarRes_grid_3s.nc 
mat_name<-"MarRes_llzh_3.mat" # low res= "OR_coast_AK_grid.mat"; high res= "MarRes_llzh_3.mat"
grid_dir<-setwd("/Volumes/TOSHIBA EXT/Desktop_backup/OR_coast_AK_grid")
	# setwd("/Users/jennifer/Desktop/OR_coast_AK_grid")
	# setwd("/Users/wongalaj/Desktop/Model/OR_coast_AK_grid") # directory of grid files
	# setwd("/home/jwongala/OR_coast_AK/OR_coast_AK_grid") # directory for grid file
	# setwd("/Users/Jenn/Desktop/Model/grid")
	# setwd("/data2/OR_coast_AK_grid") 


##########################################################
### extract ROMS model grid

makeGrid_list<-makeGrid(grid_dir=grid_dir,
						latd=latd,
						lond=lond,
						depthd=depthd,
						lat_exd=lat_exd,
						lon_exd=lon_exd,
						grid_name=grid_name,
						mat_name=mat_name)

	## unlist makeGrid_list 
	lon_rho<-makeGrid_list[[1]] # lon coordinates at rho points
	lat_rho<-makeGrid_list[[2]] # lat coordinates at rho points
	min.lon<-makeGrid_list[[3]] # minimum lon coordinate
	min.lat<-makeGrid_list[[4]] # minimum lat values
	max.lat<-makeGrid_list[[5]]
	z_rho<-makeGrid_list[[6]]
	mask.u<-makeGrid_list[[7]]
	mask.v<-makeGrid_list[[8]]
	mask.rho<-makeGrid_list[[9]]
	boundary<-makeGrid_list[[10]]


##########################################################
### load data

# # setwd('/Volumes/TOSHIBA EXT/ROMS_HiRes/2018')
	
	# # setwd('/Volumes/TOSHIBA EXT/ROMS_HiRes/2016') # high res ROMS directory
	# # setwd('/Volumes/TOSHIBA EXT/April2016_LowRes') # low res ROMS directory	

# files<-list.files('/Volumes/TOSHIBA EXT/ROMS_HiRes/2018', pattern= '*.nc', full.names=T) # list of the ROMS files names to open and loop through

# # setwd('/Volumes/TOSHIBA EXT/April2016_LowRes') # low res ROMS location

# ## convert tiff file to bathy object to be plotted

# vars<-'v'

##########################################################
### extract data for cross shelf transects 

lat_int<-1158 # where do we want to take the transect at? (in terms of grid cell location)
	# high res <-(CB (42.7N) = 312; HB (44.6) = 1158 ) 
	# low res <- (CB (42.7N) = 115 ; HB (44.6N) = 222 )

tmp<-list() # empty vector to fill data with

# tmp[[1]]<-df_tmp

for(i in 1:length(files)){
	
	print(i)
	
	nc.files<-nc_open(files[i]) # open file at i-th iteration
	dat<-ncvar_get(nc.files, varid=vars)
	fillvalue<-ncatt_get(nc.files, vars, "_FillValue") # replace the NAs with 'NA'
	dat[dat==fillvalue]<-NA # put in fill values as NA
	dat2<-0.5*(dat[2:lat_exd,,,] + dat[1:latd,,,]) 
		# v = 0.5*(v_dat[2:lat_exd,,,] + v_dat[1:latd,,,])
		# u = 0.5*(u_dat[,2:lon_exd,,] + u_dat[,1:lond,,])
		# w = 0.5*(w_dat[,,2:depth_exd,] + w_dat[,,1:depthd,])
	
	Tran_dat<-apply(dat2[,lat_int,,], MARGIN=c(1:2), FUN=mean) 
		# cross shelf = choose lat; 
		# along shore = choose lon 
		# 44.6 = HB; grid location is lat_rho[,,]
		# 42.7N = CB; grid location is lat_rho[,312,]

	
	# dim(Tran_dat) # 371  40
	
	tmp[[i]]<-Tran_dat
	
}

##########################################################
### save data 

setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing')
TranDaily_HB18_250m<-tmp
save(TranDaily_HB18_250m, file='vvel_transects_HB2018_dailyavg_250m.RData')


##########################################################
### load in data to plot

## set working directory
setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing/cross_shelf_transects')

load("vvel_transects_HB2016_dailyavg_2km.RData") 

plot_dat<-TranDaily_HB16 # data to plot

##########################################################
### variables needed to plot

xc<-1:dim(z_rho)[1] # number of dimensions in x 
yc<-1:dim(z_rho)[2] # number of dimensions in y

MS<-mapsigma(plot_dat[[1]], sigma=z_rho[,lat_int,], x=xc, resfac=2) # converts the transect matrix from sigma coordinates to depth values and increases the resolution


##########################################################
### plot single transect 

image2D(MS$var, y=MS$depth, x= MS$x, ylim=range(MS$depth), NAcol='black', ylab='depth (m)', xlab= 'x_pos (km)', main= 'EW transect', clab=c("",'m/s'), zlim=c(range(plot_dat, na.rm=T))) # plots the transect data using image2D fxn


##########################################################
### plot transect data using for loop

date<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-02 00:00:00"), by='day') # day_2016 

	## high res
	# seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-02 00:00:00"), by='day') # day_2016 
	# seq(as.POSIXct("2017-03-26 00:00:00"), as.POSIXct("2017-04-28 11:00:00"), by='day') # day_2017
	# seq(as.POSIXct("2018-03-26 00:00:00"), as.POSIXct("2018-05-01 00:00:00"), by='day') # day_2018
	
	## low res
	# seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), by='day') # day_2016 
	# seq(as.POSIXct("2017-03-27 02:00:00"), as.POSIXct("2017-05-02 00:00:00"), by='day') # day_2017
	# seq(as.POSIXct("2018-03-26 02:00:00"), as.POSIXct("2018-05-02 00:00:00"), by='day') # day_2018

ROMS_title<-'2km'

ROMS_title2<-'250m'

year<-2016 

location<-'HB' # 'CB' or 'HB'

text_title<-paste(year, ROMS_title, location, sep='_')
text_title2<-paste(year, ROMS_title2, location, sep='_') 


# set working directory where to save plots at
setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing')


for(i in 1:length(plot_dat)){
	
	print(i)
	
	MS<-mapsigma(plot_dat[[i]], sigma=z_rho[,lat_int,], x=xc, resfac=2) # converts the transect matrix from sigma coordinates to depth values and increases the resolution
	
	## plot data 
	quartz(width=6, height=6)
	dev.copy(jpeg, paste(text_title,'_', date[i], '.jpg', sep=''), height=6, width=6, res=200, units='in')
	
	image2D(MS$var, y=MS$depth, x= MS$x, ylim=range(-7.792118e+02,0), NAcol='black', ylab='depth (m)', xlab= 'x_pos (km)', main= paste(text_title,'_', date[i], sep=''), clab=c("",'m/s'), zlim=c(-0.64, 0.53), xlim=c(40,86)) # plots the transect data using image2D fxn 
	
	dev.off()
	
	## overlapping lon range = -125.1 to -123.9531
	# 2 km ROMS = add in xlim=c(40,86) 

	## zlim range
	# HB = -0.64, 0.53
	# CB = -0.79, 0.58

	
}

# max depth to plot for each transect = -7.792118e+02
##########################################################	
## plot data in panel plot

quartz(width=15, height=15)
dev.copy(jpeg, paste(text_title2,'_', date[i], '_panelplot','.jpg', sep=''), height=15, width=15, res=200, units='in')
par(mfrow=c(6,6)) # , mai=c(0.8, 0.8, 0.6, 0.6)
		
	
for(i in 1:length(plot_dat)){
	
	print(i)
	
	MS<-mapsigma(plot_dat[[i]], sigma=z_rho[,lat_int,], x=xc, resfac=2) # converts the transect matrix from 	sigma coordinates to depth values and increases the resolution
	
	image2D(MS$var, y=MS$depth, x= MS$x, ylim=range(-7.792118e+02,0), NAcol='black', ylab='depth (m)', xlab= 'x_pos', main= date[i], clab=c("",'m/s'), zlim=c(-0.64, 0.53)) # plots the transect data using image2D fxn 
	
	## overlapping lon range = -125.1 to -123.9531
	# 2 km ROMS = add in xlim=c(40,86) 

	## zlim range
	# HB = -0.64, 0.53
	# CB = -0.79, 0.58

	
}


dev.off()
















