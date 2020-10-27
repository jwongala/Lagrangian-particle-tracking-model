### Purpose: to visualize the data that will be used to groundtruth the high resolution model
### Jennifer Wong-Ala
### last worked on: 03/30/2020

#########################################################
### clear workspace

rm(list=ls())

#########################################################
### load libraries 

library(marmap) # to plot bathymetry
library(ncdf4) # to load netcdf advection files
# library(data.table) # to use fread
library(mgcv) # to do GAMs
library(sp)
library(raster)
# library(rgeos)
# library(maptools)
library(grDevices)
library(fields)
# library(pracma)
library(OceanView)
library(lattice)
library(maps)
# library(mapdata)
library(fields)
library(sgeostat)
library(rgdal)
library(colorspace)
library(RColorBrewer)
# library(sf)
library(R.matlab)
# library(rasterVis)
# library(abind)
# library(oce)
# library(mapplots)
# library(mapproj)
# library(UScensus2000tract)
library(grid)
# library(gridExtra)

# setwd("/Users/jennifer/Desktop/Modularized Model") 
# source("distance.function.R") # LC function to calculate distances between two points


#########################################################
### Load OOI Temp Data

setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/Point_velocity_OOI/Daily_Averages')
# setwd("/Users/jennifer/Desktop/Groundtruthing_Data/OOI_data_and_comparison_data/Daily_Averaged_Temperature_Data")

### CE01
CE012016_dayavg<-read.csv("Daily_Averages_CE012016.csv")
head(CE012016_dayavg)
dim(CE012016_dayavg)

CE012018_dayavg<-read.csv("Daily_Averages_CE012018.csv")
head(CE012018_dayavg)
dim(CE012018_dayavg)


### CE02
CE022016_dayavg<-read.csv("Daily_Averages_CE022016.csv")
head(CE022016_dayavg)
dim(CE022016_dayavg)

CE022017_dayavg<-read.csv("Daily_Averages_CE022017.csv")
head(CE022017_dayavg)
dim(CE022017_dayavg)

CE022018_dayavg<-read.csv("Daily_Averages_CE022018.csv")
head(CE022018_dayavg)
dim(CE022018_dayavg)



#########################################################
### Load ROMS Temp Data

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/OOI_data_and_comparison_data/Daily_Averaged_Temperature_Data")

### CE01 Comparison

ROMS_CE012016_dayavg<-read.csv("Daily_Averages_ROMS_CE012016.csv")
dim(ROMS_CE012016_dayavg)
head(ROMS_CE012016_dayavg)

ROMS_CE012018_dayavg<-read.csv("Daily_Averages_ROMS_CE012018.csv")
dim(ROMS_CE012018_dayavg)
head(ROMS_CE012018_dayavg)


### CE02 Comparison 

ROMS_CE022016_dayavg<-read.csv("Daily_Averages_ROMS_CE022016.csv")
dim(ROMS_CE022016_dayavg)
head(ROMS_CE022016_dayavg)

ROMS_CE022017_dayavg<-read.csv("Daily_Averages_ROMS_CE022017.csv")
dim(ROMS_CE022017_dayavg)
head(ROMS_CE022017_dayavg)

ROMS_CE022018_dayavg<-read.csv("Daily_Averages_ROMS_CE022018.csv")
dim(ROMS_CE022018_dayavg)
head(ROMS_CE022018_dayavg)

#########################################################
### Load OOI Salinity Data

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/OOI_data_and_comparison_data/Daily_Averaged_Salinity_Data")

### CE01
CE012016_salt_dayavg<-read.csv("Daily_salinty_Averages_CE012016.csv")
head(CE012016_salt_dayavg)
dim(CE012016_salt_dayavg)

CE012018_salt_dayavg<-read.csv("Daily_salinty_Averages_CE012018.csv")
head(CE012018_salt_dayavg)
dim(CE012018_salt_dayavg)


### CE02
CE022016_salt_dayavg<-read.csv("Daily_salinty_Averages_CE022016.csv")
head(CE022016_salt_dayavg)
dim(CE022016_salt_dayavg)

CE022017_salt_dayavg<-read.csv("Daily_salinty_Averages_CE022017.csv")
head(CE022017_salt_dayavg)
dim(CE022017_salt_dayavg)

CE022018_salt_dayavg<-read.csv("Daily_salinty_Averages_CE022018.csv")
head(CE022018_salt_dayavg)
dim(CE022018_salt_dayavg)


#########################################################
### Load ROMS Salinity Data


### CE01
ROMS_CE012016_salt_dayavg<-read.csv("Daily_salinty_Averages_ROMS_CE012016.csv")
head(ROMS_CE012016_salt_dayavg)
dim(ROMS_CE012016_salt_dayavg)

ROMS_CE012018_salt_dayavg<-read.csv("Daily_salinty_Averages_ROMS_CE012018.csv")
head(ROMS_CE012018_salt_dayavg)
dim(ROMS_CE012018_salt_dayavg)


### CE02
ROMS_CE022016_salt_dayavg<-read.csv("Daily_salinty_Averages_ROMS_CE022016.csv")
head(ROMS_CE022016_salt_dayavg)
dim(ROMS_CE022016_salt_dayavg)

ROMS_CE022017_salt_dayavg<-read.csv("Daily_salinty_Averages_ROMS_CE022017.csv")
head(ROMS_CE022017_salt_dayavg)
dim(ROMS_CE022017_salt_dayavg)

ROMS_CE022018_salt_dayavg<-read.csv("Daily_salinty_Averages_ROMS_CE022018.csv")
head(ROMS_CE022018_salt_dayavg)
dim(ROMS_CE022018_salt_dayavg)



##########################################################
### source functions

setwd('/Volumes/TOSHIBA EXT/Modularized Model_server')

source('makeGrid.R')
source('distance.function.R')

##########################################################
### input variables

latd<-371 # low res= 86; high res= 371
lond<-1446 # low res= 355; high res= 1446
depthd<-40 # both 40
depth_exd<-41
lat_exd<-372 # low res= 87; high res= 372
lon_exd<-1447 # low res= 356; high res= 1447
grid_name<-"MarRes_grid_3s.nc" # low res= 'OR_coast_AK_subdomain_grid.nc'; high res= MarRes_grid_3s.nc
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


#########################################################
### Plot OOI vs ROMS data

## Date string to plot data against

day_2018<-seq(as.POSIXct("2018-04-01 00:00:01"), as.POSIXct("2018-04-30 23:00:01"), by="day")
day_2016<-seq(as.POSIXct("2016-04-01 00:00:01"), as.POSIXct("2016-04-30 23:00:01"), by="day")
day_2017<-seq(as.POSIXct("2017-04-01 00:00:01"), as.POSIXct("2017-04-30 23:00:01"), by="day")


#############
#### CE01
plot(day_2016, CE012016_dayavg$temp)

plot(day_2018, CE012018_dayavg$temp)

#### CE02
plot(day_2016, CE022016_dayavg$temp)

plot(day_2017, CE022017_dayavg$temp)

plot(day_2018, CE022018_dayavg$temp)

#### ROMS 2016

## CE01 Comparison
plot(day_2016, ROMS_CE012016_salt_dayavg$salt, type='l', lwd=1.5, xlab='Time', ylab=NA, main="2016 Inshore Oregon Practical Salinity", ylim=c(30.5,33.5))
points(day_2016, CE012016_salt_dayavg$salt, col='purple', type='l', lwd=1.5)
legend(x='topleft', y=NULL, legend=c('ROMS', 'CE01 ISSM'), col=c('black', 'purple'), lty=1)

## CEO2 Comparison
dev.new()
plot(day_2016, ROMS_CE022016_salt_dayavg$salt, , type='l', lwd=1.5, xlab='Time', ylab=NA, main="2016 Midshelf Oregon Practical Salinity", ylim=c(30.5,33.65))
points(day_2016, CE012016_salt_dayavg$salt, col='purple', type='l', lwd=1.5)
legend(x='topleft', y=NULL, legend=c('ROMS', 'CE02 SHSM'), col=c('black', 'purple'), lty=1)

#### ROMS 2017
## CEO2 Comparison
dev.new()
plot(day_2017, ROMS_CE022017_salt_dayavg$salt, , type='l', lwd=1.5, xlab='Time', ylab=NA, main="2017 Midshelf Oregon Practical Salinity", ylim=c(31.5,33.5))
points(day_2017, CE022017_salt_dayavg$salt, col='purple', type='l', lwd=1.5)
legend(x='topleft', y=NULL, legend=c('ROMS', 'CE02 SHSM'), col=c('black', 'purple'), lty=1)

#### ROMS 2018

## CE01
dev.new()
plot(day_2018, ROMS_CE012018_salt_dayavg$salt, type='l', lwd=1.5, xlab='Time', ylab=NA, main="2018 Inshore Oregon Practical Salinity", ylim=c(30.5,33.6))
points(day_2018, CE012018_salt_dayavg$salt, col='purple', type='l', lwd=1.5)
legend(x='topleft', y=NULL, legend=c('ROMS', 'CE01 ISSM'), col=c('black', 'purple'), lty=1)


## CE02
dev.new()
plot(day_2018, ROMS_CE022018_salt_dayavg$salt, type='l', lwd=1.5, xlab='Time', ylab=NA, main="2018 Midshelf Oregon Practical Salinity", ylim=c(31.8,33.56))
points(day_2018, CE022018_salt_dayavg$salt, col='purple', type='l', lwd=1.5)
legend(x='topleft', y=NULL, legend=c('ROMS', 'CE02 SHSM'), col=c('black', 'purple'), lty=1)



#########################################################
### Load MR polygons

## set working directory

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files")
MR<-shapefile('MPA_MR_COMP_Boundaries_UTM10N.shp')
class(MR) # SpatialPolygonsDataFrame
crs(MR)
extent(MR)
crsmerc=CRS("+proj=longlat +lat_1=43 +lat_2=48 +lat_0=41 +lon_0=-117 +x_0=700000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ") # this transforms the model 
MR_transformed<-spTransform(MR, CRS=crsmerc)

## load in the etopo bedrock data as an .xyz file to plot contours
setwd("/Volumes/TOSHIBA EXT/Desktop_backup/")

bathy.dat<-read.table('etopo1_bedrock.xyz', sep='')
names(bathy.dat)<-c('lon', "lat", 'depth')

bathy.dat$depth[bathy.dat$depth>0]<-NA # Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]



#########################################################
### Plot surface currents

# NOTE: make sure to change name for ROMS vs OOI data 

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/MK_Radar_Data_lower_res/April_2018")

	mk_files<-list.files("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/MK_Radar_Data_lower_res/April_2018", pattern= '*.txt', full.names=T) 

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/ROMS_surface_curnts_April2018")

	ROMS_files<-list.files("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/ROMS_surface_curnts_April2018", pattern= '*.RData', full.names=T)


OR_bathy<-getNOAA.bathy(lon1= -124.96, lon2= -123.93, lat1= 45.28, lat2= 41.85, resolution=1) # get bathymetry from marmap (NOAA data)
# -125.1000 -123.9531 42.00000 45.24775


blues<-c( "lightsteelblue3", "lightsteelblue2", "lightsteelblue1") # "lightsteelblue4",  ,

 
for(i in 1:length(mk_files)){
	
	
	setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/")
	
	print(i)
	
	ROMS_dat<-load(ROMS_files[i])
	
	# mk_dat<-read.csv(mk_files[i])
	mk_dat<-read.table(mk_files[i], header=T, fill=T)
	# mk_dat<-mk_dat[,-6]
	# dim(mk_dat)
	colnames(mk_dat)<-c('lon', 'lat', 'prcnt','u','v')
	
	mk_dat_mat<-floor(sqrt(nrow(mk_dat)))
	
	# convert u, v, lat and lon
	
	u_rom<-daily_avg[[1]] # m/s
	# dim(u_rom)
	
	u_mk<-matrix(mk_dat$u, ncol= mk_dat_mat, nrow= mk_dat_mat)  # cm/s
	# dim(u_mk)
	
	v_rom<-daily_avg[[2]] # m/s
	# dim(v_rom)
	
	v_mk<-matrix(mk_dat$v, ncol= mk_dat_mat, nrow= mk_dat_mat) # cm/s
	# dim(v_mk)
	
	lat_mk<-matrix(mk_dat$lat, ncol= mk_dat_mat, nrow= mk_dat_mat)
	# dim(lat_mk)
	lat_rom<-lat_rho[,,1]
	
	lon_mk<-matrix(mk_dat$lon, ncol= mk_dat_mat, nrow= mk_dat_mat)
	lon_rom<-lon_rho[,,1]

	ROMS_mat<-NULL
	
	for(k in 1:dim(mk_dat)[1]){
	
		print(k)
		
		dist<-distance.function(mk_dat$lat[k], mk_dat$lon[k],lat_rom, lon_rom)
		tmp<-1*(dist==min(dist))
	
		tmp2<-which(tmp==1, arr.ind=TRUE)
	
		ROMS_mat<-rbind(ROMS_mat,tmp2)
	
	}

	row_rom<-ROMS_mat[,1]
	col_rom<-ROMS_mat[,2]
		
	lat_rom2<-lat_rom[1,ROMS_mat[,2]]
	length(lat_rom2)
		
	lon_rom2<-lon_rom[ROMS_mat[,1],1]
	length(lon_rom2)
		
	# grid_rom<-cbind(lat_rom2, lon_rom2)
	
	rom_mat_dat<-floor(sqrt(length(lat_rom2)))
	
	u_rom2<-u_rom[ROMS_mat[,1],1]
		
	v_rom2<-v_rom[1,ROMS_mat[,2]]
	
	lat_rom3<-matrix(lat_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat) # , ncol= 18, nrow= 18)
	
	lon_rom3<-matrix(lon_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat) # , ncol= 18, nrow= 18)
	
	u_rom3<-matrix(u_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat) #, ncol= 18, nrow= 18)

	v_rom3<-matrix(v_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat) # , ncol= 18, nrow= 18)
	
	v_mk2<-v_mk/100
	
	u_mk2<-u_mk/100
	
	curnt1<-sqrt(u_rom3 ^2+ v_rom3 ^2)
	curnt2<-sqrt(u_mk2 ^2+ v_mk2 ^2)
	
	ROMS_max<-round(max(curnt1, na.rm=T), digits=1)
	MK_max<-round(max(curnt2, na.rm=T), digits=1)
	
	# dim(curnt2)
	# lat<-lat_rho[,,1]
		
	setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data')
	
	# setwd("/Users/jennifer/Desktop/Groundtruthing_Data/HFRadar_and_comparison_data/reversal_example_figure")
	
	png(paste("Surface_Current_Comparison_2018_", i, ".png", sep= "" ), width=687 , height= 1790, res= 150)
	
	# png(paste("MK_radar_April2018_", i, ".png", sep= "" ), 	width=830 , height= 1885, res= 150)
	# png(paste("ROMS_surfcurnt_April2018_", i, ".png", sep= "" ), 	width=830 , height= 1885, res= 150)
	
	par(mai=c(1.1,1.1,1,0.5))
	
	# plot current vectors
	
	plot(OR_bathy, image = T, land = T, axes = T, lwd=0,  deep=0,shallow=0, step=0, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)), main='April 07, 2018',ylab=expression(paste("Latitude ("^o,'N)')), xlab='', cex.main=2, cex.axis=1.5, cex.lab=2) 

	# plot outline of continent
	plot(OR_bathy, n = 1, lwd = 0.5, add = TRUE) # plot outline of OR coast

	contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=0.5, levels=c(-200, -500, -1000),col= "gray40",add=T, lty=2)
	
	quiver2D(u= u_mk2, v= v_mk2, x= lon_mk, y= lat_mk, type= "simple",  col= 'black',  lwd=2, yaxt='n', xaxt='n', add=T, arr.max=0.2) # by=70 for only ROMS data scale=1, arr.max = 0.2, arr.min=0.2,
	 
	quiver2D(u= u_rom3, v= v_rom3, x= lon_rom3, y= lat_rom3, type= "simple",  col= '#DC3220', lwd=2, yaxt='n', xaxt='n',add=T, arr.min=0.15) # by=70 for only ROMS data 
	
	lines(MR_transformed, col= 'black', lwd= 2, lty= 1) # plot MR polygons
	
	# text(x=-124.17, y=42.7, '____', cex=1.5, lwd=1)
	# text(x=-124.15, y=42.65, paste(MK_max, "m/s - HFR", sep= " " ), cex=1, font=2)
	
	# text(x=-124.17, y=42.5, '____', cex=1.5, lwd=1,col='#DC3220')
	# text(x=-124.15, y=42.45, paste(ROMS_max, "m/s - ROMS", sep= " " ), cex=1, col='#DC3220', font=2)

	text(x=-124.25, y=42.845, 'Cape Blanco', cex=1.5)
	
	# text(x=-124.1, y=42.19, 'HFR', cex=1.6, col='black',)
	# text(x=-124.1, y=42.1, 'ROMS', cex=1.6, col='#DC3220')
		
	dev.off()
	
}

### Combine reversal examples to make 1 image
library(png)

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/reversal_example_figure")

pic1<-readPNG("Surface_Current_Comparison_2016_14_USE.png")
pic2<-readPNG("Surface_Current_Comparison_2017_07_use.png")
pic3<-readPNG("Surface_Current_Comparison_2018_07_use.png")


grid.arrange(rasterGrob(pic1),rasterGrob(pic2), rasterGrob(pic3),ncol=3, bottom=textGrob(expression(paste("Longitude ("^o,'W)')), gp=gpar(fontsize=14)), left=textGrob(expression(paste("Latitude ("^o,'N)')), gp=gpar(fontsize=14), rot=90), top=textGrob(c('April 14, 2016', 'April 07, 2017', 'April 07, 2018'), x=c(0.17, 0.51, 0.84), vjust=1.2, gp=gpar(fontsize=16)))


# # par(mfrow=c(1,3))
# layout(mat = matrix(c(1,2,3),nrow = 1, ncol = 3))  

# # plot(0:3, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# plot(NA,xlim=0:1,ylim=0:1,bty="n",axes=0,xaxs = 'i',yaxs='i')

# rasterImage(pic1, 0, 0, 1, 1)
# rasterImage(pic2,  1, 1, 2, 2)
# rasterImage(pic3,  2, 2, 1, 1)


# plot(0:2, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# rasterImage(readPNG(source="my_viz1.png"), 0, 1, 1, 2)
# rasterImage(readPNG(source="my_viz2.png"), 1, 1, 2, 2)
# rasterImage(readPNG(source="my_viz3.png"), 0, 0, 1, 1)
# rasterImage(readPNG(source="my_viz4.png"), 1, 0, 2, 1)



#########################################################
### Monthly mean for all surface currents (ROMS and HFR)

## Empty list
v_rom_list<-list()
u_rom_list<-list()

v_mk_list<-list()
u_mk_list<-list()


## Empty arrays
# v_rom_all<-array(numeric(), c(371, 1446,0))
# u_rom_all<-array(numeric(), c(371,1446,0))


## Directories for where data is
# setwd("/Users/jennifer/Desktop/Groundtruthing_Data/HFRadar_and_comparison_data/MK_Radar_Data_lower_res/April_2016")

	# mk_files<-list.files("/Users/jennifer/Desktop/Groundtruthing_Data/HFRadar_and_comparison_data/MK_Radar_Data_lower_res/April_2016", pattern= '*.txt', full.names=T) 

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/ROMS_surface_curnts_April2016")

	ROMS_files<-list.files("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/ROMS_surface_curnts_April2016", pattern= '*.RData', full.names=T)



#### to plot OR bathy
OR_bathy<-getNOAA.bathy(lon1= -125.3, lon2= -123.93, lat1= 45.28, lat2= 41.85, resolution=1) # get bathymetry from marmap (NOAA data)
# -125.1000 -123.9531 42.00000 45.24775

# color for bathy
blues<-c( "lightsteelblue3", "lightsteelblue2", "lightsteelblue1") # "lightsteelblue4",  ,



### Download monthly averaged HFR data

hf_month<-nc_open('USWC-month-LTA-6km_april2016.nc')
dim(hf_month)
length(hf_month)

time_avg<-ncvar_get(hf_month, var='time')

time_avg2<-as.POSIXct(time_avg, origin='1970-01-01 00:00:00', tz='GMT')
length(time_avg2)

lat_mk<-ncvar_get(hf_month, var='lat')
dim(lat_mk)

lon_mk<-ncvar_get(hf_month, var='lon')
dim(lon_mk)

v_avg_mk<-ncvar_get(hf_month, var='v_mean') # units m/s
dim(v_avg_mk) 

u_avg_mk<-ncvar_get(hf_month, var='u_mean') # units m/s
dim(u_avg_mk)

### subset out for v April 2016, 2017, 2018
v_avg_mk_2016<-v_avg_mk[,,1] # units m/s
dim(v_avg_mk_2016)

v_avg_mk_2017<-v_avg_mk[,,13] # units m/s
dim(v_avg_mk_2017)
 
v_avg_mk_2018<-v_avg_mk[,,25] # units m/s
dim(v_avg_mk_2018)


### subset out for u April 2016, 2017, 2018
u_avg_mk_2016<-u_avg_mk[,,1] # units m/s
dim(u_avg_mk_2016)

u_avg_mk_2017<-u_avg_mk[,,13] # units m/s
dim(u_avg_mk_2017)
 
u_avg_mk_2018<-u_avg_mk[,,25] # units m/s
dim(u_avg_mk_2018)



### Empty arrays
v_rom_all<-array(numeric(), c(371, 1446,0))
u_rom_all<-array(numeric(), c(371,1446,0))


for(i in 1:length(ROMS_files)){
	
	setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/HFRadar_and_comparison_data/")
	
	print(i)
	
	ROMS_dat<-load(ROMS_files[i])
	
	# lat and lon coordinates for roms and mk data
	
	lat_rom<-lat_rho[,,1]
	dim(lat_rom)

	lon_rom<-lon_rho[,,1]
	dim(lon_rom)
	
	# u and v data from roms and mk data 	
	u_rom<-daily_avg[[1]] 
	# dim(u_rom)

	v_rom<-daily_avg[[2]]
	# dim(v_rom)
		
	# bind all data so that it is in one array to take a monthly average from
	v_rom_all<-abind(v_rom_all, v_rom,along=3)
	
	u_rom_all<-abind(u_rom_all, u_rom, along=3)

}
		

### Take mean of u, v, and convert ROMS to cm

dim(v_rom_all)
dim(u_rom_all)

v_rom_2016<-apply(v_rom_all, c(1,2), FUN=mean, na.rm=T)
dim(v_rom_2016) # 371 1446

u_rom_2016<-apply(u_rom_all, c(1,2), FUN=mean, na.rm=T)
dim(u_rom_2016) # 15 15


save(u_rom_2017, file="u_monthly_average_2017.RData")

save(v_rom_2017, file="v_monthly_average_2017.RData")

	
### Create grid of lat and lon mk data 
mk_grid3<-expand.grid(lon_mk, lat_mk)
colnames(mk_grid3)<-c('lon', 'lat')
dim(mk_grid3)
# -125.1000 -123.9531 42.00000 45.24775

mk_grid2<-subset(mk_grid3, lat>=42.00000 & lat<=45.24775)
dim(mk_grid2)

mk_grid<-subset(mk_grid2, lon>=-125.1000 & lon<=-123.9531)
dim(mk_grid)


### Subset out lat and lon locations that are closest to mk coordinates
ROMS_mat<-NULL

for(k in 1:nrow(mk_grid)){
		
	print(k)
		
	dist<-distance.function(mk_grid$lat[k], mk_grid$lon[k], lat_rom, lon_rom)
	tmp<-1*(dist==min(dist))
		
	tmp2<-which(tmp==1, arr.ind=TRUE)
		
	ROMS_mat<-rbind(ROMS_mat,tmp2)
		
}

row_rom<-ROMS_mat[,1]
col_rom<-ROMS_mat[,2]
	
lat_rom2<-lat_rom[1,ROMS_mat[,2]]
length(lat_rom2)
		
lon_rom2<-lon_rom[ROMS_mat[,1],1]
length(lon_rom2)
		
rom_mat_dat<-floor(sqrt(length(lat_rom2)))

	
u_rom2<-u_rom_2016[ROMS_mat[,1],1]
		
v_rom2<-v_rom_2016[1,ROMS_mat[,2]]
	
	
lat_rom3<-matrix(lat_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat)
	
lon_rom3<-matrix(lon_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat)
	
u_rom_2016<-matrix(u_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat)

v_rom_2016<-matrix(v_rom2, ncol= rom_mat_dat, nrow= rom_mat_dat)
	
dim(v_rom_2018)
	
# image.plot((u_rom_2018^2)+(v_rom_2018^2))

### Data to plot
u_rom_20162<-u_rom_2016

u_rom_20172<-u_rom_2017

u_rom_20182<-u_rom_2018

v_rom_20162<-v_rom_2016

v_rom_20172<-v_rom_2017

v_rom_20182<-v_rom_2018

### subset the NaN for NAs
u_rom_20162[u_rom_20162==NaN]<-NA
head(u_rom_20162)

u_rom_20172[u_rom_20172==NaN]<-NA
head(u_rom_20172)

u_rom_20182[u_rom_20182==NaN]<-NA
head(u_rom_20182)

v_rom_20162[v_rom_20162==NaN]<-NA
head(v_rom_20162)

v_rom_20172[v_rom_20172==NaN]<-NA
head(v_rom_20172)

v_rom_20182[v_rom_20182==NaN]<-NA
head(v_rom_20182)


### Plot monthly mean for both ROMS and HFR

# # setwd('/Users/jennifer/Desktop/Groundtruthing_Data')

# png("Monthly_mean_Surface_Current_Comparison_test.png", 	width=825 , height= 1840, res= 150)

# par(mai=c(1.1,1.1,1.2,0.4))

# plot(OR_bathy, image = T, land = T, axes = T, lwd=1, deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)), cex.axis= 1.5, cex.lab=1.5,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main="Monthly Averaged Surface Currents for April 2016", cex.main=1.5) 

# # plot(OR_bathy, image = T, land = T, axes = T, lwd=1) # , deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)), cex.axis= 1.5, cex.lab=1.5,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main="Monthly Averaged Surface Currents for April 2016", cex.main=1.5) 

# # points(lon_rom3, lat_rom3, col='darkgreen', cex=2)	
# # v_avg_mk_2016
# quiver2D(u= u_avg_mk_2017, v= v_avg_mk_2017, x= lon_mk, y= lat_mk, type= "triangle",  col= 'black', scale=1,  arr.max = 0.2, lwd=3, yaxt='n', xaxt='n', add=T) # by=70 for only ROMS data 

# # points(lon_mk, lat_mk, col='blue4', cex=2)
# # map(database='usa', add=T)

# quiver2D(u= u_rom_20172, v= v_rom_20172, x= lon_rom3, y= lat_rom3, type= "triangle",  col= '#DC3220', scale=1 , arr.max = 0.2, lwd=3, yaxt='n', xaxt='n', add=T) # by=70 for only ROMS data 


# lines(MR_transformed, col= 'black', lwd= 2, lty= 1) # plot MR polygons
	
# text(x=-124.33, y=42.4, '_', cex=4.7, lwd=1)
# text(x=-124.1, y=42.4, '25 cm/s', cex=1.5)
	
# text(x=-124.3, y=42.845, 'Cape Blanco', cex=1.5)
	
# text(x=-124.1, y=42.19, 'HFR', cex=1.5, col='#005AB5')
# text(x=-124.1, y=42.1, 'ROMS', cex=1.5, col='#DC3220')

# # dev.off()

setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data')

bathy.dat<-read.table('etopo1_bedrock.xyz',sep='')
  names(bathy.dat)<-c('lon','lat','depth')
  bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
  head(bathy.dat)
  bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]

# curnt16<-sqrt((u_rom_2016^2) +(v_rom_2016^2))
# range(curnt16, na.rm=T) # 0.01309441 0.29914588

# curnt17<-sqrt((u_rom_2017^2) +(v_rom_2017^2))
# range(curnt17, na.rm=T) # 0.01974926 0.52045828

# curnt18<-sqrt((u_rom_2018^2) +(v_rom_2018^2))
# range(curnt18, na.rm=T) # 0.04879623 0.26061137

mk_curnt16<-sqrt((u_avg_mk_2016^2) +(v_avg_mk_2016^2))
range(mk_curnt16, na.rm=T) # 0.000000 0.510392

mk_curnt17<-sqrt((u_avg_mk_2017^2) +(v_avg_mk_2017^2))
range(mk_curnt17, na.rm=T) # 0.0000000 0.4159327

mk_curnt18<-sqrt((u_avg_mk_2018^2) +(v_avg_mk_2018^2))
range(mk_curnt18, na.rm=T) # 0.0000000 0.5021952





blues<-c("lightsteelblue3", "lightsteelblue2", "lightsteelblue1") # "lightsteelblue4",  , 

setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data')

png("Monthly_mean_Surface_Current_Comparison_test_final3.png", width=2110 , height= 1500, res= 150)
par(mfrow=c(1,3), mai=c(0.8,0.8,0.3,0.3))


### 2016
plot(OR_bathy, image = T, land = F, axes = T, lwd=0,  deep=0,shallow=0, step=0, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)), main='2016',ylab=expression(paste("Latitude ("^o,'N)')), xlab='', cex.main=2, cex.axis=1.5, cex.lab=2) 

quiver2D(u_rom_2016, v_rom_2016, lon_rom3, lat_rom3, col='#DC3220', add=T)
quiver2D(u_avg_mk_2016, v_avg_mk_2016, lon_mk, lat_mk,  col='black', add=T)


map(database='usa', col='gray', add=T, fill=T)

lines(MR_transformed, col= 'black', lwd= 1.5, lty= 1) # plot MR polygons

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=2, levels=-200,col= "gray40",add=T, lty=2)
# contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=1, levels=-100,col= "gray40",add=T)
	
text(x=-124.27, y=42.4, '__', cex=1.2, lwd=2, font=2)
text(x=-124.1, y=42.4, '0.5 m/s', cex=1.5)
	
text(x=-124.25, y=42.845, 'Cape Blanco', cex=1.8)
	
text(x=-124.1, y=42.19, 'HFR', cex=2, col='black')
text(x=-124.1, y=42.1, 'ROMS', cex=2, col='#DC3220')


### 2017
plot(OR_bathy, image = T, land = F, axes = T, lwd=0,  deep=0,shallow=0, step=0, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)),main='2017',xlab=expression(paste("Longitude ("^o,'W)')), ylab='',cex.main=2, cex.axis=1.5, cex.lab=2) 


quiver2D(u_rom_2017, v_rom_2017, lon_rom3, lat_rom3, col='#DC3220', add=T)
quiver2D(u_avg_mk_2017, v_avg_mk_2017, lon_mk, lat_mk,col='black', add=T)

map(database='usa', add=T, col='gray', fill=T)

lines(MR_transformed, col= 'black', lwd= 1.5, lty= 1) # plot MR polygons
text(x=-124.25, y=42.845, 'Cape Blanco', cex=1.8)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=2, levels=-200,col= "gray40",add=T, lty=2)

### 2018
plot(OR_bathy, image = T, land = F, axes = T, lwd=0,  deep=0,shallow=0, step=0, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)), main='2018', xlab='', ylab='',cex.main=2, cex.axis=1.5, cex.lab=2) 

quiver2D(u_rom_2018, v_rom_2018, lon_rom3, lat_rom3, col='#DC3220', add=T)
quiver2D(u_avg_mk_2018, v_avg_mk_2018, lon_mk, lat_mk, add=T, col='black')

map(database='usa', add=T, col='gray', fill=T)

lines(MR_transformed, col= 'black', lwd= 1.5, lty= 1) # plot MR polygons
text(x=-124.25, y=42.845, 'Cape Blanco', cex=1.8)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat, lwd=2, levels=-200,col= "gray40",add=T, lty=2)

dev.off()



#########################################################
### Download point velocity data 

setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/Point_velocity_OOI/OOI_pt_velocity_original")

ce01_2016pt<-read.csv("CE02SHSM_point_vel_2016.csv")
head(ce01_2016pt)

colnames(ce01_2016pt)<-c('lat', 'lon', 'time', 'w_vel', 'v_vel', 'u_vel', 'pitch', 'pressure', 'roll', 'temp', 'heading', 'deployment')
# rtime<-load('ROMS_time_2017.RData')
time_ceo12016<-as.POSIXct(ce01_2016pt$time, origin='1900-01-01 00:00:00', tz='GMT')

ce01_2016pt<-cbind(ce01_2016pt, time_ceo12016)
head(ce01_2016pt)

# plot(x=time_ceo12016 , y= ce01_2016pt$u_vel, typ='l', xlab='time', ylab='u velocity (m/s)')

# plot(x=time_ceo12016 , y= ce01_2016pt$v_vel, typ='l', xlab='time', ylab='v velocity (m/s)')

# class(ce01_2016pt)


### Taking Daily Averages
# ROMS_vvel2 ROMS_uvel2

dayavg_ce0116<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-01')
dim(dayavg_ce0116)
daymean_ce0116<-mean(dayavg_ce0116$u_vel, na.rm=T)

dayavg_ce01162<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-02')
dim(dayavg_ce01162)
daymean_ce01162<-mean(dayavg_ce01162$u_vel, na.rm=T)

dayavg_ce01163<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-03')
dim(dayavg_ce01163)
daymean_ce01163<-mean(dayavg_ce01163$u_vel, na.rm=T)

dayavg_ce01164<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-04')
dim(dayavg_ce01163)
daymean_ce01164<-mean(dayavg_ce01164$u_vel, na.rm=T)

dayavg_ce01165<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-05')
dim(dayavg_ce01165)
daymean_ce01165<-mean(dayavg_ce01165$u_vel, na.rm=T)

dayavg_ce01166<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-06')
dim(dayavg_ce01166)
daymean_ce01166<-mean(dayavg_ce01166$u_vel, na.rm=T)

dayavg_ce01167<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-07')
dim(dayavg_ce01167)
daymean_ce01167<-mean(dayavg_ce01167$u_vel, na.rm=T)


dayavg_ce01168<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-08')
dim(dayavg_ce01168)
daymean_ce01168<-mean(dayavg_ce01168$u_vel, na.rm=T)


dayavg_ce01169<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-09')
dim(dayavg_ce01169)
daymean_ce01169<-mean(dayavg_ce01169$u_vel, na.rm=T)


dayavg_ce011610<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-10')
dim(dayavg_ce011610)
daymean_ce011610<-mean(dayavg_ce011610$u_vel, na.rm=T)


dayavg_ce011611<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-11')
dim(dayavg_ce011611)
daymean_ce011611<-mean(dayavg_ce011611$u_vel, na.rm=T)
12

dayavg_ce011612<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-12')
dim(dayavg_ce011612)
daymean_ce011612<-mean(dayavg_ce011612$u_vel, na.rm=T)


dayavg_ce011613<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-13')
dim(dayavg_ce011613)
daymean_ce011613<-mean(dayavg_ce011613$u_vel, na.rm=T)


dayavg_ce011614<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-14')
dim(dayavg_ce011614)
daymean_ce011614<-mean(dayavg_ce011614$u_vel, na.rm=T)


dayavg_ce011615<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-15')
dim(dayavg_ce011615)
daymean_ce011615<-mean(dayavg_ce011615$u_vel, na.rm=T)


dayavg_ce011616<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-16')
dim(dayavg_ce011616)
daymean_ce011616<-mean(dayavg_ce011616$u_vel, na.rm=T)


dayavg_ce011617<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-17')
dim(dayavg_ce011617)
daymean_ce011617<-mean(dayavg_ce011617$u_vel, na.rm=T)


dayavg_ce011618<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-18')
dim(dayavg_ce011618)
daymean_ce011618<-mean(dayavg_ce011618$u_vel, na.rm=T)


dayavg_ce011619<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-19')
dim(dayavg_ce011619)
daymean_ce011619<-mean(dayavg_ce011619$u_vel, na.rm=T)


dayavg_ce011620<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-20')
dim(dayavg_ce011620)
daymean_ce011620<-mean(dayavg_ce011620$u_vel, na.rm=T)


dayavg_ce011621<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-21')
dim(dayavg_ce011621)
daymean_ce011621<-mean(dayavg_ce011621$u_vel, na.rm=T)


dayavg_ce011622<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-22')
dim(dayavg_ce011622)
daymean_ce011622<-mean(dayavg_ce011622$u_vel, na.rm=T)


dayavg_ce011623<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-23')
dim(dayavg_ce011623)
daymean_ce011623<-mean(dayavg_ce011623$u_vel, na.rm=T)


dayavg_ce011624<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-24')
dim(dayavg_ce011624)
daymean_ce011624<-mean(dayavg_ce011624$u_vel, na.rm=T)


dayavg_ce011625<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-25')
dim(dayavg_ce011625)
daymean_ce011625<-mean(dayavg_ce011625$u_vel, na.rm=T)


dayavg_ce011626<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-26')
dim(dayavg_ce011626)
daymean_ce011626<-mean(dayavg_ce011626$u_vel, na.rm=T)


dayavg_ce011627<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-27')
dim(dayavg_ce011627)
daymean_ce011627<-mean(dayavg_ce011627$u_vel, na.rm=T)


dayavg_ce011628<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-28')
dim(dayavg_ce011628)
daymean_ce011628<-mean(dayavg_ce011628$u_vel, na.rm=T)


dayavg_ce011629<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-29')
dim(dayavg_ce011629)
daymean_ce011629<-mean(dayavg_ce011629$u_vel, na.rm=T)


dayavg_ce011630<-subset(ce01_2016pt, format(time_ceo12016, '%Y-%m-%d')=='2016-04-30')
dim(dayavg_ce011630)
daymean_ce011630<-mean(dayavg_ce011630$u_vel, na.rm=T)

totdayavg_ce0218<-c(daymean_ce0116,daymean_ce01162,daymean_ce01163,daymean_ce01164,daymean_ce01165,daymean_ce01166,daymean_ce01167,daymean_ce01168,daymean_ce01169,daymean_ce011610,daymean_ce011611,daymean_ce011612,daymean_ce011613,daymean_ce011614,daymean_ce011615,daymean_ce011616,daymean_ce011617,daymean_ce011618,daymean_ce011619,daymean_ce011620,daymean_ce011621,daymean_ce011622,daymean_ce011623,daymean_ce011624,daymean_ce011625,daymean_ce011626,daymean_ce011627,daymean_ce011628,daymean_ce011629,daymean_ce011630)

day_2018<-seq(as.POSIXct("2018-04-01 00:00:01"), as.POSIXct("2018-04-30 23:00:01"), by="day")
day_2016<-seq(as.POSIXct("2016-04-01 00:00:01"), as.POSIXct("2016-04-30 23:00:01"), by="day")
day_2017<-seq(as.POSIXct("2017-04-01 00:00:01"), as.POSIXct("2017-04-30 23:00:01"), by="day")

# CE022018_dayavg<-data.frame(cbind(totdayavg_ce0218, day_2016))
# colnames(CE022018_dayavg)<-c('v_velocity', 'date')

CE022018_dayavg2<-cbind(totdayavg_ce0218, CE022018_dayavg)
colnames(CE022018_dayavg2)<-c('u_velocity','v_velocity', 'date')

setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/Point_velocity_OOI')

write.csv(CE022018_dayavg2, 'Daily_velocity_Averages_CE022016.csv')



#########################################################
### Plot point velocity data 
setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/Point_velocity_OOI/Daily_Averages')

# setwd('/Users/jennifer/Desktop/Groundtruthing_Data/Point_velocity_OOI/Daily_Averages')

day_2018<-seq(as.POSIXct("2018-04-01 00:00:01"), as.POSIXct("2018-04-30 23:00:01"), by="day")
day_2016<-seq(as.POSIXct("2016-04-01 00:00:01"), as.POSIXct("2016-04-30 23:00:01"), by="day")
day_2017<-seq(as.POSIXct("2017-04-01 00:00:01"), as.POSIXct("2017-04-30 23:00:01"), by="day")

### CE01 
ce01_2016<-read.csv("Daily_velocity_Averages_CE012016.csv")
dim(ce01_2016)

ce01_2017<-read.csv("Daily_velocity_Averages_CE012017.csv")
dim(ce01_2017)

ce01_2018<-read.csv("Daily_velocity_Averages_CE012018.csv")
dim(ce01_2018)

### ROMS - CE01
ROMS_ce01_2016<-read.csv("Daily_velocity_Averages_CE012016_ROMS.csv")
dim(ROMS_ce01_2016)

ROMS_ce01_2017<-read.csv("Daily_velocity_Averages_CE012017_ROMS.csv")
dim(ROMS_ce01_2017)

ROMS_ce01_2018<-read.csv("Daily_velocity_Averages_CE012018_ROMS.csv")
dim(ROMS_ce01_2018)


### CE02 
ce02_2016<-read.csv("Daily_velocity_Averages_CE022016.csv")
dim(ce02_2016)

ce02_2017<-read.csv("Daily_velocity_Averages_CE022017.csv")
dim(ce02_2017)

### ROMS - CE02

ROMS_ce02_2016<-read.csv("Daily_velocity_Averages_CE022016_ROMS.csv")
dim(ROMS_ce02_2016)

ROMS_ce02_2017<-read.csv("Daily_velocity_Averages_CE022017_ROMS.csv")
dim(ROMS_ce02_2017)


### Plot data

day_april<-seq(as.POSIXct("2016-04-01 00:00:01"), as.POSIXct("2016-04-30 23:00:01"), by="day")

#### Plot u on same fig (diff cols based on year and ROMS= dashed, OOI= solid)
# 2016= red, 2017= blue, 2018= purple


### Plot of all data on one figure
plot(day_april, ce01_2016$u_vel, type='l', lty=1, col='red', ylim=c(-0.25, 0.19), xlab='time', ylab='u velocity (m/s)') # , ylim=c()
lines(day_april, ce02_2016$u_vel, col='red')

lines(day_april, ce01_2017$u_vel, type='l', col='blue')
lines(day_april, ce02_2017$u_vel, col='blue')

lines(day_april, ce01_2018$u_vel, type='l', col='purple')

# ROMS
lines(day_april, ROMS_ce01_2016$u_vel, type='l', col='red', lty=2)
lines(day_april, ROMS_ce02_2016$u_vel, type='l', col='red', lty=2)


lines(day_april, ROMS_ce01_2017$u_vel, type='l', col='blue', lty=2)
lines(day_april, ROMS_ce02_2017$u_vel, type='l', col='blue', lty=2)

lines(day_april, ROMS_ce01_2018$u_vel, type='l', col='purple', lty=2)


### Plot only comparing CE01 data

# png("Monthly_mean_Surface_Current_Comparison_test2.png", 	width=825 , height= 1840, res= 150)

# # png("CE01_comparison_plot_daily_average.png", width=950 , height= 900, res= 150)

# plot(day_april, ce01_2016$u_vel, type='l', lty=1, col='#DA291CFF', ylim=c(-0.06, 0.06), xlab='time', ylab='u velocity (m/s)', lwd=2)

# lines(day_april, ce01_2017$u_vel, type='l', col='#56A8CBFF', lwd=2)
# lines(day_april, ce01_2018$u_vel, type='l', col='#53A567FF', lwd=2)


# lines(day_april, ROMS_ce01_2016$u_vel, type='l', col='#DA291CFF', lty=2, lwd=2)
# lines(day_april, ROMS_ce01_2017$u_vel, type='l', col='#56A8CBFF', lty=2, lwd=2)
# lines(day_april, ROMS_ce01_2018$u_vel, type='l', col='#53A567FF', lty=2, lwd=2)

# dev.off()



# ### Plot stacked figure for CE01 vs ROMS 

# png("CE01_comparison_stacked_plot_daily_average_u_vel.png", width=1050 , height= 1300, res= 150)

# par(mfrow=c(3,1), par(mai=c(0.4,0.7,0.2,0.4)))
# plot(day_april, ce01_2016$u_vel, type='l', lty=1, col='#DA291CFF', xlab='', ylab='u velocity (m/s)', lwd=2, ylim=c(-0.06,0.05), cex.lab=1.5, cex.axis=1.5)
# lines(day_april, ROMS_ce01_2016$u_vel, type='l', col='#DA291CFF', lty=2, lwd=2)
# mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2)
# legend(x='bottomleft', y=NULL, legend=c('CE01ISSM', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n')

# plot(day_april, ce01_2017$u_vel, type='l', col='#56A8CBFF', lwd=2, ylim=c(-0.06,0.05), xlab='', ylab='u velocity (m/s)', cex.lab=1.5, cex.axis=1.5)
# lines(day_april, ROMS_ce01_2017$u_vel, type='l', col='#56A8CBFF', lty=2, lwd=2)
# mtext('2017', side=3, line=-1.7, adj=0.01, cex=1.2)

# plot(day_april, ce01_2018$u_vel, type='l', col='#53A567FF', lwd=2, ylim=c(-0.06,0.05), xlab='time (days)', ylab='u velocity (m/s)', cex.lab=1.5, cex.axis=1.5)
# lines(day_april, ROMS_ce01_2018$u_vel, type='l', col='#53A567FF', lty=2, lwd=2)
# mtext('2018', side=3, line=-1.7, adj=0.01, cex=1.2)

# dev.off()




# ### Plot stacked figure for CE02 vs ROMS (u_vel)

# png("CE02_comparison_stacked_plot_daily_average_u_vel.png", width=1050 , height= 1300, res= 150)

# par(mfrow=c(3,1), par(mai=c(0.4,0.7,0.2,0.4)))

# plot(day_april, ce02_2016$u_vel, type='l', lty=1, col='#DA291CFF', xlab='', ylab='u velocity (m/s)', lwd=2, ylim=c(-0.23,0.2), cex.lab=1.5, cex.axis=1.5)
# lines(day_april, ROMS_ce02_2016$u_vel, type='l', col='#DA291CFF', lty=2, lwd=2)
# mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2)
# legend(x='bottomleft', y=NULL, legend=c('CE02SHSM', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n')

# plot(day_april, ce02_2017$u_vel, type='l', col='#56A8CBFF', lwd=2, ylim=c(-0.23,0.2), xlab='', ylab='u velocity (m/s)', cex.lab=1.5, cex.axis=1.5)
# lines(day_april, ROMS_ce02_2017$u_vel, type='l', col='#56A8CBFF', lty=2, lwd=2)
# mtext('2017', side=3, line=-1.7, adj=0.01, cex=1.2)




#################
# large panel plot for to plot u and v (but separately)
# dev.copy(jpeg,'Jack_pheno.jpg',height=7,width=5,res=200,units='in')
# dev.off()
# quartz(width=7,height=10)
# par(mfcol=c(3,2))

# png("comparison_stacked_plot_daily_average_ALL_u_vel_test2.png", width=1500 , height= 1300, res= 150)



quartz(width=7,height=6)
par(mfrow=c(3,2), par(mai=c(0.55,0.6,0.2,0.18)))


dev.copy(jpeg,'FINAL_VVEL_PANEL_PLOT.jpeg',width=7,height=6,res=200,units='in')

# 2016
plot(day_2016, ce01_2016$v_vel, type='l', lty=1, col='#DA291CFF', xlab='', ylab='', lwd=2, ylim=c(-0.6,0.4), cex.lab=1.5, cex.axis=1.5, xaxt='n')
lines(day_2016, ROMS_ce01_2016$v_vel, type='l', col='#DA291CFF', lty=2, lwd=2)
mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2, col='#DA291CFF', font=2)
legend(x='bottomleft', y=NULL, legend=c('OOI', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n', lwd=2)
mtext('CE01ISSM - Inshore', side=3, cex=1.2, line=-0.01)
axis(side=1,labels=F)


plot(day_2016, ce02_2016$v_vel, type='l', lty=1, col='#DA291CFF', xlab='', ylab='', lwd=2, ylim=c(-0.6,0.4), cex.lab=1.5, cex.axis=1.5, xaxt='n')
lines(day_2016, ROMS_ce02_2016$v_vel, type='l', col='#DA291CFF', lty=2, lwd=2)
# mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2)
# legend(x='bottomleft', y=NULL, legend=c('CE02SHSM', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n')
mtext('CE02SHSM - Shelf', side=3, cex=1.2)
axis(side=1,labels=F)


# 2017
plot(day_2017, ce01_2017$v_vel, type='l', col='#56A8CBFF', lwd=2, ylim=c(-0.6,0.4), xlab='', ylab='v velocity (m/s)', cex.lab=1.5, cex.axis=1.5, xaxt='n')
lines(day_2017, ROMS_ce01_2017$v_vel, type='l', col='#56A8CBFF', lty=2, lwd=2)
mtext('2017', side=3, line=-1.7, adj=0.01, cex=1.2, col='#56A8CBFF', font=2)
axis(side=1,labels=F)


plot(day_2017, ce02_2017$u_vel, type='l', col='#56A8CBFF', lwd=2, ylim=c(-0.6,0.4), xlab='Time (days)', ylab='', cex.lab=1.5, cex.axis=1.5)
lines(day_2017, ROMS_ce02_2017$u_vel, type='l', col='#56A8CBFF', lty=2, lwd=2)


# 2018
plot(day_2018, ce01_2018$v_vel, type='l', col='#53A567FF', lwd=2, ylim=c(-0.6,0.4), xlab='Time (days)', ylab='', cex.lab=1.5, cex.axis=1.5)
lines(day_2018, ROMS_ce01_2018$v_vel, type='l', col='#53A567FF', lty=2, lwd=2)
mtext('2018', side=3, line=-1.7, adj=0.01, cex=1.2, col='#53A567FF', font=2)

dev.off()



#################
# large panel plot for to plot ocean temp (but separately)


setwd("/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/OOI_data_and_comparison_data/Daily_Averaged_Temperature_Data")


# setwd('/Users/jennifer/Desktop/Groundtruthing_Data/Point_velocity_OOI/Daily_Averages')

### CE01 
ce01_2016<-read.csv("Daily_Averages_CE012016.csv")
dim(ce01_2016)

# ce01_2017<-read.csv("Daily_Averages_CE012017.csv")
# dim(ce01_2017)

ce01_2018<-read.csv("Daily_Averages_CE012018.csv")
dim(ce01_2018)

### ROMS - CE01
ROMS_ce01_2016<-read.csv("Daily_Averages_ROMS_CE012016.csv")
dim(ROMS_ce01_2016)

# ROMS_ce01_2017<-read.csv("Daily_velocity_Averages_CE012017_ROMS.csv")
# dim(ROMS_ceo1_2017)

ROMS_ce01_2018<-read.csv("Daily_Averages_ROMS_CE012018.csv")
dim(ROMS_ce01_2018)
colnames(ROMS_ce01_2018)<-c('X', 'temp')

### CE02 
ce02_2016<-read.csv("Daily_Averages_CE022016.csv")
dim(ce02_2016)

ce02_2017<-read.csv("Daily_Averages_CE022017.csv")
dim(ce02_2017)

ce02_2018<-read.csv("Daily_Averages_CE022018.csv")
dim(ce02_2018)


### ROMS - CE02

ROMS_ce02_2016<-read.csv("Daily_Averages_ROMS_CE022016.csv")
dim(ROMS_ce02_2016)

ROMS_ce02_2017<-read.csv("Daily_Averages_ROMS_CE022017.csv")
dim(ROMS_ce02_2017)
colnames(ROMS_ce02_2017)<-c('X', 'temp')

ROMS_ce02_2018<-read.csv("Daily_Averages_ROMS_CE022018.csv")
dim(ROMS_ce02_2018)
colnames(ROMS_ce02_2018)<-c('X', 'temp')

########

setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/OOI_data_and_comparison_data/')


png("comparison_stacked_plot_daily_average_temp.png", width=1500 , height= 1300, res= 150)

layout(mat = matrix(c(1,0,4,2,3,5),nrow = 3, ncol = 2))     # Widths of the two columns

par(mai=c(0.3,0.55,0.2,0.18))

# 2016
plot(day_april, ce01_2016$temp, type='l', lty=1, col='#DA291CFF', xlab='', ylab='', lwd=2, cex.lab=1.5, cex.axis=1.5, xaxt='n', ylim=c(10.3,14.5))
lines(day_april, ROMS_ce01_2016$temp, type='l', col='#DA291CFF', lty=2, lwd=2)
mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2, col='#DA291CFF', font=2)
legend(x='bottomleft', y=NULL, legend=c('OOI', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n', lwd=2)
mtext('CE01ISSM - Inshore', side=3, cex=1.2, line=-0.01)
axis(side=1,labels=F)

plot(day_april, ce02_2016$temp, type='l', lty=1, col='#DA291CFF', xlab='', ylab='', lwd=2, cex.lab=1.5, cex.axis=1.5, xaxt='n', ylim=c(10.5,14.5))
lines(day_april, ROMS_ce02_2016$temp, type='l', col='#DA291CFF', lty=2, lwd=2)
# mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2)
# legend(x='bottomleft', y=NULL, legend=c('CE02SHSM', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n')
mtext('CE02SHSM - Shelf', side=3, cex=1.2)
axis(side=1,labels=F)


# 2017
plot(day_april, ce02_2017$temp, type='l', col='#56A8CBFF', lwd=2, xlab='', ylab='degC', cex.lab=1.5, cex.axis=1.5, xaxt='n', ylim=c(10,11.6))
lines(day_april, ROMS_ce02_2017$temp, type='l', col='#56A8CBFF', lty=2, lwd=2)
mtext('2017', side=3, line=-1.7, adj=0.01, cex=1.2,col='#56A8CBFF', font=2)
axis(side=1,labels=F)


# 2018
plot(day_april, ce01_2018$temp, type='l', col='#53A567FF', lwd=2, xlab='time (days)', ylab='', cex.lab=1.5, cex.axis=1.5, ylim=c(9,11.5))
lines(day_april, ROMS_ce01_2018$temp, type='l', col='#53A567FF', lty=2, lwd=2)
mtext('2018', side=3, line=-1.7, adj=0.01, cex=1.2,col='#53A567FF',font=2)

plot(day_april, ce02_2018$temp, type='l', col='#53A567FF', lwd=2, xlab='', ylab='', cex.lab=1.5, cex.axis=1.5, ylim=c(9,11))
lines(day_april, ROMS_ce02_2018$temp, type='l', col='#53A567FF', lty=2, lwd=2)

dev.off()







#################
# large panel plot for to plot salinity (but separately)

setwd('/Volumes/TOSHIBA EXT/Desktop_backup/Groundtruthing_Data/Point_velocity_OOI/Daily_Averages')

# setwd("/Users/jennifer/Desktop/Groundtruthing_Data/OOI_data_and_comparison_data/Daily_Averaged_Salinity_Data")

# setwd('/Users/jennifer/Desktop/Groundtruthing_Data/Point_velocity_OOI/Daily_Averages')

### CE01 
ce01_2016<-read.csv("Daily_salinity_Averages_CE012016.csv")
dim(ce01_2016)

# ce01_2017<-read.csv("Daily_Averages_CE012017.csv")
# dim(ce01_2017)

ce01_2018<-read.csv("Daily_salinity_Averages_CE012018.csv")
dim(ce01_2018)

### ROMS - CE01
ROMS_ce01_2016<-read.csv("Daily_salinity_Averages_ROMS_CE012016.csv")
dim(ROMS_ce01_2016)

# ROMS_ce01_2017<-read.csv("Daily_velocity_Averages_CE012017_ROMS.csv")
# dim(ROMS_ceo1_2017)

ROMS_ce01_2018<-read.csv("Daily_salinity_Averages_ROMS_CE012018.csv")
dim(ROMS_ce01_2018)
# colnames(ROMS_ce01_2018)<-c('X', 'temp')

### CE02 
ce02_2016<-read.csv("Daily_salinity_Averages_CE022016.csv")
dim(ce02_2016)

ce02_2017<-read.csv("Daily_salinity_Averages_CE022017.csv")
dim(ce02_2017)

ce02_2018<-read.csv("Daily_salinity_Averages_CE022018.csv")
dim(ce02_2018)


### ROMS - CE02

ROMS_ce02_2016<-read.csv("Daily_salinity_Averages_ROMS_CE022016.csv")
dim(ROMS_ce02_2016)

ROMS_ce02_2017<-read.csv("Daily_salinity_Averages_ROMS_CE022017.csv")
dim(ROMS_ce02_2017)
# colnames(ROMS_ce02_2017)<-c('X', 'temp')

ROMS_ce02_2018<-read.csv("Daily_salinity_Averages_ROMS_CE022018.csv")
dim(ROMS_ce02_2018)
# colnames(ROMS_ce02_2018)<-c('X', 'temp')



########

setwd('/Users/jennifer/Desktop/Groundtruthing_Data/OOI_data_and_comparison_data/')


png("comparison_stacked_plot_daily_average_salinity.png", width=1500 , height= 1300, res= 150)

layout(mat = matrix(c(1,0,4,2,3,5),nrow = 3, ncol = 2))     # Widths of the two columns

par(mai=c(0.3,0.55,0.2,0.18))

# 2016
plot(day_2016, ce01_2016$salt, type='l', lty=1, col='#DA291CFF', xlab='', ylab='', lwd=2, cex.lab=1.5, cex.axis=1.5, xaxt='n', ylim=c(30,34))
lines(day_april, ROMS_ce01_2016$salt, type='l', col='#DA291CFF', lty=2, lwd=2)
mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2, col='#DA291CFF', font=2)
legend(x='bottomleft', y=NULL, legend=c('OOI', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n', lwd=2)
mtext('CE01ISSM - Inshore', side=3, cex=1.2, line=-0.01)
axis(side=1,labels=F)

plot(day_2016, ce02_2016$salt, type='l', lty=1, col='#DA291CFF', xlab='', ylab='', lwd=2, cex.lab=1.5, cex.axis=1.5, xaxt='n', ylim=c(30,34))
lines(day_april, ROMS_ce02_2016$salt, type='l', col='#DA291CFF', lty=2, lwd=2)
# mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2)
# legend(x='bottomleft', y=NULL, legend=c('CE02SHSM', 'ROMS'), lty=c(1,2), cex= 1.4, bty='n')
mtext('CE02SHSM - Shelf', side=3, cex=1.2)
axis(side=1,labels=F)


# 2017
plot(day_2017, ce02_2017$salt, type='l', col='#56A8CBFF', lwd=2, xlab='', ylab='PSU', cex.lab=1.5, cex.axis=1.5, xaxt='n', ylim=c(30,34))
lines(day_2017, ROMS_ce02_2017$salt, type='l', col='#56A8CBFF', lty=2, lwd=2)
mtext('2017', side=3, line=-1.7, adj=0.01, cex=1.2,col='#56A8CBFF', font=2)
axis(side=1,labels=F)


# 2018
plot(day_2018, ce01_2018$salt, type='l', col='#53A567FF', lwd=2, xlab='time (days)', ylab='', cex.lab=1.5, cex.axis=1.5, ylim=c(30,34))
lines(day_2018, ROMS_ce01_2018$salt, type='l', col='#53A567FF', lty=2, lwd=2)
mtext('2018', side=3, line=-1.7, adj=0.01, cex=1.2,col='#53A567FF',font=2)

plot(day_2018, ce02_2018$salt, type='l', col='#53A567FF', lwd=2, xlab='', ylab='', cex.lab=1.5, cex.axis=1.5, ylim=c(30,34))
lines(day_2018, ROMS_ce02_2018$salt, type='l', col='#53A567FF', lty=2, lwd=2)

dev.off()



