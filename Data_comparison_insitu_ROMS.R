### Point comparison between ROMS and in situ 
## Date: 09/30/2020


##########################################################
### clear working directory

rm(list=ls())


##########################################################
### load libraries

library(ncdf4)
# library(OceanView)
library(R.matlab)
library(marmap)
library(raster)
library(rgdal)

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



##########################################################
### load data
	
setwd('/Volumes/TOSHIBA EXT/ROMS_HiRes/2018') # high res ROMS directory
	
files<-list.files('/Volumes/TOSHIBA EXT/ROMS_HiRes/2018', pattern= '*.nc', full.names=T) # list of the ROMS files names to open adn loop through

# setwd('/Volumes/TOSHIBA EXT/April2016_LowRes') # low res ROMS location

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


##########################################################
### plot map to select regions 

# blues <- c("lightsteelblue4", "lightsteelblue3","lightsteelblue2", "lightsteelblue1")
# 
# plot.bathy(OR_bathy2, deep=-500, shallow=-10, step=20,image=T, land=T, lwd=0.1, bpal=list(c(0, max(OR_bathy2), 'grey'), c(min(OR_bathy2), 0, blues)))
# 
# contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(50, 100), labcex=0.7,add=T,col='red', lwd= 0.8)
# 
# points(-124.86, 44.1, pch=19, cex=0.7) # shelf location (HB)
# 
# points(-124.25, 44.1, pch=19, cex=0.7) # inshore location (HB)
# 
# points(-124.86, 42.82, pch=19, cex=0.7) # shelf location (CB)
# 
# points(-124.64, 42.82, pch=19, cex=0.7) # inshore location (CB)


###########################################################
### grid locations to extract vertical velocity data

## low res locations:
# shelf HB = [49,193,40,]
# inshore HB = [73,193,40,]

# shelf CB = [49,122,40,]
# inshore CB = [58,122,40,]


## high res locations
# shelf HB = [78,936,40,]
# inshore HB = [275,936,40,]

# shelf CB = [78,366,40,]
# inshore CB = [149,366,40,]

## dimensions to put into for loop
### change to high res dimensions 

lon_dim<-149

lat_dim<-366

dep_dim<-40


##########################################################
### load data

tmp<-NULL # empty vector to fill data with

for(i in 1:length(files)){
	
	print(i)
	
	nc.files<-nc_open(files[i])
	w_dat<-ncvar_get(nc.files, varid='w') # units= m/s
	fillvalue_w<-ncatt_get(nc.files, 'w', "_FillValue") # replace the NAs with 'NA'
	w_dat[w_dat==fillvalue_w]<-NA 
	w_avg2<-0.5*(w_dat[,,2:depth_exd,] + w_dat[,,1:depthd,])
	w_avg1<-w_avg2[1:latd,1:lond,,]
	w_avg1[is.na(w_avg1)]<-min(w_avg1,na.rm=T)
	
	ROMS_w<-w_avg1[lon_dim,lat_dim,dep_dim,] # extract from a single location
	
	ROMS_time<-ncvar_get(nc.files, varid='ocean_time') # extract time (no need extract because its a single vector)

	b<-cbind(ROMS_w, ROMS_time)
	
	# ROMS_2018time<-rbind(ROMS_2018time, b) 
	tmp<-rbind(tmp, b)
	
} 

setwd('/Users/wongalaj/Desktop')
inshoreCB_hires18<-tmp
save(inshoreCB_hires18, file='vert_vel_hires18_inshoreCB.RData')


##########################################################
### load in post processed data and plot figures 

setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/vertical_velocity_comparisons/')

## low res data
# 2016 = 420
# 2017 = 420
# 2018 = 444

low16<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(inshoreHB_lowres16))
length(low16) # 420
hi16<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(inshoreHB_hires16))
length(hi16)

low17<-seq(as.POSIXct("2017-03-27 00:00:00"), as.POSIXct("2017-05-01 00:00:00"), length.out=nrow(inshoreHB_lowres17))
length(low17)
hi17<-seq(as.POSIXct("2017-03-27 00:00:00"), as.POSIXct("2017-05-01 00:00:00"), length.out=nrow(inshoreHB_hires17))
length(hi17)

low18<-seq(as.POSIXct("2018-03-27 00:00:00"), as.POSIXct("2018-05-01 00:00:00"), length.out=nrow(inshoreHB_lowres18))
length(low18)
hi18<-seq(as.POSIXct("2018-03-27 00:00:00"), as.POSIXct("2018-05-01 00:00:00"), length.out=nrow(inshoreHB_hires18))
length(hi18)

##########################################################
### plot data 

title_use<-'all_plots_2018'

year<-2018
ROMS_title<-'2 km'

ROMS_title2<-'250 m'

text_title<-paste(year, ROMS_title, sep=' - ')
text_title2<-paste(year, ROMS_title2, sep=' - ')

yr_plt<-low18

yr_plt2<-hi18

# low res
HB_shelf<-shelfHB_lowres18
HB_inshore<-inshoreHB_lowres18

CB_shelf<-shelfCB_lowres18
CB_inshore<-inshoreCB_lowres18

# high res
HB_shelf2<-shelfHB_hires18
HB_inshore2<-inshoreHB_hires18

CB_shelf2<-shelfCB_hires18
CB_inshore2<-inshoreCB_hires18

# 2016 = (-6e-04, 9e-04)
# 2017 = (-8e-04, 4e-04)
# 2018 = (-2e-03, 6.6e-04)

quartz(width=9, height=7)
dev.copy(jpeg, paste(title_use, '.jpg', sep=''), height=7, width=9, res=200, units='in')
par(mfrow=c(2,2), mai=c(0.8, 0.8, 0.6, 0.6))

plot(yr_plt, HB_inshore[,1], type='l', xlab='', ylab='', main='Heceta Bank', ylim=c(-2e-03, 6.6e-04)) # 
lines(yr_plt, HB_shelf[,1], type='l', xlab='', ylab='', lty=2, col='red')
legend(x="topright", y=NULL, legend=c("inshore","shelf"), lty=c(1,2), col=c('black', 'red'), cex=0.7, bty='n')
mtext(text_title, side=3, line=0, adj=0)

plot(yr_plt, CB_inshore[,1], type='l', xlab='', ylab='', main='Cape Blanco', ylim=c(-2e-03, 6.6e-04))
lines(yr_plt, CB_shelf[,1], type='l', xlab='', ylab='',  lty=2, col='red')

plot(yr_plt2, HB_inshore2[,1], type='l', xlab='Time (hourly)', ylab='', main='', ylim=c(-2e-03, 6.6e-04)) # 
lines(yr_plt2, HB_shelf2[,1], type='l', xlab='', ylab='', lty=2, col='red')
mtext(text_title2, side=3, line=0, adj=0)


plot(yr_plt2, CB_inshore2[,1], type='l', xlab='Time (hourly)', ylab='', main='', ylim=c(-2e-03, 6.6e-04))
lines(yr_plt2, CB_shelf2[,1], type='l', xlab='', ylab='', lty=2, col='red')
mtext('Vertical velocity (m/s)', side=2, line=30, adj=4)

dev.off()






