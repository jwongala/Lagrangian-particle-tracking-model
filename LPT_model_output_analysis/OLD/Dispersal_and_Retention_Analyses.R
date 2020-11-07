### Jennnifer A.T.K. Wong-Ala
## MS Thesis
## Particle Tracking Model Code for Oregon Coast (including Marine Reserves, Yaquina Bay and Estuary)
## Last worked on: 03/26/2019 (month/day/year format)

####################################################################
##### NOTES: Read before continuing to use script!

# 1. Sections noted with (***) contain placeholders or are notes that Jenn needs to attend to
# 2. LORENZO: please change your directory so it is correct. Thanks!

####################################################################
### Prep the workspace
rm(list=ls()) # clear the workspace
options(warn=0) # warn= 0 will let me know there is an error, but won't stop the code from running; use warn=2 to stop code from running when error occurs


####################################################################
distance.function<-function(start.lat, start.lon, end.lat, end.lon)
{
med.lat <- (start.lat + end.lat)/2
rad.lat <- (pi * med.lat)/180
shrink <- cos(rad.lat)
delta.lat <- end.lat - start.lat
delta.lon <- start.lon - end.lon
mpermile <- 111195
distance <- mpermile * sqrt((delta.lon * shrink)^2 + (delta.lat)^2)
distance
}

### Load libraires
library(marmap) # to plot bathymetry
library(ncdf4) # to load netcdf advection files
library(data.table) # to use fread
library(mgcv) # to do GAMs
library(sp)
library(raster)
library(ncdf.tools)
library(rgeos)
library(maptools)
library(grDevices)
library(fields)
library(pracma)
library(OceanView)
library(animation)
library(R.matlab)
library(plot3D)
library(lattice)
library(plotly)
library(maps)
library(mapdata)
library(fields)
library(marmap)
library(sgeostat)
library(raster)
library(rgdal)
library(maptools)
library(colorspace)
library(RColorBrewer)
library(ggplot2)
library(tiff)

####################################################################
### Constants

# Radius of the Earth  
R<-6378137 # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method


####################################################################
### load RData file to plot
# I think I can only read one of these at a time because it will over write the previous ones that's why. Not totally sure. 


# set back to original directory
setwd("/Users/jennifer/Desktop/model_output_chapter1&2/used_high_res_ROMS/2018") 
	### /Users/jennifer/Desktop
	### /Users/wongalaj/Desktop/April2016_LowRes
	### /Volumes/TOSHIBA EXT/OR_coast_AK_ROMS
	### /Users/Jenn/Desktop/Model
	### /Volumes/TOSHIBA EXT/ROMS_HiRes/2016
setwd('/Users/jennifer/Desktop/RET_FO_data/FO')

	
# load("FO_high_mean.RData")

# highfo<-avg_test
# max(highfo, na.rm=T) 

# fo max= 4.130482
# ret max= 1.511111

####################################################################
### load initial locations for particles

setwd("/Users/jennifer/Desktop/model_output_chapter1&2")
load("initial_locations_2km_ROMS.RData")

# init_par - init locations for all particles

# # init_surf<-subset(init_par, ID %in% seq(1,250,1))

# init_mid<-subset(init_par, ID %in% seq(251,500,1))

# init_deep<-subset(init_par, ID %in% seq(501,750,1))

### save initial locations for particles 

# setwd("/Users/jennifer/Desktop/model_output_chapter1&2/used_low_res_ROMS/2016/OLD") 

# load('LPT_06.18.2020_JWA_6hour_2016_low.RData')

# init_par<-pop.roms[[1]][[1]]

# save(init_par, file= "initial_locations_2km_6hour_ROMS.RData")


####################################################################
### load files 

### Bathymetry

setwd('/Users/jennifer/Desktop/')
str_name<-"large_domain.tiff"
OR_bathy<-raster(str_name)

library(graphics)
OR_bathy2<-as.bathy(OR_bathy)


### Marine Reserve Outlines

setwd("/Users/jennifer/Desktop/")

MR<-readOGR(dsn=path.expand("/Users/jennifer/Desktop/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"), layer= 'MPA_MR_COMP_Boundaries_UTM10N')

class(MR) # SpatialPolygonsDataFrame
crs(MR)

extent(MR)

crsmerc=CRS("+proj=longlat +lat_1=43 +lat_2=48 +lat_0=41 +lon_0=-117 +x_0=700000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ") # this transforms the model 

MR_transformed<-spTransform(MR, CRS=crsmerc)

### Bathymetry

setwd("/Users/jennifer/Desktop/")

bathy.dat<-read.table('etopo1_bedrock.xyz', sep='')
names(bathy.dat)<-c('lon', "lat", 'depth')

bathy.dat$depth[bathy.dat$depth>0]<-NA # Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]


####################################################################
### Convert km size of MR to degrees so that I can divide up coastal ocean into grid cells that are no bigger than the largest size of the Otter Rock MR

lat1<-42 
lat2<-45.24775

lon1<-(-125.10) # -126.0909
lon2<-(-123.9263)  

Rad<-6378137 # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method

# OR.size<-3 # km2; size of Otter Rock MR

deg.lat<-lat2-lat1 # length of coastline measured in degrees 

deg.lon<-lon2-lon1

## conversion of coastline length in degrees to km

lat.len<-Rad*pi*(deg.lat/180)*0.001 # length of coastline (latitude in km)

lon.len<-(Rad*pi*(cos(pi*(lat.len/180)))*(deg.lon/180))*0.001 # width of coastline (longitude in km)


####################################################################
### Remove part of trajectory that exceeds most least common del.t for all particles (23 days=2208 del.t)

fin_dat<-subset(data_use, del.t<=1344) # only contains particle trajectories for 14 days or less (1344 del.t)

# dim(fin_dat)

fin_end<-subset(data_use, del.t==1344) # only contains end location of particles 

####################################################################
### subset particles based on depth layer release

# # # surface, middle, deep
# surf_par<-subset(fin_dat, ID %in% seq(1,250,1))
# dim(surf_par) # 179386     16
# surf_end<-subset(fin_end, ID %in% seq(1,250,1))

# mid_par<-subset(fin_dat, ID %in% seq(251,500,1))
# dim(mid_par) # 179386     16
# mid_end<-subset(fin_end, ID %in% seq(251,500,1))

# deep_par<-subset(fin_dat, ID %in% seq(501,750,1))
# dim(deep_par) # 179386     16
# deep_end<-subset(fin_end, ID %in% seq(501,750,1))


####################################################################
### Making Frequency of Occurrence Heat maps

# change what data set using to create the plots

freq_dat<-fin_dat
	# fin_dat - all par trajectories
	# surf_par
	# mid_par
	# deep_par
	
title_foo<-"2018 mean FO" # title to be used for figures 
	# "2km ROMS, all particles, initial, 23 day - FoO"
	# "2km ROMS, surface, initial, 23 day - FoO"
	# "2km ROMS, middle, initial, 23 day - FoO"
	
	
#Part 1: make a regular grid and count stations within each grid cell. Plot results
   
    nlat= 360/10 # lat.len/13 # 698/10 # determine resolution of grid
    nlon= 130/10 # lon.len/13 # 238/10
    latd=seq(lat1, lat2,length.out=nlat)
    lond=seq(lon1,lon2,length.out=nlon)
    
    
  grid.lon=data.frame(
    lon1=rep(lond[-length(lond)],(nlat-1)),
    lon2=rep(lond[-1],(nlat-1)),
    lon3=rep(lond[-1],(nlat-1)),
    lon4=rep(lond[-length(lond)],(nlat-1)))#make dataframe of just longitude
  
  grid.lat=data.frame(
    lat1=sort(rep(latd[-length(latd)],(nlon-1))),
    lat2=sort(rep(latd[-length(latd)],(nlon-1))),
    lat3=sort(rep(latd[-1],(nlon-1))),
    lat4=sort(rep(latd[-1],(nlon-1)))) # lat dataframe
   
  ### sort init_pars data
  dev.new(width=4,height=10)
  plot(freq_dat$lon, freq_dat$lat,pch='.',ylim=c(lat1, lat2),xlim=c(lon1,lon2))
  # points(maynew2$lon, maynew2$lat, col='red', pch= '.')
  n.stations=NA*(1:nrow(grid.lon))
   
  for(i in 1:length(n.stations)){
    print(i)
    tmp=in.chull(freq_dat$lon, freq_dat$lat, grid.lon[i,], grid.lat[i,])
    n.stations[i]=sum(tmp)#This decides what goes into each grid pixel (may16_3$ID*tmp)
    points(freq_dat$lon[tmp], freq_dat$lat[tmp],col=i,pch=16)
    polygon(grid.lon[i,],grid.lat[i,])
  }
  map("worldHires",fill=T,col="grey",add=T)
  
num_ts<-nrow(freq_dat)


## creating grid to plot for heatmap 
z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2 # 
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2 # 
z.mat<-((matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F))/num_ts)*100 # /uniID_nov16 # units are percentage (%)
range(z.mat, na.rm=T) # 21 117

z.matFO_24<-z.mat
save(z.matFO_24, file="FO_24.RData")
# write.csv(z.matrix_36days_FoO, "FreqofOccur_36dayrelease_lowres_20depth_HiResCompare.csv")

x<-list(z.matFO_21, z.matFO_22, z.matFO_23, z.matFO_24)
x1<-do.call(cbind,x)
x2<-array(x1, dim=c(dim(x[[1]]), length(x)))

avg_test<-apply(x2, c(1,2), mean, na.rm=T)
fo18_hi<-avg_test

save(fo18_hi, file='fo18_hi_mean.RData')

mean16<-mean(fo16_low)
mean17<-mean(fo17_low)
mean18<-mean(fo18_low)

himean16<-mean(fo16_hi)
himean17<-mean(fo17_hi)
himean18<-mean(fo18_hi)

# max(FO_2016_mean, na.rm=T) 
# 2016= 6.376574
# 2017= 3.22395
# 2018= 3.932909

# # # ### plot data 
# fo max= 4.130482
# ret max= 1.511111

pal<-colorRampPalette(c("white","cyan", "cyan4", "darkblue"))

# setwd("/Users/wongalaj/Desktop/Chapter2_Figures")

dev.new(width= 6, height=10, res=200)
par(mai=c(1.1,1.1,0.4,0.3))

image.plot(z.lon,z.lat, fo18_hi,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main= "fo18_hi", col=pal(500), cex.lab=1.5,cex.axis=1.5, cex.main=1.2, zlim=c(0,7)) 

map("worldHires",fill=T,col="grey", add=T)

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,3000,by=300),labcex=0.7,add=T,col='black', lwd= 0.1) 

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
  
text(-124.46,42.85,"CB", cex=1.5)
text(-123.97, 44.6368, 'N', cex=1.5)

####################################################################
### Import plots

# # # change what data set using to create the plots
# imp_dat<-fin_end
	# # fin_end 
	# # surf_end 
	# # mid_end 
	# # deep_end
	
# title_imp<-"2 km mean all" # title to be used for figures

	# # "2km ROMS, all pars, initial, 23 day - Imp"
	# # "2km ROMS, surface, initial, 23 day - Imp"
	
# ## create grid of only initial locations for lat and lon
    # nlat.end= 360/10 # 698/10 # determine resolution of grid
    # nlon.end= 130/10 # 238/10
    # latd.end=seq(lat1, lat2,length.out=nlat.end)
    # lond.end=seq(lon1, lon2,length.out=nlon.end) # -124.2
    
  # gridlon.end=data.frame(
    # lon1=rep(lond.end[-length(lond.end)],(nlat.end-1)),
    # lon2=rep(lond.end[-1],(nlat.end-1)),
    # lon3=rep(lond.end[-1],(nlat.end-1)),
    # lon4=rep(lond.end[-length(lond.end)],(nlat.end-1))) # make dataframe of just longitude
  
  # gridlat.end=data.frame(
    # lat1=sort(rep(latd.end[-length(latd.end)],(nlon.end-1))),
    # lat2=sort(rep(latd.end[-length(latd.end)],(nlon.end-1))),
    # lat3=sort(rep(latd.end[-1],(nlon.end-1))),
    # lat4=sort(rep(latd.end[-1],(nlon.end-1)))) # make dataframe of just latitude

# ## subset data so that it only includes last timestep for each particle because I want their end locations

# # plot what data looks like for each end location subset
# plot(imp_dat$lon, imp_dat$lat, col='red', pch= '.')
# n.stations.end=NA*(1:nrow(gridlon.end))

# # figure out what grid cell is in each 
# for(i in 1:length(n.stations.end)){
	# print(i)
	# tmp.end=in.chull(imp_dat$lon, imp_dat$lat, gridlon.end[i,], gridlat.end[i,])
	# n.stations.end[i]=sum(tmp.end)#This decides what goes into each grid pixel (may16_3$ID*tmp)
	# points(imp_dat$lon[tmp.end], imp_dat$lat[tmp.end],col=i,pch=16)
	# polygon(gridlon.end[i,], gridlat.end[i,])
	# }


# ## creating grid to plot for heatmap 
# z.lat.end<-(latd.end[1:(length(latd.end)-1)]+latd.end[2:length(latd.end)])/2
# z.lon.end<-(lond.end[1:(length(lond.end)-1)]+lond.end[2:length(lond.end)])/2
# z.mat.end<-(matrix(n.stations.end,ncol=length(z.lat.end),nrow=length(z.lon.end),byrow=F)) 
# dim(z.mat.end) # 21 117
# range(z.mat.end)


# # plot heatmap of import data

# setwd("/Users/jennifer/Desktop/Chapter2_Figures")

# pal2<-colorRampPalette(c("white","mediumorchid1", "mediumorchid3", "mediumorchid4"))

# dev.new(width= 6, height=10)
# par(mai=c(1.1,1.1,0.4,0.3))

# image.plot(z.lon.end,z.lat.end, z.mat.end ,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main= title_imp, col=pal2(500), cex.lab=1.5,cex.axis=1.5, cex.main=1.2, zlim=c(0,50)) # zlim= c(0,10)
  
# map("worldHires",fill=T,col="grey", add=T)

# contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,3000,by=300),labcex=0.7,add=T,col='black', lwd= 0.2) 

# lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
  
# text(-124.46,42.85,"CB", cex=1.5)
# text(-123.97, 44.6368, 'N', cex=1.5)
# # # range(z.matrix.IMP14, na.rm=T)

 
####################################################################
### Retention of particles

## - Def: the percent of particles that are returned to coast line retained 
# use the below dataframes as end locations for particles

lr_start<-init_par
	# init_par - init location of all pars  	
	# init_surf 
	# init_mid 
	# init_deep

lr_end<-fin_end
	# fin_end - all par 
	# surf_end 
	# mid_end 
	# deep_end

title_lr<-"2km mean all RET"
	# "2km ROMS, all pars, initial, 23 day - LR"
	# "2km ROMS, surface, initial, 23 day - LR"

### lat and lon grid needed to figure out if particles are in a certain grid cell 
 	
   	nlat.local=360/10 # grid.size.lat/2
    nlon.local=130/10 # grid.size.lat/2
    latd.local= seq(lat1, lat2,length.out= nlat.local) # 698/10 # determine resolution of grid
   	lond.local= seq(lon1, lon2,length.out= nlon.local) # 238/10
    
  gridlon.local=data.frame(
    lon1=rep(lond.local[-length(lond.local)],(nlat.local-1)),
    lon2=rep(lond.local[-1],(nlat.local-1)),
    lon3=rep(lond.local[-1],(nlat.local-1)),
    lon4=rep(lond.local[-length(lond.local)],(nlat.local-1))) # make dataframe of just longitude
  
  gridlat.local=data.frame(
    lat1=sort(rep(latd.local[-length(latd.local)],(nlon.local-1))),
    lat2=sort(rep(latd.local[-length(latd.local)],(nlon.local-1))),
    lat3=sort(rep(latd.local[-1],(nlon.local-1))),
    lat4=sort(rep(latd.local[-1],(nlon.local-1)))) # make dataframe of just latitude


# plot start and end locations on map for each day release/depth subset
plot(lr_end$lon, lr_end$lat, col='red', pch= '.', ylim=c(lat1, lat2),xlim=c(lon1, lon2))
points(lr_start$lon, lr_start$lat, col='blue', pch= '.')

n.stations.start=NA*(1:nrow(gridlon.local)) # empty vector to be filled 
n.stations.end=NA*(1:nrow(gridlon.local)) # empty vector to be filled 

### count particles for start locations (this always stays as pop_init data.frame)
for(i in 1:length(n.stations.start)){
	print(i)
	tmp.local=in.chull(lr_start$lon, lr_start$lat, gridlon.local[i,], gridlat.local[i,])
	n.stations.start[i]=sum(tmp.local) # This decides what goes into each grid pixel
	points(lr_start$lon[tmp.local], lr_start$lat[tmp.local],col=i,pch=16)
	polygon(gridlon.local[i,], gridlat.local[i,])
	} 

### count particles for end locations
for(i in 1:length(n.stations.end)){
	print(i)
	tmp.local=in.chull(lr_end$lon, lr_end$lat, gridlon.local[i,], gridlat.local[i,])
	n.stations.end[i]=sum(tmp.local) # This decides what goes into each grid pixel
	points(lr_end$lon[tmp.local], lr_end$lat[tmp.local],col=i,pch=16)
	polygon(gridlon.local[i,], gridlat.local[i,])
	} 

## creating grid to plot for start heatmap
z.lat.local<-(latd.local[1:(length(latd.local)-1)]+latd.local[2:length(latd.local)])/2 # lat grid to plot heatmap
z.lon.local<-(lond.local[1:(length(lond.local)-1)]+lond.local[2:length(lond.local)])/2 # lon grid to plot heatmap
z.matrix.start<-(matrix(n.stations.start,ncol=length(z.lat.local),nrow=length(z.lon.local),byrow=F))
dim(z.matrix.start) # 12 35
range(z.matrix.start) 


z.matrix.end2<-(matrix(n.stations.end,ncol=length(z.lat.local),nrow=length(z.lon.local),byrow=F))
dim(z.matrix.end2) # 10 58
range(z.matrix.end2)

# create retention metric to plot (difference between # of end particles/# of start particles)
z.matrix.LR<- z.matrix.end2/z.matrix.start # local recruitment heatmap
dim(z.matrix.LR) # 12 35
range(z.matrix.LR, na.rm=T) # 0 4 
z.matrix.LR[which(z.matrix.LR=='Inf')]<-NaN # make Inf to NaN
range(z.matrix.LR, na.rm=T)

z.matLR_24<-z.matrix.LR
save(z.matLR_24, file="RET_24.RData")


y<-list(z.matLR_21, z.matLR_22, z.matLR_23, z.matLR_24)
y1<-do.call(cbind,y)
y2<-array(y1, dim=c(dim(y[[1]]), length(y)))

avg_test2<-apply(y2, c(1,2), mean, na.rm=T)
ret18_hi<-avg_test2

save(ret18_hi, file='ret18_hi.RData')
mean(ret18_hi, na.rm=T)
# max(RET_2018_mean, na.rm=T) 
# 2016= 1.733333
# 2017= 0.7333333
# 2018= 0.8666667

# save(avg_test2, file='RET_high_mean.RData')
# load('RET_low_mean.RData')

# plot heat map of local recruitment metrics (plotted for different depths and release lengths)

pal3<-colorRampPalette(c('white', 'darkorange', 'darkorange1', 'darkorange2', 'darkorange3', 'darkorange4'))

setwd("/Users/jennifer/Desktop/Chapter2_Figures")

dev.new(width= 6, height=10)
par(mai=c(1.1,1.1,0.4,0.3))

image.plot(z.lon.local, z.lat.local, ret18_hi, xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main= "ret18_hi", col= pal3(500), cex.lab=1.5,cex.axis=1.5, cex.main=1.2, zlim=c(0, 3)) 

map("worldHires",fill=T,col="grey", add=T)

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,3000,by=300),labcex=0.7,add=T,col='black', lwd= 0.2) 

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons

text(-124.46,42.85,"CB", cex=1.5)
text(-123.97, 44.6368, 'N', cex=1.5)
# # range(z.matrix.IMP14, na.rm=T)


####################################################################
### Averaged distance travelled

fin_trav<-fin_dat # entire data set to be looking at 
init_trav<-init_par # only initial locations of data set

end_trav<-fin_end # only end locations of data set

empty<-NA*1:length(unique(end_trav$ID))
empty2<-NA*1:length(unique(end_trav$ID))
par_ids<-unique(end_trav$ID)

for(i in 1:length(empty)){
	begin<-subset(init_trav, ID %in% par_ids[i])
	begin_lat<-begin[1,1]
	begin_lon<-begin[1,2]
	
	end<-subset(end_trav, ID %in% par_ids[i]) 
	end_lat<-end[1,1] 
	end_lon<-end[1,2] 
	
	dlat<-begin_lat-end_lat
	dlon<-begin_lon-end_lon
	
	empty[i]<-dlat # distance in degrees
	empty2[i]<-dlon 
}

empty<-3.24774
empty2<-1.6

# convert degrees to km (negative= travel north, positive=travel south) 

test_lat<-(pi*empty/180)

shrink<-cos(test_lat)
mpmile<-111195

test_lon<-mpmile*sqrt((empty2*shrink)^2 + (test_lat)^2) # meters

dist_lon<-test_lon/1000 # lon dist travelled in km 
dist_lat<-Rad*pi*(empty/180)*0.001

tot_dist<-as.data.frame(cbind(par_ids, dist_lat, dist_lon))

##
avg_lat<-mean(tot_dist$dist_lat)
sd_lat<-sd(tot_dist$dist_lat)

err_lat<-qnorm(0.95)*sd_lat/sqrt(length(par_ids))
err_lat

l_lat<-avg_lat-err_lat
r_lat<-avg_lat+err_lat

##
avg_lon<-mean(tot_dist$dist_lon)
sd_lon<-sd(tot_dist$dist_lon)

err_lon<-qnorm(0.95)*sd_lon/sqrt(length(par_ids))
err_lon

l_lon<-avg_lon-err_lon
r_lon<-avg_lon + err_lon


### subset surf, mid, deep
surf_dist<-subset(tot_dist, par_ids %in% seq(1,250,1))
surf_ind<-rep(1, length.out=nrow(surf_dist))
surf_combo<-cbind(surf_dist, surf_ind)
names(surf_combo)<-c("par_ids", "dist_lat", "dist_lon", "ind")

mid_dist<-subset(tot_dist, par_ids %in% seq(251,500,1))
mid_ind<-rep(2, length.out=nrow(mid_dist))
mid_combo<-cbind(mid_dist,mid_ind )
names(mid_combo)<-c("par_ids", "dist_lat", "dist_lon", "ind")

deep_dist<-subset(tot_dist, par_ids %in% seq(501,750,1))
deep_ind<-rep(3, length.out=nrow(deep_dist))
deep_combo<-cbind(deep_dist, deep_ind)
names(deep_combo)<-c("par_ids", "dist_lat", "dist_lon", "ind")

dist_plot<-rbind(surf_combo, mid_combo, deep_combo)

# boxplot(dist_lat ~ ind, data=dist_plot)

### averages for each subset
avg_surflat<-mean(surf_dist$dist_lat)
avg_surflon<-mean(surf_dist$dist_lon)

avg_midlat<-mean(mid_dist$dist_lat)
avg_midlon<-mean(mid_dist$dist_lon)

avg_deeplat<-mean(deep_dist$dist_lat)
avg_deeplon<-mean(deep_dist$dist_lon)

avg_dist_lat<-rbind(avg_lat, avg_surflat, avg_midlat, avg_deeplat)
avg_dist_lon<-rbind(avg_lon, avg_surflon, avg_midlon, avg_deeplon)

ind<-1:4

dim_name<-c('surface', 'middle', 'deep')
all_dist_dat<-as.data.frame(cbind(avg_dist_lat, avg_dist_lon, ind))
names(all_dist_dat)<-c("avg_lat", "avg_lon", "particles")

####################################################################
### Site source potential metric 

init_use<-init_par[,12:14]


# calculate end lat and lon indexes for end_par dataframe 

for(k in 1:nrow(fin_end)){
	dist<-distance.function(fin_end $lat[k], fin_end $lon[k],lat_rho[,,1],lon_rho[,,1])
	tmp<-1*(dist==min(dist))
	lat.index<-(1:dim(lat_rho)[1])[apply(tmp,1,sum)==1]
	lon.index<-(1:dim(lat_rho)[2])[apply(tmp,2,sum)==1]
	fin_end$lat.index[k]<-lat.index
	fin_end$lon.index[k]<-lon.index
}

fin_end_use<-fin_end[,12:14]

fin_id<-fin_end_use$ID

# subset out same par IDs as in fin_end_use
init_use<-init_par[,12:14]

init_use2<-subset(init_use, ID %in% fin_id)


# lat_empty<-matrix(data=NA, nrow=nrow(fin_end), ncol=nrow(init_use))
# lon_empty<-matrix(data=NA, nrow=nrow(fin_end), ncol=nrow(init_use))

# for(i in 1:nrow(init_use)){
	
	# lat_empty[,i]<-init_use$lat.index[i]-fin_end_use$lat
	# lon_empty[,i]<-init_use$lon.index[i]-fin_end_use$lon
	
# }

# lat_empty<-as.data.frame(lat_empty)

# colnames(lat_empty)<-1:750

####################################################################
### find out what is the max del.t to determine the time at which the particle left the model at
data_use2<-subset(data_use, del.t<=2208)
dim(data_use2)

# all particles in model simulation
par_ID<-unique(data_use2$ID)
tmp<-NA*1:length(par_ID)


# only surface released particles
surf_par<-subset(data_use2, ID %in% seq(1,250,1))
dim(surf_par) # 179386     16
par_ID2<-unique(surf_par$ID)
tmp2<-NA*1:length(par_ID2)

# only middle released particles
mid_par<-subset(data_use2, ID %in% seq(251,500,1))
dim(mid_par) # 177727     16
par_ID3<-unique(mid_par$ID)
tmp3<-NA*1:length(par_ID3)

# only deep released particles
deep_par<-subset(data_use2, ID %in% seq(501,750,1))
dim(deep_par) # 180791     16
par_ID4<-unique(deep_par$ID)
tmp4<-NA*1:length(par_ID4)

## for loop to find which is max del.t for each particle and then put into empty vector for each group

# all particles
for(i in 1:length(par_ID)){
	print(i)
	tmp[i]<-max(data_use2$del.t[which(data_use2$ID %in% par_ID[i])])
}

# surface particles
for(i in 1:length(par_ID2)){
	print(i)
	tmp2[i]<-max(surf_par$del.t[which(surf_par$ID %in% par_ID2[i])])
}

# middle particles
for(i in 1:length(par_ID3)){
	print(i)
	tmp3[i]<-max(mid_par$del.t[which(mid_par$ID %in% par_ID3[i])])
}

# deep particles
for(i in 1:length(par_ID4)){
	print(i)
	tmp4[i]<-max(deep_par$del.t[which(deep_par$ID %in% par_ID4[i])])
}

# convert del.t into days
tmp_tot<-tmp/96
# max(tmp_tot)
tmp2<-tmp2/96
tmp3<-tmp3/96
tmp4<-tmp4/96

# add NA at the end if lengths of tmp's are different
tmp2<-c(tmp2,NA)
# tmp3<-c(tmp3,NA)
tmp4<-c(tmp4,NA)

# combine all into one matrix and then paste into excel to post process 
df<-cbind(tmp2, tmp3, tmp4)

####################################################################
### load variables

# colors for bathymetry
blues<-c("lightsteelblue4",  "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")  


####################################################################
### bathymetry plots

dev.new(width=5.26, height=10) 
par(mai=c(0.8,0.9,0.2,0.3)) # , mfrow=c(1,2)
plot(OR_bathy2, image = T, land = T, axes = T, lwd=1, deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy2), "grey"), c(min(OR_bathy2),0,blues)),xlab=expression(paste("Longitude ("^o,'W)')), ylab= expression(paste("Latitude ("^o,'N)')), cex.lab=1.3, cex.axis=1.2)  #, xlim= c(-125.5,-123.5), ylim=c(40.2, 47) 

rect(-126.0909, 40.65895, -123.9263, 46.99730, lwd=2) # 2km domain
rect(-125.1, 42.00001, -123.5, 45.24775, lty=2, lwd=2) # 250m domain

points(-124.52,42.84, pch=21, cex= 1.2, col='black', bg='aquamarine3') # cape blanco
points(-123.997,44.6368, pch=21, cex= 1.2, col='black', bg='aquamarine3') # newport

points(-124.095, 44.6598, pch=22, cex= 1.5, col='black', bg='orange')
points(-124.304, 44.6393, pch=22, cex= 1.5, col='black', bg='orange')
points(-124.067, 44.613, pch=23, col='black', bg='gold', cex=1.2)

text(-124.25, 42.85, "CB", cex=1.2)
text(-123.84, 44.6368, 'N', cex=1.2)

# text(-124.2, 45.75, "CF", col= 'black')
# # text(-124.2, 45.5, "CA", col='darkblue')

# text(-124.25, 45, "CH", col= 'black')
# # text(-124.2, 45.5, "CA", col='darkblue')

# text(-124.25, 44.75, "OR", col= 'black')
# text(-124.35, 44.2, "CP", col= 'black')
# text(-124.79, 42.8, 'RR', col= 'black')

# lines(MR_transformed, col= 'darkorchid4', lwd= 1.4, lty= 1) # plot MR polygons

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,2000,by=500),labcex=0.7,add=T,col='black', lwd= 0.07) 

points(init_pop$lon, init_pop$lat, col='coral3', cex=0.5, pch=20)

##########
# release locations

setwd('/Users/jennifer/Desktop/')
str_name<-"exportImage.tiff"
OR_bathy3<-raster(str_name)

library(graphics)
OR_bathy4<-as.bathy(OR_bathy3)

# init_pop<-pop.roms[[1]][[1]] # only lat and lon values

# dev.new(width=7, height=10) 
# par(mai=c(0.8,0.9,0.2,0.3))

plot(OR_bathy4, image = T, land = T, axes = T, lwd=1, deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy4), "grey"), c(min(OR_bathy4),0,blues)),xlab=expression(paste("Longitude ("^o,'W)')), ylab='', cex.lab=1.3, cex.axis=1.2)  #, xlim= c(-125.5,-123.5), ylim=c(40.2, 47) 

# rect(-126.0909, 40.65895, -123.9263, 46.99730, lwd=2) # 2km domain
# rect(-125.1, 42.00001, -123.5, 45.5, lty=2, lwd=2) # 250m domain

lines(MR_transformed, col= 'darkorchid4', lwd= 1.4, lty= 1) # plot MR polygons

points(init_pop$lon, init_pop$lat, col='brown4', cex=0.8, pch=20)


####################################################################
### Calculations 

setwd('/Users/jennifer/Desktop/model_output_chapter1&2/used_low_res_ROMS')

# load('LPT_data_pp_06.09.2020_JWA_init_2016_low.RData.RData')

#### Depth range (release location)
range(pop.mat$depth[1:250])
range(pop.mat$depth[251:500])
range(pop.mat$depth[501:750])

#### Depth range (entire simulation)
range(surf_par$depth, na.rm=T)
range(mid_par$depth, na.rm=T)
range(deep_par$depth, na.rm=T)

#### % of particles advected out of domain after 10 days
d<-deep_par

tot_len<-250 # total number of particles subset is already contains (each layer has 250 at beginning of simulation)

tenday<-subset(d, del.t==960)
len<-length(unique(tenday$ID)) # number of particles left at 10 days
per_left<-(len/tot_len)*100 # percent of particles left at 10 days
per_out<-100-per_left  # subtract percent left in model from 100 to get % of particles advected out after 10 days
per_out
                                                      
#### % of particles with simulation lengths >=35 days
tmp_use<-tmp4
# par_use<-mid_par

a<-length(which(tmp_use>=23)) 
a
b<-250
b
c<-a/b*100 
c

####################################################################
# Pau 


