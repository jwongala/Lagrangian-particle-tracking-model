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
	
# load("LPT_07.05.2020_JWA_6hour_2018_high.RData")

# pop.fin.a16.lay1<-read.csv("JWA_20200123_April2016_depth01.csv")
# dim(pop.fin) # 4362240      28
# head(pop.fin)


####################################################################
### create dataframe to plot with 

out.row<-nrow(pop.roms[[1]][[1]])
out.col<-1+ncol(pop.roms[[1]][[1]])
pop.fish<-as.data.frame(matrix(1:(out.col*out.row),ncol=out.col,nrow=out.row)*NA)
pop.roms2<-list()
pop.roms2[[1]]<-pop.fish
pop.mat<-NULL
pop.mat2<-NULL # tides 
pop.roms4<-pop.roms[-1] # use this for data after 1st day () tides


# for(l in 1:(length(pop.roms))){ 
	# tmp<-pop.roms[[l]]
	# for(j in 1:length(tmp)){
	# tmp[[j]]$del.t<-j+length(tmp)*(l-1)
	# pop.fish<-tmp[[j]]
	# pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
# }}

l=1
tmp<-pop.roms[[l]]

for(j in 1:length(tmp)){
	tmp[[j]]$del.t<-j+length(tmp)*(l-1)
	pop.fish<-tmp[[j]]
	pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
}

pop.roms3<-pop.roms2 # rename tide data

rm(pop.roms2)

pop.fish<-as.data.frame(matrix(1:(out.col*out.row),ncol=out.col,nrow=out.row)*NA)
pop.roms2<-list()
pop.roms2[[1]]<-pop.fish


# tides data 
for(l in 1:(length(pop.roms4))){ 
	tmp<-pop.roms4[[l]]
	for(j in 1:length(tmp)){
	tmp[[j]]$del.t<-j+length(tmp)*(l-1)
	pop.fish<-tmp[[j]]
	pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
}}

# 
for(m in 1:(length(pop.roms3))){
	print(m)
	nw<-pop.roms3[[m]]
	nw$timestep<-m
	pop.mat2<-rbind(pop.mat2, nw)
}

for(m in 1:(length(pop.roms2))){
	print(m)
	nw<-pop.roms2[[m]]
	nw$timestep<-m
	pop.mat<-rbind(pop.mat, nw)
}

pop.mat<-rbind(pop.mat, pop.mat2) # combine day 1 ROMS and rest of days of data

dim(pop.mat) # 2880   16
head(pop.mat)
tail(pop.mat)
                            
# setwd('/Users/jennifer/Desktop')
# save(pop.mat, file='LPT_data_05.09.2020_6hour_high_01.RData')


# pop.fin.a16.lay40.hires<-pop.mat
# head(pop.fin.a16.lay40.hires)
# tail(pop.fin.a16.lay40)
# dim(pop.fin.a16.lay40) # 92943615       16


# ### date split date hour to get individual year, month, day, and hour
# datehour_levels<-levels(pop.mat$datehour)

# datehoursplt2<-strsplit(as.character(pop.mat$datehour), "-")

# year<-length(datehoursplt2)*NA
# month<-length(datehoursplt2)*NA
# dayhour<-length(datehoursplt2)*NA
# day<-length(datehoursplt2)*NA
# hour<-length(datehoursplt2)*NA

# for (z in 1:length(datehoursplt2)){
	# year[z]<-datehoursplt2[[z]][1]
	# month[z]<-datehoursplt2[[z]][2]
	# dayhour[z]<-datehoursplt2[[z]][3]
# }

# dayhoursplt<-strsplit(dayhour, " ")
	
# for (v in 1:length(dayhour)){
	# day[v]<-dayhoursplt[[v]][1]
	# hour[v]<-dayhoursplt[[v]][2]
# }


# ### convert u and v from degrees back to m/s and then convert units (m/s) to (km/day)
# m2km<-1000
# secday<-1.15741e-5
# Rad<-6378137 # radius of the Earth in meters

# v_ms<-(Rad*pi*pop.mat$v_vel[pop.mat$ID])/180
# v_kmday<-v_ms*((secday)/m2km)

# u_ms<-((pop.mat$u_vel[pop.mat$ID]*Rad*pi)+cos(pi*v_ms))/180
# u_kmday<-u_ms*((secday)/m2km)

# curnt_vel<-sqrt((u_kmday^2)+(v_kmday^2))

# # column bind all new variables curnt_vel, u_ms, v_ms, u_kmday, and v_kmday to pop.fin
# pop.fin<-cbind(pop.mat, u_ms, u_kmday, v_ms, v_kmday, curnt_vel) # day, year, month,  hour,
# pop.fin$timestep<-pop.fin$del.t*15
# pop.fin$timestep_hour<-pop.fin$timestep/60
# pop.fin$timestep_day<-pop.fin$timestep_hour/24
# pop.fin$timestep_2hr<-pop.fin$timestep_hour*2

# # write.csv(pop.fin, "JWA_20200123_April2016_depth01.csv")

# head(pop.fin)


####################################################################
### load files 

### Territorial Sea outline

TS<-readOGR(dsn=path.expand('/Users/wongalaj/Desktop/BASE_Oregon_Territorial_Sea_3NMLine_OCMP'), layer= 'BASE_Oregon_Territorial_Sea_3NMLine_OCMP')

class(TS) # SpatialLinesDataFrame
crs(TS) # this tells us how it is projected so that I can convert it below

extent(TS) # this shows what the extent in based on the units it is using

TS_transformed<-spTransform(TS,CRS=crsmerc) # this transforms TS so that is is now projected with lat and lon coordinates 

### Marine Reserve Outlines

setwd("/Users/jennifer/Desktop/")

MR<-readOGR(dsn=path.expand("/Users/jennifer/Desktop/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"), layer= 'MPA_MR_COMP_Boundaries_UTM10N')

class(MR) # SpatialPolygonsDataFrame
crs(MR)

extent(MR)

crsmerc=CRS("+proj=longlat +lat_1=43 +lat_2=48 +lat_0=41 +lon_0=-117 +x_0=700000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ") # this transforms the model 

MR_transformed<-spTransform(MR, CRS=crsmerc)


### Bathymetry
bathy.dat<-read.table('etopo1_bedrock.xyz', sep='')
names(bathy.dat)<-c('lon', "lat", 'depth')

bathy.dat$depth[bathy.dat$depth>0]<-NA # Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]



####################################################################
### Post processing of data (removing trajectories that go past domain and return)

## Lat and lon boundaries

lat1<-42 
lat2<-45.24775


lon1<-(-125.10) # -126.0909
lon2<-(-123.9263)  

# setwd("/Users/wongalaj/Desktop/Model/Data")  

## Load Data
may2018<-pop.mat # pop.fin.a16.lay40.hires # read.csv('ORcoast_20191010_May2016_thesis.csv')  
dim(may2018) # 4362240      28

### remove the inital positions for all particles
# may2018<-subset(may2018, del.t!=1) # 4360820      28
# dim(may2018) # 92943260       24

## unique ID numbes that have lon values less than lon1 
lon.uni<-sort(unique(may2018$ID[which(may2018$lon<lon1)]))
length(lon.uni) # 4
# lon.uni<-levels(lon.uni) # uncomment this when I there are multiple particles being released from the same location

lat.uni2<-sort(unique(may2018$ID[which(may2018$lat>=lat2)])) 
length(lat.uni2) # 68
lat.uni1<-sort(unique(may2018$ID[which(may2018$lat<=lat1)])) 
length(lat.uni1) # 700

lat.uni<-sort(c(lat.uni1,lat.uni2)) # levels(lat.uni1) # only included lat.uni beceause both include all of the particles released, particles are the same for both 
length(lat.uni) # 23958

### these are the ID numbers that have either lat1, lat2, or lon1 as a value
#### now I need to subset out each df for all of the individual ID numbers from uni_ids and remove the trajectory once lat or lon is equal to a boundary condition value
uni_ids<-sort(c(lat.uni, lon.uni)) # only included lat.uni because both lat.uni and lon.uni include all particle IDs
length(uni_ids) # 24,024
uni_ids<-as.numeric(uni_ids)
### new dataframe without particle IDs that need to be edited to remove trajectory that touched domain border

## Commented this out because all particles have part of trajectory outside of model domain. 
maynew<-may2018[!may2018$ID %in% uni_ids,] 
dim(maynew) # 39301769       24

may2018$ID<-as.numeric(may2018$ID)

par_fill<-NULL
# i=i+1
for(i in 1:length(uni_ids)){
	
	print(i)
	
	# subset out particles with specific particle ID
	par<-subset(may2018, ID==uni_ids[i]) 
	dim(par)
	# need to find the row number where the particle touched the boundary
	num_remove<-which(par$lat<=lat1 | par$lat>=lat2 | par$lon<lon1)[1] 
	
	# go to df and remove the part of the trajectory that touched the boundary
	par_new<-par[-(num_remove:nrow(par)),]
	dim(par_new)
	
	# rbind new dataset to par_fill	
	par_fill<-rbind(par_fill, par_new)	
	
}

dim(par_fill)
head(par_fill)

### Remove NAs from maynew use this data set to plot
## has NA's removed
## removed particles with trajectories that leave model domain (need to fix this once I have another data set to use)
## removed initial location
# maynew2<-maynew[,-(17:20)]
# dim(maynew) # 39301769  20
# head(maynew)

# maynew2<-na.omit(maynew2)
# dim(maynew2)

# # ### plot to check that no more particles are touching the borders
# plot(may2018$lon, may2018$lat, pch= '.', col= 'green')
# points(maynew2$lon, maynew2$lat, pch= '.', col= 'black')

uniID_may18<-length(unique(maynew$ID)) # 28213

### Combine maynew2 df and par_fill df
# par_fill2<-par_fill[,-(17:20)] # remove col 17-20 from par_fill

data_use<-rbind(maynew, par_fill)
dim(data_use) # 179386     16


save(data_use, file='LPT_data_pp_07.05.2020_JWA_6hour_2018_high.RData')

# plot comparing using AKv and AKt
# plot(data_use$del.t, data_use$lat, col=data_use$ID, ylab='latitude', xlab='model time step', main='AKt')
# plot(data_use$del.t, data_use$lat, col=data_use$ID, ylab='latitude', xlab='model time step', main='AKv')

## subset out out for 44N and 45N

# surf_par<-seq(1,250)
# seq_45<-seq(1,300, by=1)

# seq_44<-seq(301,600, by=1)

# seq_45mid<-seq(601,900,by=1)

# seq_44mid<-seq(901,1200, by=1)

# seq_45deep<-seq(1201,1500,by=1)

# seq_44deep<-seq(1501,1800,by=1)

## test plots at 44N and 45N

# # surfpar_44<-subset(data_use, ID %in% seq_44)
# plot(surfpar_44$del.t, surfpar_44$lat)

# surfpar_45<-subset(data_use, ID %in% seq_45)
# plot(surfpar_45$del.t, surfpar_45$lat)

# midpar_44<-subset(data_use, ID %in% seq_44mid)
# plot(midpar_44$del.t, midpar_44$lat)

# midpar_45<-subset(data_use, ID %in% seq_45mid)
# plot(midpar_45$del.t, midpar_45$lat)

# deeppar_44<-subset(data_use, ID %in% seq_44deep)
# plot(deeppar_44$del.t, deeppar_44$lat)

# deeppar_45<-subset(data_use, ID %in% seq_45deep)
# plot(deeppar_45$del.t, deeppar_45$lat)

# d16_1<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-03-31 00:00:00"), length.out=nrow(surfpar_44), format='%d %b')

# d16_2<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(surfpar_45))

# d16_3<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(midpar_44))

# d16_4<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(midpar_45))

# d16_5<-as.Date(seq(as.POSIXct("2016-03-27"), as.POSIXct("2016-05-01"), length.out=nrow(deeppar_44)))

# d16_6<-as.Date(seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(deeppar_45)), format='%b %d')
# # , format='%b %d'

# text(-124.34,42.85,"CB", cex=1.5)

# par(mfrow=c(3,2), mai=c(0.2, 0.6, 0.5, 0.1))
# plot(d16_1, 1:length(d16_1)) # , 1:length(d16_1), 'n')

# plot(surfpar_44$del.t, surfpar_44$lat, pch='.', main='44N particle comparison', col='black', xlab='', ylab='', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(42,45.24770), xaxt='n')
# text(x=450, y=45.24770, 'surface', col='black', cex=2)

# plot(surfpar_45$del.t, surfpar_45$lat, pch='.', main='45N particle comparison', col='black', xlab='', ylab='', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(42,45.24770), xaxt='n')

# plot(midpar_44$del.t, midpar_44$lat, pch='.', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(42,45.24770), ylab='Latitude', col='blue', xaxt='n', xlab='')
# text(x=450, y=45.24770, 'middle', col='blue', cex=2)

# plot(midpar_45$del.t, midpar_45$lat, pch='.', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(42,45.24770), col='blue', ylab='', xaxt='n', xlab='')

# plot(deeppar_44$del.t, deeppar_44$lat, pch='.', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(42,45.24770), col='red', ylab='', xaxt='n', xlab='')
# text(x=450, y=45.2, 'deep', col='red', cex=2)
# # axis(1, at=deeppar_44$del.t, labels=d16_5)

# plot(deeppar_45$del.t, deeppar_45$lat, pch='.',  cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(42,45.24770), col='red', ylab='', xlab='', xaxt='n') # 


# #### calculations 
# tenday<-subset(deeppar_45, del.t==480)
# dim(tenday)

# avg<-mean(tenday$lat) # 43.07116
# avg

# sd_10<-sd(tenday$lat)

# err<-qnorm(0.95)*sd_10/sqrt(length(unique(tenday$ID)))
# err

# left<-avg-err

# right<-avg +err

# # 44.71533 44.74026
# n<-300

# samp_avg2<-NA*(1:100)

# for(i in 1:100){
	# samp<-sample(tenday $ID, size=n, replace=F)

	# datsub<-subset(tenday, ID %in% samp)

	# sampavg<-mean(datsub$lat)
	# samp_avg2[i]<-sampavg
# }

# (1:100)[samp_avg2>left & samp_avg2<right]

# points(deeppar_45$lat)
# points(deeppar_45$del.t, deeppar_45$lat)
# plot(as.factor(deeppar_45$datehour), 1:nrow(deeppar_45), type='n')
# axis(1,at=deeppar_45$del.t, d16_6)
# axis(side=1, at=deeppar_45$del.t, labels= d16_6)

# plot(d16_6,1:length(d16_6), type="n", xlab="", ylab="", cex.lab=2, cex.axis=1.5, yaxt='n') # plotting to show date on axis

# plot(deeppar_45$lat~as.Date(deeppar_45$datehour), pch='.', ylim=c(42,45.24770))
# day<-deeppar_45$

### 
surf_par<-subset(data_use, ID %in% seq(1,250,1))
dim(surf_par) # 179386     16

mid_par<-subset(data_use, ID %in% seq(251,500,1))
dim(mid_par) # 179386     16

deep_par<-subset(data_use, ID %in% seq(501,750,1))
dim(deep_par) # 179386     16
## surface comparison plot
# par(mfrow=c(1,3), mai=c(0.7, 0.8, 0.8, 0.2))

# text(-124.34,42.85,"CB", cex=1.5)

# # plot(surf_par1$del.t, surf_par1$depth, pch='.', main='surface particle comparison', col='black', xlab='', ylab='depth (m)', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(-100,0)) # latitude depth (m)
# points(surf_par2$del.t, surf_par2$depth, pch='.', col='red')
# legend(x='bottomleft', y=NULL, legend=c('group 1', 'group 2'), col=c('black', 'red'), pch=20, bty='n', cex=1.5)

# ## mid depth comparison plot
# plot(mid_par1$del.t, mid_par1$depth, pch='.', main='mid depth particle comparison', col='black', xlab='model time step', ylab='', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(-100,0)) # latitude depth (m)
# points(mid_par2$del.t, mid_par2$depth, pch='.', col='green')
# legend(x='bottomleft', y=NULL, legend=c('group 1', 'group 2'), col=c('black', 'green'), pch=20, bty='n', cex=1.5)

# ## deep depth
# plot(deep_par1$del.t, deep_par1$depth, pch='.', main='deep depth particle comparison', col='black', xlab='', ylab='', cex.main=2, cex.lab=2, cex.axis=1.5, ylim=c(-100,0)) # latitude depth (m)
# points(deep_par2$del.t, deep_par2$depth, pch='.', col='blue')
# legend(x='bottomleft', y=NULL, legend=c('group 1', 'group 2'), col=c('black', 'blue'), pch=20, bty='n', cex=1.5)



# mid_par<-subset(data_use, ID %in% seq(101,200,1))
# dim(mid_par) # 177727     16

# deep_par<-subset(data_use, ID %in% seq(201,300,1))
# dim(deep_par) # 180791     16


# par_500<-subset(surf_par, del.t==500)
# dim(par_500)

# id_500<-unique(par_500$ID)

# par_45<-subset(pop.mat, ID==45)
# head(par_45)
# tail(par_45)

# plot(par_45$del.t, par_45$lat)

## test plots to send to Lorenzo
# par2<-subset(data_use, ID==2)
# dim(par2)



# plot(data_use$del.t, data_use$lat, col=data_use$ID, pch='.', main= "2km ROMS - April 2016 - Initial release", xlab="Time",ylab=expression(paste("Latitude ("^o,'N)')))

# plot(data_use$del.t, data_use$depth, pch=19, col=data_use$ID, cex=0.5, main= "2km ROMS - April 2016 - Initial release", xlab="Time",ylab="Depth (m)")
setwd('/Users/jennifer/Desktop/model_output_chapter1&2/used_high_res_ROMS')
load("LPT_data_pp_06.11.2020_JWA_init_2016_high.RData")


#53A567FF #56A8CBFF #DA291CFF
#### lat
plot(surf_par$del.t, surf_par$lat, col="#53A567FF", pch='.', main= "2km ROMS - April 2016 - initial release", xlab="Time (model time step)",ylab=expression(paste("Latitude ("^o,'N)')))

points(mid_par$del.t, mid_par$lat, col='#56A8CBFF', pch='.')

points(deep_par$del.t, deep_par$lat, col='#DA291CFF', pch='.')

# legend(x='topright', y=NULL, legend=c('surface', 'middle', 'deep'), bty="n", lty=1, col=c('#53A567FF', '#56A8CBFF', '#DA291CFF'), cex=0.8)

#### depth
plot(surf_par$del.t, surf_par$depth, col="#53A567FF", pch='.', main= "2km ROMS - April 2016 - 6 hour release", xlab="Time (model time step)",ylab='Depth (m)', ylim=c(-800,0))

points(mid_par$del.t, mid_par$depth, col='#56A8CBFF', pch='.')

points(deep_par$del.t, deep_par$depth, col='#DA291CFF', pch='.')

# legend('bottomleft', y=NULL, legend=c('surface', 'middle', 'deep'), bty="n", lty=1, lwd=3, col=c('#53A567FF', '#56A8CBFF', '#DA291CFF'))



####################################################################
### Load marmap bathy locations

# ROMS model domain + more land so labels can fit
# OR_bathyOG<-getNOAA.bathy(lon1= -125.1, lon2= -123.5, lat1= 45.5, lat2= 42.00001, resolution=1) # get bathymetry from marmap (NOAA data) (40.65895 - 47.01520)

# ROMS model domain
# OR_bathy<-getNOAA.bathy(lon1= -126.0909, lon2= -123.9263, lat1= 46.99730, lat2= 40.65895, resolution=1) # get bathymetry from marmap (NOAA data)

# OR_bathy<-getNOAA.bathy(lon1= -127, lon2= -123.5, lat1= 47.5, lat2= 40, resolution=1) # get bathymetry from marmap (NOAA data)

# zoomed in version 
# OR_bathy2<-getNOAA.bathy(lon1= -125, lon2= -124.3, lat1= 42.7, lat2= 42.1, resolution=1) # get bathymetry from marmap (NOAA data)

# lat_seq<-c(42.00001, 45.5,45.5, 42.00001)

# lon_seq<-c(125.1,125.1,123.5,123.5)*(-1)



library(tiff)
setwd('/Users/jennifer/Desktop/')
str_name<-"exportImage.tiff"
OR_bathy<-raster(str_name)

library(graphics)
OR_bathy2<-as.bathy(OR_bathy)

blues<-c("lightsteelblue4",  "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")

dev.new(width=5, height=8.5) 
par(mai=c(0.8,0.9,0.2,0.3))
plot(OR_bathy2, image = T, land = T, axes = T, lwd=1, deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy2), "grey"), c(min(OR_bathy2),0,blues)),xlab=expression(paste("Longitude ("^o,'W)')), ylab= expression(paste("Latitude ("^o,'N)')), cex.lab=1.3, cex.axis=1.2)  # , xlim= c(-125,-124)


text(-124.44,42.85,"CB", cex=1.5)
text(-123.97, 44.6368, 'N', cex=1.5)

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
 
# # segments(-125.15,42.5,-124.42,42.5) 
# segments(-125.15,43,-124.45,43) # make lines to put sections apart
# segments(-125.15,43.5,-124.25,43.5)
# segments(-125.15,44,-124.14,44)
# segments(-125.15,44.5,-124.1,44.5)
# segments(-125.15,45,-124.01,45)


# # text("A", x=-125.06, y=45.06, cex=2)
# text("B", x=-125.06, y=44.56, cex=2)
# text("C", x=-125.06, y=44.06, cex=2)
# text("D", x=-125.06, y=43.56, cex=2)
# text("E", x=-125.06, y=43.06, cex=2)
# text("F", x=-125.06, y=42.56, cex=2)
# text("G", x=-125.06, y=42.06, cex=2)



# rect(-126.0909, 40.65895, -123.9263, 46.99730) # 2km domain

# rect(-125.1, 42.00001, -123.5, 45.5, lty=2) # 2km domain

# points(-124.52,42.84, pch=21, cex= 1.5, col='black', bg='green') # cape blanco

# points(-123.997,44.6368, pch=21, cex= 1.5, col='black', bg='green') # newport

# points(-124.095, 44.6598, pch=22, cex= 1.5, col='black', bg='orange')
# points(-124.304, 44.6393, pch=22, cex= 1.5, col='black', bg='orange')



###############################################################
### Select out particles from each grid cell along the 10m bathy contour

setwd('/Users/jennifer/Desktop/')
str_name<-"OR_bathy_marmap_small.tiff"
OR_bathy<-raster(str_name)

library(graphics)
OR_bathy2<-as.bathy(OR_bathy)

blues<-c("lightsteelblue4",  "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")

# take subset of surf data

#53A567FF #56A8CBFF #DA291CFF
surf_samp<-sample(unique(surf_par$ID), size=10, replace=F)
surf_par_sub<-subset(surf_par, ID %in% surf_samp)
unique(surf_par_sub$ID)

mid_samp<-sample(unique(mid_par$ID), size=10, replace=F)
mid_par_sub<-subset(mid_par, ID %in% mid_samp)
unique(mid_par_sub$ID)

deep_samp<-sample(unique(deep_par$ID), size=10, replace=F)
deep_par_sub<-subset(deep_par, ID %in% deep_samp)
unique(deep_par_sub$ID)

# # star_surf<-subset(surf_par_sub, del.t==1)
# star_mid<-subset(mid_par_sub, del.t==1)
# star_deep<-subset(deep_par_sub, del.t==1)

star_surf<-subset(surf_par, del.t==1)
star_mid<-subset(mid_par, del.t==1)
star_deep<-subset(deep_par, del.t==1)

# par25<-subset(data_use, ID==45)
# start<-subset(par25, del.t==1)

dev.new(width=5.12, height=8.9)
par(mai=c(0.9,0.9,0.25,0.3))
plot(OR_bathy2, image = T, land = T, axes = T, lwd=1, deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy2), "grey"), c(min(OR_bathy2),0,blues)),xlab=expression(paste("Longitude ("^o,'W)')), ylab= expression(paste("Latitude ("^o,'N)')), cex.lab=1.3, cex.axis=1.2, main='250m - 2016 - 6 hour release - deep')  # , xlim= c(-125,-124)

# points(par25$lon, par25$lat)
# points(start$lon, start$lat, col="red")
# points(data_use$lon, data_use$lat, col= 'black', pch='.', cex= 0.3)

# points(surf_par$lon, surf_par$lat, col= '#53A567FF', pch='.', cex= 0.3)
# points(mid_par$lon, mid_par$lat, col= '#56A8CBFF', pch='.', cex= 0.3)
points(deep_par$lon, deep_par$lat, col= '#DA291CFF', pch='.', cex= 0.3)

# points(surf_par_sub$lon, surf_par_sub$lat, col= '#53A567FF', pch=19, cex= 0.3)
# points(mid_par_sub$lon, mid_par_sub$lat, col= "#56A8CBFF", pch=19, cex= 0.3)
# points(deep_par_sub$lon, deep_par_sub$lat, col= "#DA291CFF", pch=19, cex= 0.3)

# points(star_surf$lon, star_surf$lat,col='black', bg='#53A567FF', pch=21)
# points(star_mid $lon, star_mid $lat, col='black', bg='#56A8CBFF', pch=21)
points(star_deep$lon, star_deep $lat, col='black', bg='#DA291CFF', pch=21)


# lines(MR_transformed, col= 'black', lwd= 2, lty= 1) # plot MR polygons
  
text(-124.25,42.85,"Cape Blanco", cex=1.2)
text(-123.85, 44.6368, 'Newport', cex=1.2)

# legend('bottomright', y=NULL, legend=c('surface', 'middle', 'deep'), bty="n", lty=1, col=c('#53A567FF', '#56A8CBFF', '#DA291CFF'), lwd=3)

####################################################################
### Histogram of time it takes particle to be advected out of model
# surface= '#53A567FF' mid= "#56A8CBFF" deep= "#DA291CFF"

# setwd('/Users/jennifer/Desktop/model_output_chapter1&2/used_low_res_ROMS/2016')

# load('LPT_data_pp_06.10.2020_JWA_6hour_2016_low.RData')
# load('LPT_data_pp_05.20.2020_6hour_low_01.RData')

setwd('/Users/jennifer/Desktop/model_output_chapter1&2/used_high_res_ROMS')
load('LPT_data_pp_06.13.2020_JWA_6hour_2016_high.RData')
# load('LPT_data_pp_05.10.2020_init_high_01.RData')
# load('LPT_data_pp_05.09.2020_6hour_high_01.RData')

# lowinit<-data_use
# lowsix<-data_use

# hiinit<-data_use
# hisix<-data_use

dim(data_use)
par_ID<-unique(data_use $ID)

tmp<-NA*1:length(par_ID)

surf_par<-subset(data_use, ID %in% seq(1,250,1))
dim(surf_par) # 179386     16
par_ID2<-unique(surf_par$ID)
tmp2<-NA*1:length(par_ID2)

mid_par<-subset(data_use, ID %in% seq(251,500,1))
dim(mid_par) # 177727     16
par_ID3<-unique(mid_par$ID)
tmp3<-NA*1:length(par_ID3)

deep_par<-subset(data_use, ID %in% seq(501,750,1))
dim(deep_par) # 180791     16
par_ID4<-unique(deep_par$ID)
tmp4<-NA*1:length(par_ID4)


# find out what is the max del.t that the particle left the model at? 
# for(i in 1:length(par_ID)){
	
	# print(i)
	# tmp[i]<-max(lowsix $del.t[which(lowsix $ID %in% par_ID[i])])
	
# }

for(i in 1:length(par_ID2)){
	
	print(i)
	tmp2[i]<-max(surf_par$del.t[which(surf_par$ID %in% par_ID2[i])])
	
}

for(i in 1:length(par_ID3)){
	
	print(i)
	tmp3[i]<-max(mid_par$del.t[which(mid_par$ID %in% par_ID3[i])])
	
}

for(i in 1:length(par_ID4)){
	
	print(i)
	tmp4[i]<-max(deep_par$del.t[which(deep_par$ID %in% par_ID4[i])])
	
}

tmp_tot<-tmp/96
max(tmp_tot)

tmp2<-tmp2/96
tmp3<-tmp3/96
tmp4<-tmp4/96

tmp2<-c(tmp2,NA)
# tmp3<-c(tmp3,NA)
tmp4<-c(tmp4,NA)

df<-cbind(tmp2, tmp3, tmp4)

# adjustcolor('#53A567FF', alpha.f=0.98) # "#53A56733" green
# adjustcolor('#56A8CBFF', alpha.f=0.2) # "#56A8CBB3" blue
# adjustcolor('#DA291CFF', alpha.f=0.1) # "#DA291CB3" red

setwd('/Users/jennifer/Desktop')

init_low_use<-read.csv('lowinit_displengths.csv')
# # dim(init_low_use)

six_low_use<-read.csv('lowsix_displengths.csv')
# head(six_low_use)

high_init_use<-read.csv('hiinit_displengths.csv')

high_six_use<-read.csv('hisix_displengths.csv')


### stacked barplot

ggplot(data = init_low_use, aes(x = disp_length, fill = release)) + geom_histogram(colour = 'black', binwidth=1) + labs(title='2km initial release', y='counts', x='Dispersal lengths (Days)') + theme(text=element_text(size=15)) + ylim(c(0,85)) + xlim(c(0,35))

dev.new()
ggplot(data = six_low_use, aes(x = disp_length, fill = release)) + geom_histogram(colour = 'black', binwidth=1) + labs(title='2km 6 hour release', y='counts', x='Dispersal lengths (Days)') + theme(text=element_text(size=15)) + ylim(c(0,85)) + xlim(c(0,35))

dev.new()
ggplot(data = high_init_use, aes(x = disp_length, fill = release)) + geom_histogram(colour = 'black', binwidth=1) + labs(title='250m initial release', y='counts', x='Dispersal lengths (Days)') + theme(text=element_text(size=15)) + ylim(c(0,85)) + xlim(c(0,35))

dev.new()
ggplot(data = high_six_use, aes(x = disp_length, fill = release)) + geom_histogram(colour = 'black', binwidth=1) + labs(title='250m 6 hour release', y='counts', x='Dispersal lengths (Days)') + theme(text=element_text(size=15)) + ylim(c(0,85)) + xlim(c(0,35))


####################################################################
### Calculations 
setwd('/Users/jennifer/Desktop/model_output_chapter1&2/used_low_res_ROMS')

load('LPT_data_pp_06.09.2020_JWA_init_2016_low.RData.RData')

#### Depth range (release location)
range(pop.mat$depth[1:250])
range(pop.mat$depth[251:500])
range(pop.mat$depth[501:750])

#### Depth range (entire simulation)
range(surf_par$depth)
range(mid_par$depth)
range(deep_par$depth)

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

a<-length(which(tmp_use>=35)) 
a
b<-250
b
c<-a/b*100 
c


####################################################################
### Convert km size of MR to degrees so that I can divide up coastal ocean into grid cells that are no bigger than the largest size of the Otter Rock MR

Rad<-6378137 # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method

OR.size<-3 # km2; size of Otter Rock MR

deg.lat<-6.33835 # length of coastline measured in degrees 

deg.lon<-2.1646

## conversion of coastline length in degrees to km

lat.len<-Rad*pi*(deg.lat/180)*0.001 # length of coastline (latitude in km)

lon.len<-(Rad*pi*(cos(pi*(lat.len/180)))*(deg.lon/180))*0.001

lon.state<-5.556 # km (width of state waters) 

grid.size.lat<-lat.len/OR.size # number of grid cells latitudinally

grid.size.lon<-lon.len/lon.state # number of grid cells longitudinally

####################################################################
### Subset out same particles as high res ROMS 

# # release_par2<-seq(360, 1440, 360)

# lon_release2<-rep(-123.9263, length.out= 100)

# lat_release2<-seq(40.65895, 47.01520, length.out= 100)

# layer_release<-40

### use makePop.R function to make a pop that will give initial locations for high res model. And then subset that info out (only do this for LOW RES MODEL RUNS)

# # init.lon.hr<-pop[[1]]$lon
# length(init.lon.hr)

# init.lat.hr<-pop[[1]]$lat
# length(init.lat.hr)

# init_datause<-subset(data_use, del.t==1)
# unique(init_datause$ID)

# for(i in 1:nrow(init_datause)){
	
	# print(i)
	# dist<-distance.function(init_datause$lat[i], init_datause$lon[i], init.lat.hr, init.lon.hr)  
	# tmp<-1*(dist==min(dist))
	# print(paste0("ID=", init_datause$ID[i]))
	# print(paste0("minimum distance in meters ", min(dist)))
	# print(paste0('closest to HR initial location  ', which(tmp==1)))
	
# }

# # ### convert ID in data_use to numeric
# data_use$ID<-as.numeric(data_use$ID)
# class(data_use$ID)

# IDs_subset<-c(5,7,12,15,19,22,26,30,33,37,40,44,47,51,55,58,62,65,69,73,76,80,83,87,90,94,98,101,105,108,112,115,119,123,126,130,134,137,140,144,148,152,155,158,162,166,169,173,176,180,183,187,191,194,198,201,205,208,212,216,219,223,226,230,235,237,241,244,248,251,255,258,262,266,269,273,276,280,283,287,291,294,298,301,305,308,312,316,319,323,326,330,333,337, 341, 344, 348, 351,353)

# data_init<-subset(data_use, ID %in% IDs_subset)
# dim(data_init) # 192651     15

# IDs_360<-c(5.360,7.360,12.360,15.360,19.360,22.360,26.360,30.360,33.360,37.360,40.360,44.360,47.360,51.360,55.360,58.360,62.360,65.360,69.360,73.360,76.360,80.360,83.360,87.360,90.360,94.360,98.360,101.360,105.360,108.360,112.360,115.360,119.360,123.360,126.360,130.360,134.360,137.360,140.360,144.360,148.360,152.360,155.360,158.360,162.360,166.360,169.360,173.360,176.360,180.360,183.360,187.360,191.360,194.360,198.360,201.360,205.360,208.360,212.360,216.360,219.360,223.360,226.360,230.360,235.360,237.360,241.360,244.360,248.360,251.360,255.360,258.360,262.360,266.360,269.360,273.360,276.360,280.360,283.360,287.360,291.360,294.360,298.360,301.360,305.360,308.360,312.360,316.360,319.360,323.360,326.360,330.360,333.360,337.360, 341.360, 344.360, 348.360, 351.360,353.360)


# data_360<-subset(data_use, ID %in% IDs_360)
# dim(data_360) # 450376     15


# IDs_720<-c(5.720,7.720,12.720,15.720,19.720,22.720,26.720,30.720,33.720,37.720,40.720,44.720,47.720,51.720,55.720,58.720,62.720,65.720,69.720,73.720,76.720,80.720,83.720,87.720,90.720,94.720,98.720,101.720,105.720,108.720,112.720,115.720,119.720,123.720,126.720,130.720,134.720,137.720,140.720,144.720,148.720,152.720,155.720,158.720,162.720,166.720,169.720,173.720,176.720,180.720,183.720,187.720,191.720,194.720,198.720,201.720,205.720,208.720,212.720,216.720,219.720,223.720,226.720,230.720,235.720,237.720,241.720,244.720,248.720,251.720,255.720,258.720,262.720,266.720,269.720,273.720,276.720,280.720,283.720,287.720,291.720,294.720,298.720,301.720,305.720,308.720,312.720,316.720,319.720,323.720,326.720,330.720,333.720,337.720, 341.720, 344.720, 348.720, 351.720,353.720)

# data_720<-subset(data_use, ID %in% IDs_720)
# dim(data_720) # 341076     15



# IDs_1080<-c(5.1080,7.1080,12.1080,15.1080,19.1080,22.1080,26.1080,30.1080,33.1080,37.1080,40.1080,44.1080,47.1080,51.1080,55.1080,58.1080,62.1080,65.1080,69.1080,73.1080,76.1080,80.1080,83.1080,87.1080,90.1080,94.1080,98.1080,101.1080,105.1080,108.1080,112.1080,115.1080,119.1080,123.1080,126.1080,130.1080,134.1080,137.1080,140.1080,144.1080,148.1080,152.1080,155.1080,158.1080,162.1080,166.1080,169.1080,173.1080,176.1080,180.1080,183.1080,187.1080,191.1080,194.1080,198.1080,201.1080,205.1080,208.1080,212.1080,216.1080,219.1080,223.1080,226.1080,230.1080,235.1080,237.1080,241.1080,244.1080,248.1080,251.1080,255.1080,258.1080,262.1080,266.1080,269.1080,273.1080,276.1080,280.1080,283.1080,287.1080,291.1080,294.1080,298.1080,301.1080,305.1080,308.1080,312.1080,316.1080,319.1080,323.1080,326.1080,330.1080,333.1080,337.1080, 341.1080, 344.1080, 348.1080, 351.1080,353.1080)

# data_1080<-subset(data_use, ID %in% IDs_1080)
# dim(data_1080) # 352280     15


# IDs_1440<-c(5.1440,7.1440,12.1440,15.1440,19.1440,22.1440,26.1440,30.1440,33.1440,37.1440,40.1440,44.1440,47.1440,51.1440,55.1440,58.1440,62.1440,65.1440,69.1440,73.1440,76.1440,80.1440,83.1440,87.1440,90.1440,94.1440,98.1440,101.1440,105.1440,108.1440,112.1440,115.1440,119.1440,123.1440,126.1440,130.1440,134.1440,137.1440,140.1440,144.1440,148.1440,152.1440,155.1440,158.1440,162.1440,166.1440,169.1440,173.1440,176.1440,180.1440,183.1440,187.1440,191.1440,194.1440,198.1440,201.1440,205.1440,208.1440,212.1440,216.1440,219.1440,223.1440,226.1440,230.1440,235.1440,237.1440,241.1440,244.1440,248.1440,251.1440,255.1440,258.1440,262.1440,266.1440,269.1440,273.1440,276.1440,280.1440,283.1440,287.1440,291.1440,294.1440,298.1440,301.1440,305.1440,308.1440,312.1440,316.1440,319.1440,323.1440,326.1440,330.1440,333.1440,337.1440, 341.1440, 344.1440, 348.1440, 351.1440,353.1440)

# data_1440<-subset(data_use, ID %in% IDs_1440)
# dim(data_1440) # 327458     15

####################################################################
### Separate datasets based on release time 


hour6<-360 # 360 mins
stop_par_release<-20160 # 20160 is correct unit for 2 weeks in minutes

release_par<-seq(from= hour6, to= stop_par_release, by= hour6)
head(release_par)

data_use$ID<-as.numeric(data_use$ID) # convert ID to numeric

### subset out particles with ID tags from 360 - 20160

ID_init<-seq(1,100,1) # seq(1,355,1)
ID_360<-seq(1.3600, 100.3600, 1) # seq(1.360, 355.360, 1)
ID_720<-seq(1.1155, 100.1155, 1) # seq(1.720, 355.720, 1)
ID_1080<-seq(1.3750, 100.3750,1) # seq(1.1080, 355.1080, 1)
ID_1440<-seq(1.7650, 100.7650,1) # seq(1.1440, 355.1440, 1)


# ID_1800<-seq(1.1800, 355.1800, 1)
# ID_2160<-seq(1.2160, 355.2160, 1)
# ID_2520<-seq(1.2520, 355.2520, 1)
# ID_2880<-seq(1.2880, 355.2880, 1)
# ID_3240<-seq(1.3240, 355.3240, 1)
# ID_3600<-seq(1.3600, 355.3600, 1)
# ID_3960<-seq(1.3960, 355.3960, 1)
# ID_4320<-seq(1.4320, 355.4320, 1)
# ID_4680<-seq(1.4680, 355.4680, 1)
# ID_5040<-seq(1.5040, 355.5040, 1)
# ID_5400<-seq(1.5400, 355.5400, 1)
# ID_5760<-seq(1.5760, 355.5760, 1)
# ID_6120<-seq(1.6120, 355.6120, 1)
# ID_6480<-seq(1.6480, 355.6480, 1)
# ID_6840<-seq(1.6840, 355.6840, 1)
# ID_7200<-seq(1.7200, 355.7200, 1)
# ID_7560<-seq(1.7560, 355.7560, 1)
# ID_7920<-seq(1.7920, 355.7920, 1)
# ID_8280<-seq(1.8280, 355.8280, 1)
# ID_8640<-seq(1.8640, 355.8640, 1)
# ID_9000<-seq(1.9000, 355.9000, 1)
# ID_9360<-seq(1.9360, 355.9360, 1)
# ID_9720<-seq(1.9720, 355.9720, 1)
# ID_10080<-seq(1.10080, 355.10080, 1)
# ID_10440<-seq(1.10440, 355.10440, 1)
# ID_10800<-seq(1.10800, 355.10800, 1)
# ID_11160<-seq(1.11160, 355.11160, 1)
# ID_11520<-seq(1.11520, 355.11520, 1)
# ID_11880<-seq(1.11880, 355.11880, 1)
# ID_12240<-seq(1.12240, 355.12240, 1)
# ID_12600<-seq(1.12600, 355.12600, 1)
# ID_12960<-seq(1.12960, 355.12960, 1)
# ID_13320<-seq(1.13320, 355.13320, 1)
# ID_13680<-seq(1.13680, 355.13680, 1)
# ID_14040<-seq(1.14040, 355.14040, 1)
# ID_14400<-seq(1.14400, 355.14400, 1)
# ID_14760<-seq(1.14760, 355.14760, 1)
# ID_15120<-seq(1.15120, 355.15120, 1)
# ID_15480<-seq(1.15480, 355.15480, 1)
# ID_15840<-seq(1.15840, 355.15840, 1)
# ID_16200<-seq(1.16200, 355.16200, 1)
# ID_16560<-seq(1.16560, 355.16560, 1)
# ID_16920<-seq(1.16920, 355.16920, 1)
# ID_17280<-seq(1.17280, 355.17280, 1)
# ID_17640<-seq(1.17640, 355.17640, 1)
# ID_18000<-seq(1.18000, 355.18000, 1)
# ID_18360<-seq(1.18360, 355.18360, 1)
# ID_18720<-seq(1.18720, 355.18720, 1)
# ID_19080<-seq(1.19080, 355.19080, 1)
# ID_19440<-seq(1.19440, 355.19440, 1)
# ID_19800<-seq(1.19800, 355.19800, 1)
# ID_20160<-seq(1.20160, 355.20160, 1)

# release_ID<-c(ID_init, ID_360,ID_720,ID_1080,ID_1440,ID_1800,ID_2160,ID_2520,ID_2880,ID_3240,ID_3600,ID_3960,ID_4320,ID_4680,ID_5040,ID_5400,ID_5760,ID_6120,ID_6480,ID_6840,ID_7200,ID_7560,ID_7920,ID_8280,ID_8640,ID_9000,ID_9360,ID_9720,ID_10080,ID_10440,ID_10800,ID_11160,ID_11520,ID_11880,ID_12240,ID_12600,ID_12960,ID_13320,ID_13680,ID_14040,ID_14400,ID_14760,ID_15120,ID_15480,ID_15840,ID_16200,ID_16560,ID_16920,ID_17280,ID_17640,ID_18000,ID_18360,ID_18720,ID_19080,ID_19440,ID_19800,ID_20160)
# length(release_ID)

release_ID<-c(ID_init, ID_360, ID_720, ID_1080, ID_1440)

# create df of just the particles with above ID
data_use2<-data_use[data_use$ID %in% release_ID,]
dim(data_use2) # 1663841      15
length(unique(data_use2$ID)) # 18080 particles 

head(unique(data_use2$ID))


traj_time<-7 # unit in days. 

### create subset for all particles that have/can travel 7 days (and only have subset contain 7 day trajectory)
init_par<-subset(data_use2, ID %in% ID_init) # only particles from 1st timestep
dim(init_par) # 344544     20
rnum<-which(init_par$sim_length==15*96*traj_time)[length(which(init_par$sim_length==15*96*traj_time))[1]]

init_par7<-init_par[1:rnum,]
dim(init_par7)
head(init_par7)
### subset out end time for particles at the end of 7 days
end_par_init_par7<-init_par[rnum,3]
par_init_end<-subset(init_par7, sim_length==end_par_init_par7)
head(par_init_end)

par_360<-subset(data_use2, ID %in% ID_360)
dim(par_360) # 950015     20
rnum<-which(par_360$sim_length==release_par2[1]+(15*96*traj_time))[length(which(par_360$sim_length==release_par2[1]+(15*96*traj_time)))[1]]

par_360_7<-par_360[1:rnum,]
dim(par_360_7) 
### subset out end time for particles at the end of 7 days
end_par360<-par_360[rnum,3]
par_360_end<-subset(par_360, sim_length==end_par360)
head(par_360_end)

par_720<-subset(data_use2, ID %in% IDs_720)
dim(par_720)
rnum<-which(par_720$sim_length==release_par2[2]+(15*96*traj_time))[length(which(par_720$sim_length==release_par2[2]+(15*96*traj_time)))[1]]

par_720_7<-par_720[1:rnum,]
dim(par_720_7)
### subset out end time for particles at the end of 7 days
end_par720<-par_720[rnum,3]
par_720_end<-subset(par_720, sim_length==end_par720)
head(par_720_end)

par_1080<-subset(data_use2, ID %in% ID_1080)
dim(par_1080)
rnum<-which(par_1080$sim_length==release_par2[3]+(15*96*traj_time))[length(which(par_1080$sim_length==release_par2[3]+(15*96*traj_time)))[1]]

par_1080_7<-par_1080[1:rnum,]
dim(par_1080_7)
### subset out end time for particles at the end of 7 days
end_par1080<-par_1080[rnum,3]
par_1080_end<-subset(par_1080, sim_length==end_par1080)
head(par_1080_end)

par_1440<-subset(data_use2, ID %in% ID_1440)
dim(par_1440)
rnum<-which(par_1440$sim_length==release_par2[4]+(15*96*traj_time))[length(which(par_1440$sim_length==release_par2[4]+(15*96*traj_time)))[1]]

par_1440_7<-par_1440[1:rnum,]
dim(par_1440_7)
### subset out end time for particles at the end of 7 days
end_par1440<-par_1440[rnum,3]
par_1440_end<-subset(par_1440, sim_length==end_par1440)
head(par_1440_end)

# par_1800<-subset(data_use2, ID %in% ID_1800)
# dim(par_1800)
# rnum<-which(par_1800$sim_length==release_par[5]+(15*96*traj_time))[length(which(par_1800$sim_length==release_par[5]+(15*96*traj_time)))[1]]

# par_1800_7<-par_1800[1:rnum,]
# dim(par_1800_7)
# ### subset out end time for particles at the end of 7 days
# end_par1800<-par_1800[rnum,3]
# par_1800_end<-subset(par_1800, sim_length==end_par1800)
# head(par_1800_end)

# par_2160<-subset(data_use2, ID %in% ID_2160)
# dim(par_2160)
# rnum<-which(par_2160$sim_length==release_par[6]+(15*96*traj_time))[length(which(par_2160$sim_length==release_par[6]+(15*96*traj_time)))[1]]

# par_2160_7<-par_2160[1:rnum,]
# dim(par_2160_7)
# ### subset out end time for particles at the end of 7 days
# end_par2160<-par_2160[rnum,3]
# par_2160_end<-subset(par_2160, sim_length==end_par2160)
# head(par_2160_end)

# par_2520<-subset(data_use2, ID %in% ID_2520)
# dim(par_2520)
# rnum<-which(par_2520$sim_length==release_par[7]+(15*96*traj_time))[length(which(par_2520$sim_length==release_par[7]+(15*96*traj_time)))[1]]

# par_2520_7<-par_2520[1:rnum,]
# dim(par_2520_7)
# ### subset out end time for particles at the end of 7 days
# end_par2520<-par_2520[rnum,3]
# par_2520_end<-subset(par_2520, sim_length==end_par2520)
# head(par_2520_end)

# par_2880<-subset(data_use2, ID %in% ID_2880)
# dim(par_2880)
# rnum<-which(par_2880$sim_length==release_par[8]+(15*96*traj_time))[length(which(par_2880$sim_length==release_par[8]+(15*96*traj_time)))[1]]

# par_2880_7<-par_2880[1:rnum,]
# dim(par_2880_7)
# ### subset out end time for particles at the end of 7 days
# end_par2880<-par_2880[rnum,3]
# par_2880_end<-subset(par_2880, sim_length==end_par2880)
# head(par_2880_end)

# par_3240<-subset(data_use2, ID %in% ID_3240)
# dim(par_3240)
# rnum<-which(par_3240$sim_length==release_par[9]+(15*96*traj_time))[length(which(par_3240$sim_length==release_par[9]+(15*96*traj_time)))[1]]

# par_3240_7<-par_3240[1:rnum,]
# dim(par_3240_7)
# ### subset out end time for particles at the end of 7 days
# end_par3240<-par_3240[rnum,3]
# par_3240_end<-subset(par_3240, sim_length==end_par3240)
# head(par_3240_end)

# par_3600<-subset(data_use2, ID %in% ID_3600)
# dim(par_3600)
# rnum<-which(par_3600$sim_length==release_par[10]+(15*96*traj_time))[length(which(par_3600$sim_length==release_par[10]+(15*96*traj_time)))[1]]

# par_3600_7<-par_3600[1:rnum,]
# dim(par_3600_7)
# ### subset out end time for particles at the end of 7 days
# end_par3600<-par_3600[rnum,3]
# par_3600_end<-subset(par_3600, sim_length==end_par3600)
# head(par_3600_end)

# par_3960<-subset(data_use2, ID %in% ID_3960)
# dim(par_3960)
# rnum<-which(par_3960$sim_length==release_par[11]+(15*96*traj_time))[length(which(par_3960$sim_length==release_par[11]+(15*96*traj_time)))[1]]

# par_3960_7<-par_3960[1:rnum,]
# dim(par_3960_7)
# ### subset out end time for particles at the end of 7 days
# end_par3960<-par_3960[rnum,3]
# par_3960_end<-subset(par_3960, sim_length==end_par3960)
# head(par_3960_end)

# par_4320<-subset(data_use2, ID %in% ID_4320)
# dim(par_4320)
# rnum<-which(par_4320$sim_length==release_par[12]+(15*96*traj_time))[length(which(par_4320$sim_length==release_par[12]+(15*96*traj_time)))[1]]

# par_4320_7<-par_4320[1:rnum,]
# dim(par_4320_7)
# ### subset out end time for particles at the end of 7 days
# end_par4320<-par_4320[rnum,3]
# par_4320_end<-subset(par_4320, sim_length==end_par4320)
# head(par_4320_end)

# par_4680<-subset(data_use2, ID %in% ID_4680)
# dim(par_4680)
# rnum<-which(par_4680$sim_length==release_par[13]+(15*96*traj_time))[length(which(par_4680$sim_length==release_par[13]+(15*96*traj_time)))[1]]

# par_4680_7<-par_4680[1:rnum,]
# dim(par_4680_7)
# ### subset out end time for particles at the end of 7 days
# end_par4680<-par_4680[rnum,3]
# par_4680_end<-subset(par_4680, sim_length==end_par4680)
# head(par_4680_end)

# par_5040<-subset(data_use2, ID %in% ID_5040)
# dim(par_5040)
# rnum<-which(par_5040$sim_length==release_par[14]+(15*96*traj_time))[length(which(par_4680$sim_length==release_par[14]+(15*96*traj_time)))[1]]

# par_5040_7<-par_5040[1:rnum,]
# dim(par_5040_7)
# ### subset out end time for particles at the end of 7 days
# end_par5040<-par_5040[rnum,3]
# par_5040_end<-subset(par_5040, sim_length==end_par5040)
# head(par_5040_end)

# par_5400<-subset(data_use2, ID %in% ID_5400)
# dim(par_5400)
# rnum<-which(par_5400$sim_length==release_par[15]+(15*96*traj_time))[length(which(par_5400$sim_length==release_par[15]+(15*96*traj_time)))[1]]

# par_5400_7<-par_5400[1:rnum,]
# dim(par_5400_7)
# ### subset out end time for particles at the end of 7 days
# end_par5400<-par_5400[rnum,3]
# par_5400_end<-subset(par_5400, sim_length==end_par5400)
# head(par_5400_end)

# par_5760<-subset(data_use2, ID %in% ID_5760)
# dim(par_5760)
# rnum<-which(par_5760$sim_length==release_par[16]+(15*96*traj_time))[length(which(par_5760$sim_length==release_par[16]+(15*96*traj_time)))[1]]

# par_5760_7<-par_5760[1:rnum,]
# dim(par_5760_7)
# ### subset out end time for particles at the end of 7 days
# end_par5760<-par_5760[rnum,3]
# par_5760_end<-subset(par_5760, sim_length==end_par5760)
# head(par_5760_end)

# par_6120<-subset(data_use2, ID %in% ID_6120)
# dim(par_6120)
# rnum<-which(par_6120$sim_length==release_par[17]+(15*96*traj_time))[length(which(par_6120$sim_length==release_par[17]+(15*96*traj_time)))[1]]

# par_6120_7<-par_6120[1:rnum,]
# dim(par_6120_7)
# ### subset out end time for particles at the end of 7 days
# end_par6120<-par_6120[rnum,3]
# par_6120_end<-subset(par_6120, sim_length==end_par6120)
# head(par_6120_end)

# par_6480<-subset(data_use2, ID %in% ID_6480)
# dim(par_6480)
# rnum<-which(par_6480$sim_length==release_par[18]+(15*96*traj_time))[length(which(par_6480$sim_length==release_par[18]+(15*96*traj_time)))[1]]

# par_6480_7<-par_6480[1:rnum,]
# dim(par_6480_7)
# ### subset out end time for particles at the end of 7 days
# end_par6480<-par_6480[rnum,3]
# par_6480_end<-subset(par_6480, sim_length==end_par6480)
# head(par_6480_end)

# par_6840<-subset(data_use2, ID %in% ID_6840)
# dim(par_6840)
# rnum<-which(par_6840$sim_length==release_par[19]+(15*96*traj_time))[length(which(par_6840$sim_length==release_par[19]+(15*96*traj_time)))[1]]

# par_6840_7<-par_6840[1:rnum,]
# dim(par_6840_7)
# ### subset out end time for particles at the end of 7 days
# end_par6840<-par_6840[rnum,3]
# par_6840_end<-subset(par_6840, sim_length==end_par6840)
# head(par_6840_end)

# par_7200<-subset(data_use2, ID %in% ID_7200)
# dim(par_7200)
# rnum<-which(par_7200$sim_length==release_par[20]+(15*96*traj_time))[length(which(par_7200$sim_length==release_par[20]+(15*96*traj_time)))[1]]

# par_7200_7<-par_7200[1:rnum,]
# dim(par_7200_7)
# ### subset out end time for particles at the end of 7 days
# end_par7200<-par_7200[rnum,3]
# par_7200_end<-subset(par_7200, sim_length==end_par7200)
# head(par_7200_end)

# par_7560<-subset(data_use2, ID %in% ID_7560)
# dim(par_7560)
# rnum<-which(par_7560$sim_length==release_par[21]+(15*96*traj_time))[length(which(par_7560$sim_length==release_par[21]+(15*96*traj_time)))[1]]

# par_7560_7<-par_7560[1:rnum,]
# dim(par_7560_7)
# ### subset out end time for particles at the end of 7 days
# end_par7560<-par_7560[rnum,3]
# par_7560_end<-subset(par_7560, sim_length==end_par7560)
# head(par_7560_end)

# par_7920<-subset(data_use2, ID %in% ID_7920)
# dim(par_7920)
# rnum<-which(par_7920$sim_length==release_par[22]+(15*96*traj_time))[length(which(par_7920$sim_length==release_par[22]+(15*96*traj_time)))[1]]

# par_7920_7<-par_7920[1:rnum,]
# dim(par_7920_7)
# ### subset out end time for particles at the end of 7 days
# end_par7920<-par_7920[rnum,3]
# par_7920_end<-subset(par_7920, sim_length==end_par7920)
# head(par_7920_end)

# par_8280<-subset(data_use2, ID %in% ID_8280)
# dim(par_8280)
# rnum<-which(par_8280$sim_length==release_par[23]+(15*96*traj_time))[length(which(par_8280$sim_length==release_par[23]+(15*96*traj_time)))[1]]

# par_8280_7<-par_8280[1:rnum,]
# dim(par_8280_7)
# ### subset out end time for particles at the end of 7 days
# end_par8280<-par_8280[rnum,3]
# par_8280_end<-subset(par_8280, sim_length==end_par8280)
# head(par_8280_end)

# par_8640<-subset(data_use2, ID %in% ID_8640)
# dim(par_8640)
# rnum<-which(par_8640$sim_length==release_par[24]+(15*96*traj_time))[length(which(par_8640$sim_length==release_par[24]+(15*96*traj_time)))[1]]

# par_8640_7<-par_8640[1:rnum,]
# dim(par_8640_7)
# ### subset out end time for particles at the end of 7 days
# end_par8640<-par_8640[rnum,3]
# par_8640_end<-subset(par_8640, sim_length==end_par8640)
# head(par_8640_end)

# par_9000<-subset(data_use2, ID %in% ID_9000)
# dim(par_9000)
# rnum<-which(par_9000$sim_length==release_par[25]+(15*96*traj_time))[length(which(par_9000$sim_length==release_par[25]+(15*96*traj_time)))[1]]

# par_9000_7<-par_9000[1:rnum,]
# dim(par_9000_7)
# ### subset out end time for particles at the end of 7 days
# end_par9000<-par_9000[rnum,3]
# par_9000_end<-subset(par_9000, sim_length==end_par9000)
# head(par_9000_end)

# par_9360<-subset(data_use2, ID %in% ID_9360)
# dim(par_9360)
# rnum<-which(par_9360$sim_length==release_par[26]+(15*96*traj_time))[length(which(par_9360$sim_length==release_par[26]+(15*96*traj_time)))[1]]

# par_9360_7<-par_9360[1:rnum,]
# dim(par_9360_7)
# ### subset out end time for particles at the end of 7 days
# end_par9360<-par_9360[rnum,3]
# par_9360_end<-subset(par_9360, sim_length==end_par9360)
# head(par_9360_end)

# par_9720<-subset(data_use2, ID %in% ID_9720)
# dim(par_9720)
# rnum<-which(par_9720$sim_length==release_par[27]+(15*96*traj_time))[length(which(par_9720$sim_length==release_par[27]+(15*96*traj_time)))[1]]

# par_9720_7<-par_9720[1:rnum,]
# dim(par_9720_7)
# ### subset out end time for particles at the end of 7 days
# end_par9720<-par_9720[rnum,3]
# par_9720_end<-subset(par_9720, sim_length==end_par9720)
# head(par_9720_end)

# par_10080<-subset(data_use2, ID %in% ID_10080)
# dim(par_10080)
# rnum<-which(par_10080$sim_length==release_par[28]+(15*96*traj_time))[length(which(par_10080$sim_length==release_par[28]+(15*96*traj_time)))[1]]

# par_10080_7<-par_10080[1:rnum,]
# dim(par_10080_7)
# ### subset out end time for particles at the end of 7 days
# end_par10080<-par_10080[rnum,3]
# par_10080_end<-subset(par_10080, sim_length==end_par10080)
# head(par_10080_end)

# par_10440<-subset(data_use2, ID %in% ID_10440)
# dim(par_10440)
# rnum<-which(par_10440$sim_length==release_par[29]+(15*96*traj_time))[length(which(par_10440$sim_length==release_par[29]+(15*96*traj_time)))[1]]

# par_10440_7<-par_10440[1:rnum,]
# dim(par_10440_7)
# ### subset out end time for particles at the end of 7 days
# end_par10440<-par_10440[rnum,3]
# par_10440_end<-subset(par_10440, sim_length==end_par10440)
# head(par_10440_end)

# par_10800<-subset(data_use2, ID %in% ID_10800)
# dim(par_10800)
# rnum<-which(par_10800$sim_length==release_par[30]+(15*96*traj_time))[length(which(par_10800$sim_length==release_par[30]+(15*96*traj_time)))[1]]

# par_10800_7<-par_10800[1:rnum,]
# dim(par_10800_7)
# ### subset out end time for particles at the end of 7 days
# end_par10800<-par_10800[rnum,3]
# par_10800_end<-subset(par_10800, sim_length==end_par10800)
# head(par_10800_end)

# par_11160<-subset(data_use2, ID %in% ID_11160)
# dim(par_11160)
# rnum<-which(par_11160$sim_length==release_par[31]+(15*96*traj_time))[length(which(par_11160$sim_length==release_par[31]+(15*96*traj_time)))[1]]

# par_11160_7<-par_11160[1:rnum,]
# dim(par_11160_7)
# ### subset out end time for particles at the end of 7 days
# end_par11160<-par_11160[rnum,3]
# par_11160_end<-subset(par_11160, sim_length==end_par11160)
# head(par_11160_end)

# par_11520<-subset(data_use2, ID %in% ID_11520)
# dim(par_11520)
# rnum<-which(par_11520$sim_length==release_par[32]+(15*96*traj_time))[length(which(par_11520$sim_length==release_par[32]+(15*96*traj_time)))[1]]

# par_11520_7<-par_11520[1:rnum,]
# dim(par_11520_7)
# ### subset out end time for particles at the end of 7 days
# end_par11520<-par_11520[rnum,3]
# par_11520_end<-subset(par_11520, sim_length==end_par11520)
# head(par_11520_end)

# par_11880<-subset(data_use2, ID %in% ID_11880)
# dim(par_11880)
# rnum<-which(par_11880$sim_length==release_par[33]+(15*96*traj_time))[length(which(par_11880$sim_length==release_par[33]+(15*96*traj_time)))[1]]

# par_11880_7<-par_11880[1:rnum,]
# dim(par_11880_7)
# ### subset out end time for particles at the end of 7 days
# end_par11880<-par_11880[rnum,3]
# par_11880_end<-subset(par_11880, sim_length==end_par11880)
# head(par_11880_end)

# par_12240<-subset(data_use2, ID %in% ID_12240)
# dim(par_12240)
# rnum<-which(par_12240$sim_length==release_par[34]+(15*96*traj_time))[length(which(par_12240$sim_length==release_par[34]+(15*96*traj_time)))[1]]

# par_12240_7<-par_12240[1:rnum,]
# dim(par_12240_7)
# ### subset out end time for particles at the end of 7 days
# end_par12240<-par_12240[rnum,3]
# par_12240_end<-subset(par_12240, sim_length==end_par12240)
# head(par_12240_end)

# par_12600<-subset(data_use2, ID %in% ID_12600)
# dim(par_12600)
# rnum<-which(par_12600$sim_length==release_par[35]+(15*96*traj_time))[length(which(par_12600$sim_length==release_par[35]+(15*96*traj_time)))[1]]

# par_12600_7<-par_12600[1:rnum,]
# tail(par_12600_7)
# ### subset out end time for particles at the end of 7 days
# end_par12600<-par_12600[rnum,3]
# par_12600_end<-subset(par_12600, sim_length==end_par20160)
# dim(par_12600_end)

# par_12960<-subset(data_use2, ID %in% ID_12960)
# dim(par_12960)
# rnum<-which(par_12960$sim_length==release_par[36]+(15*96*traj_time))[length(which(par_12960$sim_length==release_par[36]+(15*96*traj_time)))[1]]

# par_12960_7<-par_12960[1:rnum,]
# dim(par_12960_7)
# ### subset out end time for particles at the end of 7 days
# end_par12960<-par_12960[rnum,3]
# par_12960_end<-subset(par_12960, sim_length==end_par12960)
# dim(par_12960_end)

# par_13320<-subset(data_use2, ID %in% ID_13320)
# dim(par_13320)
# rnum<-which(par_13320$sim_length==release_par[37]+(15*96*traj_time))[length(which(par_13320$sim_length==release_par[37]+(15*96*traj_time)))[1]]

# par_13320_7<-par_13320[1:rnum,]
# dim(par_13320_7)
# ### subset out end time for particles at the end of 7 days
# end_par13320<-par_13320[rnum,3]
# par_13320_end<-subset(par_13320, sim_length==end_par13320)
# head(par_13320_end)

# par_13680<-subset(data_use2, ID %in% ID_13680)
# dim(par_13680)
# rnum<-which(par_13680$sim_length==release_par[38]+(15*96*traj_time))[length(which(par_13680$sim_length==release_par[38]+(15*96*traj_time)))[1]]

# par_13680_7<-par_13680[1:rnum,]
# dim(par_13680_7)
# ### subset out end time for particles at the end of 7 days
# end_par13680<-par_13680[rnum,3]
# par_13680_end<-subset(par_13680, sim_length==end_par13680)
# head(par_13680_end)

# par_14040<-subset(data_use2, ID %in% ID_14040)
# dim(par_14040)
# rnum<-which(par_14040$sim_length==release_par[39]+(15*96*traj_time))[length(which(par_14040$sim_length==release_par[39]+(15*96*traj_time)))[1]]

# par_14040_7<-par_14040[1:rnum,]
# dim(par_14040_7)
# ### subset out end time for particles at the end of 7 days
# end_par14040<-par_14040[rnum,3]
# par_14040_end<-subset(par_14040, sim_length==end_par14040)
# head(par_14040_end)

# par_14400<-subset(data_use2, ID %in% ID_14400)
# dim(par_14400)
# rnum<-which(par_14400$sim_length==release_par[40]+(15*96*traj_time))[length(which(par_14400$sim_length==release_par[40]+(15*96*traj_time)))[1]]

# par_14400_7<-par_14400[1:rnum,]
# dim(par_14400_7)
# ### subset out end time for particles at the end of 7 days
# end_par14400<-par_14400[rnum,3]
# par_14400_end<-subset(par_14400, sim_length==end_par14400)
# head(par_14400_end)

# par_14760<-subset(data_use2, ID %in% ID_14760)
# dim(par_14760)
# rnum<-which(par_14760$sim_length==release_par[41]+(15*96*traj_time))[length(which(par_14760$sim_length==release_par[41]+(15*96*traj_time)))[1]]

# par_14760_7<-par_14760[1:rnum,]
# dim(par_14760_7)
# ### subset out end time for particles at the end of 7 days
# end_par14760<-par_14760[rnum,3]
# par_14760_end<-subset(par_14760, sim_length==end_par14760)
# head(par_14760_end)

# par_15120<-subset(data_use2, ID %in% ID_15120)
# dim(par_15120)
# rnum<-which(par_15120$sim_length==release_par[42]+(15*96*traj_time))[length(which(par_15120$sim_length==release_par[42]+(15*96*traj_time)))[1]]

# par_15120_7<-par_15120[1:rnum,]
# dim(par_15120_7)
# ### subset out end time for particles at the end of 7 days
# end_par15120<-par_15120[rnum,3]
# par_15120_end<-subset(par_15120, sim_length==end_par20160)
# dim(par_15120_end)

# par_15480<-subset(data_use2, ID %in% ID_15480)
# dim(par_15480)
# rnum<-which(par_15480$sim_length==release_par[43]+(15*96*traj_time))[length(which(par_15480$sim_length==release_par[43]+(15*96*traj_time)))[1]]

# par_15480_7<-par_15480[1:rnum,]
# dim(par_15480_7)
# ### subset out end time for particles at the end of 7 days
# end_par15480<-par_15480[rnum,3]
# par_15480_end<-subset(par_15480, sim_length==end_par15480)
# head(par_15480_end)

# par_15840<-subset(data_use2, ID %in% ID_15840)
# dim(par_15840)
# rnum<-which(par_15840$sim_length==release_par[44]+(15*96*traj_time))[length(which(par_15840$sim_length==release_par[44]+(15*96*traj_time)))[1]]

# par_15840_7<-par_15840[1:rnum,]
# dim(par_15840_7)
# ### subset out end time for particles at the end of 7 days
# end_par15840<-par_15840[rnum,3]
# par_15840_end<-subset(par_15840, sim_length==end_par15840)
# head(par_15840_end)

# par_16200<-subset(data_use2, ID %in% ID_16200)
# dim(par_16200)
# rnum<-which(par_16200$sim_length==release_par[45]+(15*96*traj_time))[length(which(par_16200$sim_length==release_par[45]+(15*96*traj_time)))[1]]

# par_16200_7<-par_16200[1:rnum,]
# dim(par_16200_7)
# ### subset out end time for particles at the end of 7 days
# end_par16200<-par_16200[rnum,3]
# par_16200_end<-subset(par_16200, sim_length==end_par16200)
# dim(par_20160_end)

# par_16560<-subset(data_use2, ID %in% ID_16560)
# dim(par_16560)
# rnum<-which(par_16560$sim_length==release_par[46]+(15*96*traj_time))[length(which(par_16560$sim_length==release_par[46]+(15*96*traj_time)))[1]]

# par_16560_7<-par_16560[1:rnum,]
# dim(par_16560_7)
# ### subset out end time for particles at the end of 7 days
# end_par16560<-par_16560[rnum,3]
# par_16560_end<-subset(par_16560, sim_length==end_par16560)
# dim(par_16560_end)

# par_16920<-subset(data_use2, ID %in% ID_16920)
# dim(par_16920)
# rnum<-which(par_16920$sim_length==release_par[47]+(15*96*traj_time))[length(which(par_16920$sim_length==release_par[47]+(15*96*traj_time)))[1]]

# par_16920_7<-par_16920[1:rnum,]
# dim(par_16920_7)
# ### subset out end time for particles at the end of 7 days
# end_par16920<-par_16920[rnum,3]
# par_16920_end<-subset(par_16920, sim_length==end_par16920)
# head(par_16920_end)

# par_17280<-subset(data_use2, ID %in% ID_17280)
# dim(par_17280)
# rnum<-which(par_17280$sim_length==release_par[48]+(15*96*traj_time))[length(which(par_17280$sim_length==release_par[48]+(15*96*traj_time)))[1]]

# par_17280_7<-par_17280[1:rnum,]
# dim(par_17280_7)
# ### subset out end time for particles at the end of 7 days
# end_par17280<-par_17280[rnum,3]
# par_17280_end<-subset(par_17280, sim_length==end_par17280)
# head(par_17280_end)

# par_17640<-subset(data_use2, ID %in% ID_17640)
# dim(par_17640)
# rnum<-which(par_17640$sim_length==release_par[49]+(15*96*traj_time))[length(which(par_17640$sim_length==release_par[49]+(15*96*traj_time)))[1]]

# par_17640_7<-par_17640[1:rnum,]
# dim(par_17640_7)
# ### subset out end time for particles at the end of 7 days
# end_par17640<-par_17640[rnum,3]
# par_17640_end<-subset(par_17640, sim_length==end_par17640)
# head(par_17640_end)

# par_18000<-subset(data_use2, ID %in% ID_18000)
# dim(par_ID_18000)
# rnum<-which(par_ID_18000$sim_length==release_par[50]+(15*96*traj_time))[length(which(par_ID_18000$sim_length==release_par[50]+(15*96*traj_time)))[1]]

# par_18000_7<-par_18000[1:rnum,]
# dim(par_18000_7)
# ### subset out end time for particles at the end of 7 days
# end_par18000<-par_18000[rnum,3]
# par_18000_end<-subset(par_18000, sim_length==end_par18000)
# head(par_18000_end)

# par_18360<-subset(data_use2, ID %in% ID_18360)
# dim(par_18360)
# rnum<-which(par_18360$sim_length==release_par[51]+(15*96*traj_time))[length(which(par_18360$sim_length==release_par[51]+(15*96*traj_time)))[1]]

# par_18360_7<-par_18360[1:rnum,]
# dim(par_18360_7)
# ### subset out end time for particles at the end of 7 days
# end_par18360<-par_18360[rnum,3]
# par_18360_end<-subset(par_18360, sim_length==end_par18360)
# head(par_18360_end)

# par_18720<-subset(data_use2, ID %in% ID_18720)
# dim(par_18720)
# rnum<-which(par_18720$sim_length==release_par[52]+(15*96*traj_time))[length(which(par_18720$sim_length==release_par[52]+(15*96*traj_time)))[1]]

# par_18720_7<-par_18720[1:rnum,]
# dim(par_18720_7)
# ### subset out end time for particles at the end of 7 days
# end_par18720<-par_18720[rnum,3]
# par_18720_end<-subset(par_18720, sim_length==end_par18720)
# head(par_18720_end)

# par_19080<-subset(data_use2, ID %in% ID_19080)
# dim(par_19080)
# rnum<-which(par_19080$sim_length==release_par[53]+(15*96*traj_time))[length(which(par_19080$sim_length==release_par[53]+(15*96*traj_time)))[1]]

# par_19080_7<-par_19080[1:rnum,]
# dim(par_19080_7)
# ### subset out end time for particles at the end of 7 days
# end_par19080<-par_19080[rnum,3]
# par_19080_end<-subset(par_19080, sim_length==end_par19080)
# head(par_19080_end)

# par_19440<-subset(data_use2, ID %in% ID_19440)
# dim(par_19440)
# rnum<-which(par_19440$sim_length==release_par[54]+(15*96*traj_time))[length(which(par_19440$sim_length==release_par[54]+(15*96*traj_time)))[1]]

# par_19440_7<-par_19440[1:rnum,]
# dim(par_19440_7)
# ### subset out end time for particles at the end of 7 days
# end_par19440<-par_19440[rnum,3]
# par_19440_end<-subset(par_19440, sim_length==end_par19440)
# head(par_19440_end)


# par_19800<-subset(data_use2, ID %in% ID_19800)
# dim(par_19800)
# rnum<-which(par_19800$sim_length==release_par[55]+(15*96*traj_time))[length(which(par_19800$sim_length==release_par[55]+(15*96*traj_time)))[1]]

# par_19800_7<-par_19800[1:rnum,]
# dim(par_19800_7)

# ### subset out end time for particles at the end of 7 days
# end_par19800<-par_19800[rnum,3]
# par_19800_end<-subset(par_19800, sim_length==end_par19800)
# head(par_19800_end)

# par_20160<-subset(data_use2, ID %in% ID_20160)
# dim(par_20160)
# rnum20160<-which(par_20160$sim_length==release_par[56]+(15*96*traj_time))[length(which(par_20160$sim_length==release_par[56]+(15*96*traj_time)))[1]]

# par_20160_7<-par_20160[1:rnum,]
# dim(par_20160_7)

# ### subset out end time for particles at the end of 7 days
# end_par20160<-par_20160[rnum,3]
# par_20160_end<-subset(par_20160, sim_length==end_par20160)
# head(par_20160_end)

#### FREQUENCEY OF OCCURRENCE DATAFRAMES
### dataframe for particles that travelled 7 days in the model domain
# par_7days<-rbind(init_par7,par_360_7,par_720_7,par_1080_7,par_1440_7) 

# ,par_1800_7,par_2160_7,par_2520_7,par_2880_7,par_3240_7,par_3600_7,par_3960_7,par_4320_7,par_4680_7,par_5040_7,par_5400_7,par_5760_7,par_6120_7,par_6480_7,par_6840_7,par_7200_7,par_7560_7,par_7920_7,par_8280_7,par_8640_7,par_9000_7,par_9360_7,par_9720_7,par_10080_7,par_10440_7,par_10800_7,par_11160_7,par_11520_7,par_11880_7,par_12240_7,par_12600_7,par_12960_7,par_13320_7,par_13680_7,par_14040_7,par_14400_7,par_14760_7,par_15120_7,par_15480_7,par_15840_7,par_16200_7,par_16560_7, par_16920_7, par_17280_7, par_17640_7, par_18000_7, par_18360_7, par_18720_7, par_19080_7, par_19440_7, par_19800_7, par_20160_7

dim(par_7days) # 1661224      15



# par_14days<- rbind(init_par7,par_360_7,par_720_7,par_1080_7,par_1440_7) 
#rbind(init_par7,par_360_7,par_720_7,par_1080_7,par_1440_7,par_1800_7,par_2160_7,par_2520_7,par_2880_7,par_3240_7,par_3600_7,par_3960_7,par_4320_7,par_4680_7,par_5040_7,par_5400_7,par_5760_7,par_6120_7,par_6480_7,par_6840_7,par_7200_7,par_7560_7,par_7920_7,par_8280_7,par_8640_7,par_9000_7,par_9360_7,par_9720_7,par_10080_7,par_10440_7,par_10800_7,par_11160_7,par_11520_7,par_11880_7,par_12240_7,par_12600_7,par_12960_7,par_13320_7,par_13680_7,par_14040_7,par_14400_7,par_14760_7,par_15120_7,par_15480_7,par_15840_7,par_16200_7,par_16560_7, par_16920_7, par_17280_7, par_17640_7, par_18000_7, par_18360_7, par_18720_7, par_19080_7, par_19440_7, par_19800_7, par_20160_7)

dim(par_14days) # 1588200      15


par_21days<-rbind(init_par7,par_360_7,par_720_7,par_1080_7,par_1440_7) # rbind(init_par7,par_360_7,par_720_7,par_1080_7,par_1440_7,par_1800_7,par_2160_7,par_2520_7,par_2880_7,par_3240_7,par_3600_7,par_3960_7,par_4320_7,par_4680_7,par_5040_7,par_5400_7,par_5760_7,par_6120_7,par_6480_7,par_6840_7,par_7200_7,par_7560_7,par_7920_7,par_8280_7,par_8640_7,par_9000_7,par_9360_7,par_9720_7,par_10080_7,par_10440_7,par_10800_7,par_11160_7,par_11520_7,par_11880_7,par_12240_7,par_12600_7,par_12960_7,par_13320_7,par_13680_7,par_14040_7,par_14400_7,par_14760_7,par_15120_7,par_15480_7,par_15840_7,par_16200_7,par_16560_7, par_16920_7, par_17280_7, par_17640_7, par_18000_7, par_18360_7, par_18720_7, par_19080_7, par_19440_7, par_19800_7, par_20160_7)

dim(par_21days) # 1591595      15


par_36days<-init_par
dim(par_36days)


#### IMPORT DATA FRAMES FOR DIFFERENT TIME PERIODS

par_7end<-rbind(par_init_end, par_360_end,par_720_end,par_1080_end,par_1440_end) # ,par_1800_end,par_2160_end,par_2520_end,par_2880_end,par_3240_end,par_3600_end,par_3960_end,par_4320_end,par_4680_end,par_5040_end,par_5400_end,par_5760_end,par_6120_end,par_6480_end,par_6840_end,par_7200_end,par_7560_end,par_7920_end,par_8280_end,par_8640_end,par_9000_end,par_9360_end,par_9720_end,par_10080_end,par_10440_end,par_10800_end,par_11160_end,par_11520_end,par_11880_end,par_12240_end,par_12600_end,par_12960_end,par_13320_end,par_13680_end,par_14040_end,par_14400_end,par_14760_end,par_15120_end,par_15480_end,par_15840_end,par_16200_end,par_16560_end,par_16920_end,par_17280_end,par_17640_end,par_18000_end,par_18360_end,par_18720_end,par_19080_end,par_19440_end,par_19800_end,par_20160_end)

dim(par_7end) # 745  15

par_14end<-rbind(par_init_end, par_360_end,par_720_end,par_1080_end,par_1440_end) # rbind(par_init_end, par_360_end,par_720_end,par_1080_end,par_1440_end,par_1800_end,par_2160_end,par_2520_end,par_2880_end,par_3240_end,par_3600_end,par_3960_end,par_4320_end,par_4680_end,par_5040_end,par_5400_end,par_5760_end,par_6120_end,par_6480_end,par_6840_end,par_7200_end,par_7560_end,par_7920_end,par_8280_end,par_8640_end,par_9000_end,par_9360_end,par_9720_end,par_10080_end,par_10440_end,par_10800_end,par_11160_end,par_11520_end,par_11880_end,par_12240_end,par_12600_end,par_12960_end,par_13320_end,par_13680_end,par_14040_end,par_14400_end,par_14760_end,par_15120_end,par_15480_end,par_15840_end,par_16200_end,par_16560_end,par_16920_end,par_17280_end,par_17640_end,par_18000_end,par_18360_end,par_18720_end,par_19080_end,par_19440_end,par_19800_end,par_20160_end)

dim(par_14end) # 442  15

par_21end<-rbind(par_init_end, par_360_end,par_720_end,par_1080_end,par_1440_end) # rbind(par_init_end, par_360_end,par_720_end,par_1080_end,par_1440_end,par_1800_end,par_2160_end,par_2520_end,par_2880_end,par_3240_end,par_3600_end,par_3960_end,par_4320_end,par_4680_end,par_5040_end,par_5400_end,par_5760_end,par_6120_end,par_6480_end,par_6840_end,par_7200_end,par_7560_end,par_7920_end,par_8280_end,par_8640_end,par_9000_end,par_9360_end,par_9720_end,par_10080_end,par_10440_end,par_10800_end,par_11160_end,par_11520_end,par_11880_end,par_12240_end,par_12600_end,par_12960_end,par_13320_end,par_13680_end,par_14040_end,par_14400_end,par_14760_end,par_15120_end,par_15480_end,par_15840_end,par_16200_end,par_16560_end,par_16920_end,par_17280_end,par_17640_end,par_18000_end,par_18360_end,par_18720_end,par_19080_end,par_19440_end,par_19800_end,par_20160_end)

dim(par_21end) # 365   15

end_sim_init<-max(init_par$sim_length, na.rm=T)
par_36end<-subset(init_par, sim_length==end_sim_init)
dim(par_36end) # 30  15




####################################################################
### Making Frequency of Occurrence Heat maps

#Part 1: make a regular grid and count stations within each grid cell. Plot results
   
    nlat= grid.size.lat/2# 698/10 # determine resolution of grid
    nlon= grid.size.lon/2 # 238/10
    latd=seq(40.65895, 46.99730,length.out=nlat)
    lond=seq(-126.0909,-123.9263,length.out=nlon)
    
    
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
  # dev.new(width=4,height=10)
  plot(init_par$lon, init_par$lat,pch='.',ylim=c(40.65895, 46.99730),xlim=c(-126.0909,-123.9263))
  # points(maynew2$lon, maynew2$lat, col='red', pch= '.')
  n.stations=NA*(1:nrow(grid.lon))
   
  for(i in 1:length(n.stations)){
    print(i)
    tmp=in.chull(init_par$lon, init_par$lat, grid.lon[i,], grid.lat[i,])
    n.stations[i]=sum(tmp)#This decides what goes into each grid pixel (may16_3$ID*tmp)
    points(init_par$lon[tmp], init_par$lat[tmp],col=i,pch=16)
    polygon(grid.lon[i,],grid.lat[i,])
  }
  map("worldHires",fill=T,col="grey",add=T)
  
num_ts<-nrow(init_par)


## creating grid to plot for heatmap 
z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
z.matrix_36days_FoO<-((matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F))/num_ts)*100 # /uniID_nov16 # units are percentage (%)
range(z.matrix_36days_FoO, na.rm=T) # 21 117

write.csv(z.matrix_36days_FoO, "FreqofOccur_36dayrelease_lowres_20depth_HiResCompare.csv")


### plot data 

pal<-colorRampPalette(c("white","cyan", "cyan4", "darkblue"))

setwd("/Users/wongalaj/Desktop")
png('freqofoccurrence_7days_depth40_lowres', width=800, height=2000, res=150)
par(mai=c(1.1,1.1,0.3,0.3))

image(z.lon,z.lat, z.matrix_7days_FoO,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main="7 day release", col=pal(500), zlim= c(0,2), cex.lab=2,cex.axis=2)

map("worldHires",fill=T,col="grey", add=T)

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,3000,by=300),labcex=0.7,add=T,col='black', lwd= 0.1) 

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
  
# points(-124.52,42.84, pch=20, cex= 1.2)
text(-124.49,42.95,"Cape Blanco", adj=c(0,1.2))

# add in zlim legend for map (only for 7 day release figure) 
# par(new=T,mai=c(0,0.1,5,24))
# image.plot(legend.only=T,col=pal(500),zlim= c(0,2),legend.width=2,legend.shrink= 0.5)

# range(z.matrix_7days_FoO, na.rm=T)

####################################################################
### Import plots

## create grid of only initial locations for lat and lon
    nlat.end= grid.size.lat/2 # 698/10 # determine resolution of grid
    nlon.end= grid.size.lon/2 # 238/10
    latd.end=seq(40.65895 ,47.01520,length.out=nlat.end)
    lond.end=seq(-124.5884 ,-123.9263,length.out=nlon.end) # -124.2
    
  gridlon.end=data.frame(
    lon1=rep(lond.end[-length(lond.end)],(nlat.end-1)),
    lon2=rep(lond.end[-1],(nlat.end-1)),
    lon3=rep(lond.end[-1],(nlat.end-1)),
    lon4=rep(lond.end[-length(lond.end)],(nlat.end-1))) # make dataframe of just longitude
  
  gridlat.end=data.frame(
    lat1=sort(rep(latd.end[-length(latd.end)],(nlon.end-1))),
    lat2=sort(rep(latd.end[-length(latd.end)],(nlon.end-1))),
    lat3=sort(rep(latd.end[-1],(nlon.end-1))),
    lat4=sort(rep(latd.end[-1],(nlon.end-1)))) # make dataframe of just latitude

## subset data so that it only includes last timestep for each particle because I want their end locations



# plot what data looks like for each end location subset
plot(par_7end$lon, par_7end$lat, col='red', pch= '.')
n.stations.end=NA*(1:nrow(gridlon.end))

# figure out what grid cell is in each 
for(i in 1:length(n.stations.end)){
	print(i)
	tmp.end=in.chull(par_36end$lon, par_36end$lat, gridlon.end[i,], gridlat.end[i,])
	n.stations.end[i]=sum(tmp.end)#This decides what goes into each grid pixel (may16_3$ID*tmp)
	points(par_36end$lon[tmp.end], par_36end$lat[tmp.end],col=i,pch=16)
	polygon(gridlon.end[i,], gridlat.end[i,])
	}


## creating grid to plot for heatmap 
z.lat.end<-(latd.end[1:(length(latd.end)-1)]+latd.end[2:length(latd.end)])/2
z.lon.end<-(lond.end[1:(length(lond.end)-1)]+lond.end[2:length(lond.end)])/2
z.matrix.end<-(matrix(n.stations.end,ncol=length(z.lat.end),nrow=length(z.lon.end),byrow=F)) 
dim(z.matrix.end) # 21 117



# z.matrix.IMP21<-z.matrix.end # create this matrix so it can be used below
# z.matrix.IMP36<-z.matrix.end # create this matrix so it can be used below
# z.matrix.IMP14<-z.matrix.end # create this matrix so it can be used below
# z.matrix.IMP7<-z.matrix.end # create this matrix so it can be used below
# range(z.matrix.IMP7, na.rm=T)

write.csv(z.matrix.IMP36, 'Import_36dayrelease_lowres_20depth_hirescomparison.csv')

# plot heatmap of import data

pal2<-colorRampPalette(c("white","mediumorchid1", "mediumorchid3", "mediumorchid4"))

dev.new(width=4.5, height=10)

image(z.lon.end,z.lat.end, z.matrix.IMP36,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main="36 day release", col=pal2(500), zlim= c(0,10)) #, xlim=c(-126.0909,-123.9263)
  
map("worldHires",fill=T,col="grey", add=T)

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,3000,by=300),labcex=0.7,add=T,col='black', lwd= 0.2) 

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
  
text(-124.49,42.95,"Cape Blanco", adj=c(0,1.2))

# add in zlim legend for map (this only needs to be in the first figure) 
par(new=T,mai=c(0,0.1,3,5))
image.plot(legend.only=T,col=pal2(500),zlim=c(0,10),legend.width=1.3, legend.shrink= 0.3 )
  
# range(z.matrix.IMP14, na.rm=T)

 

####################################################################
### Local Retention of particles

## - Def: the percent of particles that are locally retained (aka the reigon that they were released from)
# use the below dataframes as end locations for particles
dim(par_7end) # par_14end, par21_end, par_init_end

range(pop.roms[[1]][[1]]$lon[1:355], na.rm=T) # this will give us the starting latitudes for all the particles released because they are all released at the same locations just over different time scales

### subset initial pop for comparison
pop_init<-subset(may2018, del.t==1) # pop.roms[[1]][[1]]
dim(pop_init) # 355 14 
head(pop_init)
range(pop_init$lon)

### lat and lon grid needed to figure out if particles are in a certain grid cell 
 	
    latd.local=seq(40.65895, 47.01520,length.out=grid.size.lat/2) # grid.size.lat/2
    lond.local=seq(-124.5884 ,-123.9263, length.out= grid.size.lon/2) # grid.size.lat/2
    nlat.local= length(latd.local) # 698/10 # determine resolution of grid
    nlon.local= length(lond.local) # 238/10
    
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
plot(par_36end$lon, par_36end$lat, col='red', pch= '.', ylim=c(40.65895, 46.99730),xlim=c(-126.0909,-123.9263))
points(pop_init$lon, pop_init$lat, col='blue', pch= '.')

n.stations.start=NA*(1:nrow(gridlon.local)) # empty vector to be filled 
n.stations.end=NA*(1:nrow(gridlon.local)) # empty vector to be filled 

### count particles for start locations (this always stays as pop_init data.frame)
for(i in 1:length(n.stations.start)){
	print(i)
	tmp.local=in.chull(pop_init$lon, pop_init$lat, gridlon.local[i,], gridlat.local[i,])
	n.stations.start[i]=sum(tmp.local) # This decides what goes into each grid pixel
	points(pop_init$lon[tmp.local], pop_init$lat[tmp.local],col=i,pch=16)
	polygon(gridlon.local[i,], gridlat.local[i,])
	} 

# ### count particles for end locations
for(i in 1:length(n.stations.end)){
	print(i)
	tmp.local=in.chull(par_36end$lon, par_36end$lat, gridlon.local[i,], gridlat.local[i,])
	n.stations.end[i]=sum(tmp.local) # This decides what goes into each grid pixel
	points(par_36end$lon[tmp.local], par_36end$lat[tmp.local],col=i,pch=16)
	polygon(gridlon.local[i,], gridlat.local[i,])
	} 

## creating grid to plot for start heatmap
z.lat.local<-(latd.local[1:(length(latd.local)-1)]+latd.local[2:length(latd.local)])/2 # lat grid to plot heatmap
z.lon.local<-(lond.local[1:(length(lond.local)-1)]+lond.local[2:length(lond.local)])/2 # lon grid to plot heatmap
z.matrix.start<-(matrix(n.stations.start,ncol=length(z.lat.local),nrow=length(z.lon.local),byrow=F))
dim(z.matrix.start) # 10 58
range(z.matrix.start, na.rm=T)

# ## creating grid to plot for start heatmap
z.lat.local<-(latd.local[1:(length(latd.local)-1)]+latd.local[2:length(latd.local)])/2 # lat grid to plot heatmap
z.lon.local<-(lond.local[1:(length(lond.local)-1)]+lond.local[2:length(lond.local)])/2 # lon grid to plot heatmap
z.matrix.end2<-(matrix(n.stations.end,ncol=length(z.lat.local),nrow=length(z.lon.local),byrow=F))
dim(z.matrix.end2) # 10 58
range(z.matrix.end2)

# create local retention metric to plot (difference between # of end particles/# of start particles)
z.matrix.LR<- z.matrix.end2/z.matrix.start # local recruitment heatmap
dim(z.matrix.LR) # 21 117
range(z.matrix.LR, na.rm=T) # 0 4 
z.matrix.LR[which(z.matrix.LR=='Inf')]<-1 # remove Inf to maxed out number (make sure zlim is less than 500)
hist(z.matrix.LR) # ALWAYS CHECK HISTOGRAM BEFORE MAKING INF A NUMBER

# z.matrix.LR21<-z.matrix.LR # create this matrix so it can be used below
# z.matrix.LR36<-z.matrix.LR # create this matrix so it can be used below
# z.matrix.LR14<-z.matrix.LR # create this matrix so it can be used below
# z.matrix.LR7<-z.matrix.LR # create this matrix so it can be used below

write.csv(z.matrix.LR36, 'LocalRetention_36dayrelease_lowres_20depth_comparison.csv')

# plot heat map of local recruitment metrics (plotted for different depths and release lengths)

pal3<-colorRampPalette(c('white', 'darkolivegreen1', 'darkolivegreen2','darkolivegreen3', 'darkolivegreen'))

dev.new(width=4.5, height=10)

image(z.lon.local, z.lat.local, z.matrix.LR36,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main="36 day release", col= pal3(500), zlim=c(0,7))  

map("worldHires",fill=T,col="grey", add=T)

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,3000,by=300),labcex=0.7,add=T,col='black', lwd= 0.2) 

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
  
text(-124.49,42.95,"Cape Blanco", adj=c(0,1.2))

# add in zlim legend for map (this only needs to be in the first figure) 
par(new=T,mai=c(0,0.1,3,5))
image.plot(legend.only=T,col=pal3(500),zlim=c(0,7),legend.width=1.3, legend.shrink= 0.3 )
####################################################################
### Self Recruitment Parameter

## Def: this parameter is defined as Local Retention divided by Import (SR=LR/Import)

## local recruitment datasets for different release times 
# z.matrix.LR21 # 21 day release
# z.matrix.LR36 # 36 day release
# z.matrix.LR14 # 14 day release
# z.matrix.LR7 # 7 day release

## import datasets for different release times
# z.matrix.IMP21 # 21 day release
# z.matrix.IMP36 # 36 day release
# z.matrix.IMP14 # 14 day release
# z.matrix.IMP7 # 7 day release

## self recruitment calculation for each release period
z.matrix.SR36<-z.matrix.LR36/z.matrix.IMP36
dim(z.matrix.SR36) # 42 235 
range(z.matrix.SR36, na.rm=T) # just a bunch of NaNs 
z.matrix.SR36[which(z.matrix.SR36=='Inf')]<-2
hist(z.matrix.SR36)

z.matrix.SR21<-z.matrix.LR21/z.matrix.IMP21
dim(z.matrix.SR21)
range(z.matrix.SR21, na.rm=T)
hist(z.matrix.SR21)
z.matrix.SR21[which(z.matrix.SR21=='Inf')]<-10

z.matrix.SR14<-z.matrix.LR14/z.matrix.IMP14
dim(z.matrix.SR14)
range(z.matrix.SR14, na.rm=T)
hist(z.matrix.SR14)
z.matrix.SR14[which(z.matrix.SR14=='Inf')]<-10

z.matrix.SR7<-z.matrix.LR7/z.matrix.IMP7
dim(z.matrix.SR7) # 21 117
range(z.matrix.SR7, na.rm=T)
z.matrix.SR7[which(z.matrix.SR7=='Inf')]<-10
hist(z.matrix.SR7) # CHECK HISTOGRAM FIRST TO SEE WHAT IS MAX VALUE AND THEN ASSIGN INF TO BE A NUMBER

## save data
write.csv(z.matrix.SR36, 'SelfRecruitment_36dayrelease_lowres_20depth_highrescomparison.csv')



#### Plot self recruitment metric 

pal4<-colorRampPalette(c('white', 'orange1', 'orange2', 'orange3', 'orange4'))

dev.new(width=4.5, height=10)

image(z.lon.local, z.lat.local, z.matrix.SR36,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main="36 day release", col=pal4(500), zlim=c(0,5)) 
 
map("worldHires",fill=T,col="grey", add=T)

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,3000,by=300),labcex=0.7,add=T,col='black', lwd= 0.2) 

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
  
text(-124.49,42.95,"Cape Blanco", adj=c(0,1.2))

# add in zlim legend for map (this only needs to be in the first figure) 
par(new=T,mai=c(0,0.1,3,5))
image.plot(legend.only=T,col=pal4(500),zlim=c(0,5),legend.width=1.3, legend.shrink= 0.3 )





# range(pop[[1]]$lon.index) # 1 355


plot(((data_init$del.t)*15)/(1440), data_init$lat, pch=19, cex= 0.2, col= data_init$ID, xlab= 'Time (days)', ylab=expression(paste("Latitude ("^o,'N)')))

plot(((data_360$del.t)*15)/(1440), data_360$lat, pch=19, cex= 0.2, col= data_init$ID, xlab= 'Time (days)', ylab=expression(paste("Latitude ("^o,'N)')))


##########################################
### Particle dispersal distances figure

### Plot of PARTICLE DISPERSAL TRAJECTORIES (TIME VS. DISTANCE)
# plot(pop.use$del.t,pop.use$lat,col=pop.use$ID)


plot(pop.fin$timestep_day[pop.fin$ID==1],pop.fin$lat[pop.fin$ID==1],col=pop.fin$ID[pop.fin$ID==1],ylim=c(40.65895, 46.99730),ylab= expression(paste("Latitude ("^o,'N)')), xlab= "timestep (days)", type= "l", lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.360"],pop.fin$lat[pop.fin$ID=="1.360"],col="darkred", lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.720"],pop.fin$lat[pop.fin$ID=="1.720"],col='darkorange3', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.1080"],pop.fin$lat[pop.fin$ID=="1.1080"],col='goldenrod4', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.1440"],pop.fin$lat[pop.fin$ID=="1.1440"],col='chartreuse3', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.1800"],pop.fin$lat[pop.fin$ID=="1.1800"],col='forestgreen', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.2160"],pop.fin$lat[pop.fin$ID=="1.2160"],col='darkslategray4', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.2520"],pop.fin$lat[pop.fin$ID=="1.2520"],col='cadetblue3', lwd= 3)
# lines(pop.fin$timestep_day[pop.fin$ID==9],pop.fin$lat[pop.fin$ID==9],col='darkorchid3', lwd= 3)
# lines(pop.fin$timestep_day[pop.fin$ID==10],pop.fin$lat[pop.fin$ID==10],col='brown', lwd= 3)
# lines(pop.fin$timestep_day[pop.fin$ID==11],pop.fin$lat[pop.fin$ID==11],col='gray', lwd= 3)
legend(x= 'topright', y=NULL, legend= seq(1,4,1), lty= 1, col = c('black', "darkred", 'darkorange3', 'goldenrod4', 'chartreuse3', 'forestgreen', 'darkslategray4', 'cadetblue3', 'darkorchid3', 'brown', 'gray'), title= 'Particle IDs', lwd= 3)

plot(pop.fin$timestep_day, pop.fin$lat, type='l', ylim=c(40, 47)) 


##########################################
## Bathy with trajectories and vector plots
# setwd("/Users/Jenn/Desktop/UW_Cascadia/OR_coast_AK/OR_coast_AK_grid")
# pop.fin<-read.csv("OR_coast_20190527_server_TEST1.csv")

## colors
blues<-c("lightsteelblue4",  "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))
rainbow<-rainbow(220)


## subset each particle ID
pop1<-subset(pop.fin, ID==1)
pop21<-subset(pop.fin, ID=="1.360")
pop41<-subset(pop.fin, ID=="1.720")
pop61<-subset(pop.fin, ID=="1.1080")
pop81<-subset(pop.fin, ID=="1.1440")
pop101<-subset(pop.fin, ID=="1.1800")
pop121<-subset(pop.fin, ID=="1.2160")
pop141<-subset(pop.fin, ID=="1.2520")
# pop161<-subset(pop.fin, ID==161)
# pop181<-subset(pop.fin, ID==181)
# pop201<-subset(pop.fin, ID==201)
# pop221<-subset(pop.fin, ID==221)
# pop241<-subset(pop.fin, ID==241)
# pop261<-subset(pop.fin, ID==261)
# pop281<-subset(pop.fin, ID==281)
# pop301<-subset(pop.fin, ID==301)
# pop321<-subset(pop.fin, ID==321)
# pop341<-subset(pop.fin, ID==341)
# pop361<-subset(pop.fin, ID==361)
# pop381<-subset(pop.fin, ID==381)
# pop401<-subset(pop.fin, ID==401)
# pop421<-subset(pop.fin, ID==421)
# pop441<-subset(pop.fin, ID==441)
# pop461<-subset(pop.fin, ID==461)
# pop481<-subset(pop.fin, ID==481)
# pop501<-subset(pop.fin, ID==501)
# pop521<-subset(pop.fin, ID==521)
# pop541<-subset(pop.fin, ID==541)
# pop561<-subset(pop.fin, ID==561)
# pop581<-subset(pop.fin, ID==581)
# pop601<-subset(pop.fin, ID==601)
# pop621<-subset(pop.fin, ID==621)
# pop641<-subset(pop.fin, ID==641)
# pop661<-subset(pop.fin, ID==661)
# pop681<-subset(pop.fin, ID==681)
# pop701<-subset(pop.fin, ID==701)



pop1new<-pop1[seq(1, nrow(pop1), 8),]
plt_del.t<-seq(1, dim(pop1new)[1] ,1) # every 8th iteration of 30 day run is 360 
pop1new<-cbind(plt_del.t,pop1new)
# head(pop1new)
# head(pop1new)
pop1new[is.na(pop1new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop21new<-pop21[seq(1, nrow(pop21), 8),]
dim(pop21new)
plt_del.t<-seq(1, dim(pop21new)[1] ,1)
pop21new<-cbind(plt_del.t,pop21new)
head(pop21new)
pop21new[is.na(pop21new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop41new<-pop41[seq(1, nrow(pop41), 8),]
dim(pop41new)
plt_del.t<-seq(1, dim(pop41new)[1] ,1)
pop41new<-cbind(plt_del.t,pop41new)
head(pop41new)
pop41new[is.na(pop41new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop61new<-pop61[seq(1, nrow(pop61), 8),]
dim(pop61new)
plt_del.t<-seq(1, dim(pop61new)[1] ,1)
pop61new<-cbind(plt_del.t,pop61new)
head(pop61new)
dim(pop61new)
pop61new[is.na(pop61new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop81new<-pop81[seq(1, nrow(pop81), 8),]
dim(pop81new)
plt_del.t<-seq(1, dim(pop81new)[1] ,1)
pop81new<-cbind(plt_del.t,pop81new)
head(pop81new)
dim(pop81new)
pop81new[is.na(pop81new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop101new<-pop101[seq(1, nrow(pop101), 8),]
dim(pop101new)
plt_del.t<-seq(1, dim(pop101new)[1] ,1)
pop101new<-cbind(plt_del.t,pop101new)
head(pop101new)
dim(pop101new)
pop101new[is.na(pop101new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop121new<-pop121[seq(1, nrow(pop121), 8),]
dim(pop121new)
plt_del.t<-seq(1, dim(pop121new)[1] ,1)
pop121new<-cbind(plt_del.t,pop121new)
head(pop121new)
dim(pop121new)
pop121new[is.na(pop121new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop141new<-pop141[seq(1, nrow(pop141), 8),]
dim(pop141new)
plt_del.t<-seq(1, dim(pop141new)[1] ,1)
pop141new<-cbind(plt_del.t,pop141new)
head(pop141new)
dim(pop141new)
pop141new[is.na(pop141new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop161new<-pop161[seq(1, nrow(pop161), 8),]
plt_del.t<-seq(1, dim(pop161new)[1] ,1)
pop161new<-cbind(plt_del.t,pop161new)
head(pop161new)
dim(pop161new)
pop161new[is.na(pop161new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop181new<-pop181[seq(1, nrow(pop181), 8),]
plt_del.t<-seq(1, dim(pop181new)[1] ,1)
pop181new<-cbind(plt_del.t,pop181new)
head(pop181new)
dim(pop181new)
pop181new[is.na(pop181new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop201new<-pop201[seq(1, nrow(pop201), 8),]
plt_del.t<-seq(1, dim(pop201new)[1] ,1)
pop201new<-cbind(plt_del.t,pop201new)
head(pop201new)
dim(pop201new)
pop201new[is.na(pop201new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop221new<-pop221[seq(1, nrow(pop221), 8),]
plt_del.t<-seq(1, dim(pop221new)[1] ,1)
pop221new<-cbind(plt_del.t,pop221new)
head(pop221new)
dim(pop221new)
pop221new[is.na(pop221new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop241new<-pop241[seq(1, nrow(pop241), 8),]
plt_del.t<-seq(1, dim(pop241new)[1] ,1)
pop241new<-cbind(plt_del.t,pop241new)
head(pop241new)
dim(pop241new)
pop241new[is.na(pop241new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop261new<-pop261[seq(1, nrow(pop261), 8),]
plt_del.t<-seq(1, dim(pop261new)[1] ,1)
pop261new<-cbind(plt_del.t,pop261new)
head(pop261new)
dim(pop261new)
pop261new[is.na(pop261new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop281new<-pop281[seq(1, nrow(pop281), 8),]
plt_del.t<-seq(1, dim(pop281new)[1] ,1)
pop281new<-cbind(plt_del.t,pop281new)
head(pop281new)
dim(pop281new)
pop281new[is.na(pop281new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop301new<-pop301[seq(1, nrow(pop301), 8),]
plt_del.t<-seq(1, dim(pop301new)[1] ,1)
pop301new<-cbind(plt_del.t,pop301new)
head(pop301new)
dim(pop301new)
pop301new[is.na(pop301new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop321new<-pop321[seq(1, nrow(pop321), 8),]
plt_del.t<-seq(1, dim(pop321new)[1] ,1)
pop321new<-cbind(plt_del.t,pop321new)
head(pop321new)
dim(pop321new)
pop321new[is.na(pop321new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop341new<-pop341[seq(1, nrow(pop341), 8),]
plt_del.t<-seq(1, dim(pop341new)[1] ,1)
pop341new<-cbind(plt_del.t,pop341new)
head(pop341new)
dim(pop341new)
pop341new[is.na(pop341new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop361new<-pop361[seq(1, nrow(pop361), 8),]
plt_del.t<-seq(1, dim(pop361new)[1] ,1)
pop361new<-cbind(plt_del.t,pop361new)
head(pop361new)
dim(pop361new)
pop361new[is.na(pop361new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop381new<-pop381[seq(1, nrow(pop381), 8),]
plt_del.t<-seq(1, dim(pop381new)[1] ,1)
pop381new<-cbind(plt_del.t,pop381new)
head(pop381new)
dim(pop381new)
pop381new[is.na(pop381new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop401new<-pop401[seq(1, nrow(pop401), 8),]
plt_del.t<-seq(1, dim(pop401new)[1] ,1)
pop401new<-cbind(plt_del.t,pop401new)
head(pop401new)
dim(pop401new)
pop401new[is.na(pop401new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop421new<-pop421[seq(1, nrow(pop421), 8),]
plt_del.t<-seq(1, dim(pop421new)[1] ,1)
pop421new<-cbind(plt_del.t,pop421new)
head(pop421new)
dim(pop421new)
pop421new[is.na(pop421new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop441new<-pop441[seq(1, nrow(pop441), 8),]
plt_del.t<-seq(1, dim(pop441new)[1] ,1)
pop441new<-cbind(plt_del.t,pop441new)
head(pop441new)
dim(pop441new)
pop441new[is.na(pop441new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop461new<-pop461[seq(1, nrow(pop461), 8),]
plt_del.t<-seq(1, dim(pop461new)[1] ,1)
pop461new<-cbind(plt_del.t,pop461new)
head(pop461new)
dim(pop461new)
pop461new[is.na(pop461new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop481new<-pop481[seq(1, nrow(pop481), 8),]
plt_del.t<-seq(1, dim(pop481new)[1] ,1)
pop481new<-cbind(plt_del.t,pop481new)
head(pop481new)
dim(pop481new)
pop481new[is.na(pop481new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs

pop501new<-pop501[seq(1, nrow(pop501), 8),]
plt_del.t<-seq(1, dim(pop501new)[1] ,1)
pop501new<-cbind(plt_del.t,pop501new)
head(pop501new)
dim(pop501new)
pop501new[is.na(pop501new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop521new<-pop521[seq(1, nrow(pop521), 8),]
plt_del.t<-seq(1, dim(pop521new)[1] ,1)
pop521new<-cbind(plt_del.t,pop521new)
head(pop521new)
dim(pop521new)
pop521new[is.na(pop521new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop541new<-pop541[seq(1, nrow(pop541), 8),]
plt_del.t<-seq(1, dim(pop541new)[1] ,1)
pop541new<-cbind(plt_del.t,pop541new)
head(pop541new)
dim(pop541new)
pop541new[is.na(pop541new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop561new<-pop561[seq(1, nrow(pop561), 8),]
plt_del.t<-seq(1, dim(pop561new)[1] ,1)
pop561new<-cbind(plt_del.t,pop561new)
head(pop561new)
dim(pop561new)
pop561new[is.na(pop561new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop581new<-pop581[seq(1, nrow(pop581), 8),]
plt_del.t<-seq(1, dim(pop581new)[1] ,1)
pop581new<-cbind(plt_del.t,pop581new)
head(pop581new)
dim(pop581new)
pop581new[is.na(pop581new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop601new<-pop441[seq(1, nrow(pop601), 8),]
plt_del.t<-seq(1, dim(pop601new)[1] ,1)
pop601new<-cbind(plt_del.t,pop601new)
head(pop601new)
dim(pop601new)
pop601new[is.na(pop601new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop621new<-pop621[seq(1, nrow(pop621), 8),]
plt_del.t<-seq(1, dim(pop621new)[1] ,1)
pop621new<-cbind(plt_del.t,pop621new)
head(pop621new)
dim(pop621new)
pop621new[is.na(pop621new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop641new<-pop641[seq(1, nrow(pop641), 8),]
plt_del.t<-seq(1, dim(pop641new)[1] ,1)
pop641new<-cbind(plt_del.t,pop641new)
head(pop641new)
dim(pop641new)
pop641new[is.na(pop641new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop661new<-pop661[seq(1, nrow(pop661), 8),]
plt_del.t<-seq(1, dim(pop661new)[1] ,1)
pop661new<-cbind(plt_del.t,pop661new)
head(pop661new)
dim(pop661new)
pop661new[is.na(pop661new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop681new<-pop681[seq(1, nrow(pop681), 8),]
plt_del.t<-seq(1, dim(pop681new)[1] ,1)
pop681new<-cbind(plt_del.t,pop681new)
head(pop681new)
dim(pop681new)
pop681new[is.na(pop681new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs


pop701new<-pop701[seq(1, nrow(pop701), 8),]
plt_del.t<-seq(1, dim(pop701new)[1] ,1)
pop701new<-cbind(plt_del.t,pop701new)
head(pop701new)
dim(pop701new)
pop701new[is.na(pop701new$sim_length),]<-NA # if sim_length==NA, then make entire row into NAs




# setwd("/Users/Jenn/Desktop/UW_Cascadia/OR_coast_AK")


numROMS<-range(pop.fin$ROMS_file, na.rm=T)[2] # 187 ROMS files ran for this particular model run
# numdays<-ceiling(numROMS/2.5) # 75 days (when particles begin to enter back into the domain)
# numdays<-numROMS/1.87 # 100 days
# numdays<-floor(numROMS/6)
# numdays2<-floor(numROMS/1.55) # 20 days (when particles leave the domain)


for (p in 1:numROMS){
	
	print(paste0("p= ", p))
	
	setwd("/Users/wongalaj/Desktop/Model/May2016_ROMSlowres")
	
	### Open netcdf files
	nc.files<-nc_open(files[p])
	# print(nc.files)
	# attributes(nc.files$var)$names # look at the names of the variables for each variable 

	### Extract v
	v_dat<-ncvar_get(nc.files, varid= 'v') # save v variable as object
	fillvalue_v<-ncatt_get(nc.files, 'v', '_FillValue') # replace the NAs with 'NA'
	v_dat[v_dat==fillvalue_v]<-NA
	v_avg<-0.5*(v_dat[2:87,,,] + v_dat[1:86,,,]) # only surface layer

	### Extract u
	u_dat<-ncvar_get(nc.files, varid= 'u') # save u variable as object
	# dim(u_dat) #  86 356  40  12 (each object is 5 dimensions that has: 	u, lon, lat, depth and time (hour))
	fillvalue_u<-ncatt_get(nc.files, 'u', "_FillValue") # replace the NAs with 'NA'
	u_dat[u_dat==fillvalue_u]<-NA 
	u_avg<-0.5*(u_dat[,2:356,,] + u_dat[,1:355,,]) # only surface layer

	### Extract time/dates
	time_dat<-ncvar_get(nc.files, varid= 'ocean_time') # save time variable as an object; units= seconds since 01-01-2005 (gregorian calendar)
	# dim(time_dat) # 12
	date<-as.Date(as.POSIXct(time_dat, origin= "2005-01-01 00:00:00", tz= "GMT"))
	datehour2<-as.POSIXct(time_dat, origin= "2005-01-01 00:00:00", tz= "GMT")
	datehour<-as.character(datehour2)
	
	for(r in (1:dim(u_avg)[4])){
		
		setwd("/Users/wongalaj/Desktop/Model/figs_for_gif")
		
		### select only surface currents
		
		uavg2<-u_avg[,,40,r]  
		vavg2<-v_avg[,,40,r]
	
		png(paste("PTMgif_20191107_surface_", r+(dim(u_avg)[4]*(p-1)), ".png", sep= "" ), 	width=1350 , height= 2350, res= 150)
	
		print(paste0("r= ", r))

	
		plot(OR_bathy, image = T, land = T, axes = T, lwd=0.1, 	deep=-3000,shallow=0, step=100, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)),xlab=expression(paste("Longitude ("^o,'W)')), ylab= expression(paste("Latitude ("^o,'N)')), main= paste("released on 05/01/2016 for 31 days", r+(dim(u_avg)[4]*(p-1)), p, sep="_"))  
	
		# plot outline of continent
		plot(OR_bathy, n = 1, lwd = 0.5, add = TRUE) # plot outline of OR coast		
		
		
		#plot current vectors
		quiver2D(u= uavg2, v= vavg2, x=lon_rho[,,40], y=lat_rho[,,40], scale= 1, by= c(5,10), type= "triangle", add= T)

		# plot MR lines
		## Otter Rock (NESTED)
		# lines(otter$Longitude, otter$Latitude, col= "red", lwd= 1.5)

		# ## Cape Perpetua (NESTED)
		# lines(NMPA$Longitude, NMPA$Latitude, col= "red", lwd= 1.5)
		# lines(capeMR$Longitude, capeMR$Latitude, col= "red", lwd= 1.5) 
		# lines(seabird$Longitude, seabird$Latitude, col= "red", lwd= 1.5)
		# lines(sempa$Longitude, sempa$Latitude, col= "red", lwd= 1.5)

		# ## Redfish Rocks (NESTED)
		# lines(RR$Longitude, RR$Latitude, col= "red", lwd= 1.5)

		# ## Cape Falcon
		# lines(cfmr$Longitude, cfmr$Latitude, col= "red", lwd= 1.5)
		# lines(cfmpa$Longitude, cfmpa$Latitude, col= "red", lwd= 1.5)

		# ## Cascade Head
		# lines(chnmpa$Longitude, chnmpa$Latitude, col= "red", lwd= 1.5)
		# lines(chwmpa$Longitude, chwmpa$Latitude, col= "red", lwd= 1.5)
		# lines(chmr$Longitude, chmr$Latitude, col= "red", lwd= 1.5)
		# lines(chsmpa$Longitude, chsmpa$Latitude, col= "red", lwd= 1.5)

		## MR labels
		# Cape Blanco
		points(-124.52, 42.836294, pch= 21, cex= 1.2, bg= "darkmagenta", col="black")
		text(-124.15, 42.836294, "Cape Blanco", cex=1.3)
		# Redfish Rocks
		text(-124.02, 42.72,"Redfish Rock MR", cex=1.3)
		# points(-124.463736, 42.740, pch= 17,cex= 1.5, col= "coral3")
		# Newport 
		points(-124.04, 44.621265, pch= 21, cex= 1.2, bg= "darkmagenta", col="black")
		text(-123.81, 44.621745, "Newport", cex=1.3)
		# Cape Perpetua
		text(-123.65, 44.247974, "Cape Perpetua MR", cex=1.3)
		# points(-124.10, 44.247974, pch= 17, cex= 1.5, col= "coral3")
		# Otter Rock 
		# points(-124.02, 44.744694, pch= 17, cex= 1.5, col= "coral3")
		text(-123.69, 44.77, "Otter Rock MR", cex=1.3)
		# Cascade Head 
		# points(-123.95, 45.012449, pch= 17, cex= 1.5, col= "coral3")
		text(-123.54, 45.012449, "Cascade Head MR", cex=1.3)
		# Astoria
		points(-123.611840, 46.237110, pch= 21, cex= 1.2, bg= "darkmagenta", col="black")
		text(-123.6, 46.12, "Astoria", cex=1.3)
		# Cape Falcon
		# points(-123.89, 45.766817, pch= 17, cex= 1.5, col= "coral3")
		text(-123.54, 45.766817, "Cape Falcon MR", cex=1.3)

		## Legend and scale
		legend(x=-124.1, y= 42.2, legend= c("Other Locations", "Marine Reserves"),  pch= c(21, NA), lty= c(NA,1), col= c("black", "red"), box.col= "black" ,pt.bg= c("darkmagenta"), bg="white", pt.cex= 1.2, lwd=1.5)
		scaleBathy(OR_bathy, x= "bottomright", deg= 0.5)

	
		### Plot surface particles at each timestep
		points(pop1$lon[pop1new$del.t==r+(dim(u_avg)[4]*(p-1))], pop1$lat[pop1new$del.t==r+(dim(u_avg)[4]*(p-1))],col='black',cex= 1.75, pch= 19)
		points(pop21$lon[pop21new$del.t==r+(dim(u_avg)[4]*(p-1))], pop21$lat[pop21new$del.t==r+(dim(u_avg)[4]*(p-1))],col="darkred",cex= 1.75, pch= 19)
		points(pop41new$lon[pop41new$del.t==r+(dim(u_avg)[4]*(p-1))], pop41new$lat[pop41new$del.t==r+(dim(u_avg)[4]*(p-1))],col= 'darkorange3',cex=1.75, pch= 19)
		points(pop61new$lon[pop61new$del.t==r+(dim(u_avg)[4]*(p-1))], pop61new$lat[pop61new$del.t ==r+(dim(u_avg)[4]*(p-1))],col='goldenrod4',cex= 1.75, pch= 19)
		points(pop81new$lon[pop81new$del.t==r+(dim(u_avg)[4]*(p-1))], pop81new$lat[pop81new$del.t==r+(dim(u_avg)[4]*(p-1))],col='chartreuse3',cex= 1.75, pch= 19)
		points(pop101new$lon[pop101new$del.t==r+(dim(u_avg)[4]*(p-1))], pop101new$lat[pop101new$del.t==r+(dim(u_avg)[4]*(p-1))],col='forestgreen',cex= 1.75, pch= 19)
		points(pop121new$lon[pop121new$del.t==r+(dim(u_avg)[4]*(p-1))], pop121new$lat[pop121new$del.t==r+(dim(u_avg)[4]*(p-1))],col='darkslategray4',cex= 1.75, pch= 19)
		points(pop141new$lon[pop141new$del.t==r+(dim(u_avg)[4]*(p-1))], pop141new$lat[pop141new$del.t==r+(dim(u_avg)[4]*(p-1))],col='cadetblue3',cex= 1.75, pch= 19)
		# points(pop161new$lon[pop161new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop161new$lat[pop161new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkorchid3',cex= 1.75, pch= 19)
		# points(pop181new$lon[pop181new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop181new$lat[pop181new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='brown',cex= 1.75, pch= 19)
		# points(pop201new$lon[pop201new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop201new$lat[pop201new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='gray',cex= 1.75, pch= 19)
		# points(pop221new$lon[pop221new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop221new$lat[pop221new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='black',cex= 1.75, pch= 19)
		# points(pop241new$lon[pop241new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop241new$lat[pop241new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col="darkred",cex= 1.75, pch= 19)
		# points(pop261new$lon[pop261new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop261new$lat[pop261new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkorange3',cex= 1.75, pch= 19)
		# points(pop281new$lon[pop281new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop281new$lat[pop281new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='goldenrod4',cex= 1.75, pch= 19)
		# points(pop301new$lon[pop301new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop301new$lat[pop301new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='chartreuse3',cex= 1.75, pch= 19)
		# points(pop321new$lon[pop321new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop321new$lat[pop321new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='forestgreen',cex= 1.75, pch= 19)
		# points(pop341new$lon[pop341new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop341new$lat[pop341new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col= 'darkslategray4',cex= 1.75, pch= 19)
	
	
	### particles released at depth
	
		# points(pop361new$lon[pop361new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop361new$lat[pop361new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='black',cex= 1.75, pch= 19)
		# points(pop381new$lon[pop381new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop381new$lat[pop381new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkred',cex= 1.75, pch= 19)
		# points(pop401new$lon[pop401new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop401new$lat[pop401new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkorange3',cex= 1.75, pch= 19)
		# points(pop421new$lon[pop421new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop421new$lat[pop421new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='goldenrod4',cex= 1.75, pch= 19)	
		# points(pop441new$lon[pop441new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop441new$lat[pop441new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='chartreuse3',cex= 1.75, pch= 19)
		# points(pop461new$lon[pop461new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop461new$lat[pop461new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='forestgreen',cex= 1.75, pch= 19)
		# points(pop481new$lon[pop481new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop481new$lat[pop481new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkslategray4',cex= 1.75, pch= 19)
		# points(pop501new$lon[pop501new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop501new$lat[pop501new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='cadetblue3',cex= 1.75, pch= 19)
		# points(pop521new$lon[pop521new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop521new$lat[pop521new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkorchid3',cex= 1.75, pch= 19)
		# points(pop541new$lon[pop541new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop541new$lat[pop541new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='brown',cex= 1.75, pch= 19)
		# points(pop561new$lon[pop561new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop561new$lat[pop561new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='gray',cex= 1.75, pch= 19)
		# points(pop581new$lon[pop581new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop581new$lat[pop581new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='black',cex= 1.75, pch= 19)
		# points(pop601new$lon[pop601new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop601new$lat[pop601new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkred',cex= 1.75, pch= 19)
		# points(pop621new$lon[pop621new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop621new$lat[pop621new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkorange3',cex= 1.75, pch= 19)
		# points(pop641new$lon[pop641new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop641new$lat[pop641new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='goldenrod4',cex= 1.75, pch= 19)
		# points(pop661new$lon[pop661new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop661new$lat[pop661new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='chartreuse3',cex= 1.75, pch= 19)
		# points(pop681new$lon[pop681new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop681new$lat[pop681new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='forestgreen',cex= 1.75, pch= 19)
		# points(pop701new$lon[pop701new$plt_del.t==r+(dim(u_avg)[4]*(p-1))], pop701new$lat[pop701new$plt_del.t==r+(dim(u_avg)[4]*(p-1))],col='darkslategray4',cex= 1.75, pch= 19)		
		
	### Close device
	dev.off()
	
	}

}



##########################################
## Contour plots of current velocities with vector plots on top to show direction and trajectories


#### Plot of depths particle travelled over time. 

plot(pop.fin$timestep_day[pop.fin$ID==1], pop.fin$lat[pop.fin$ID==1], xlab= "time step (days)", ylab= "depth (m)", type= 'l',lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.360"], pop.fin$lat[pop.fin$ID=="1.360"], type= 'l', col= 'darkred', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.720"], pop.fin$lat[pop.fin$ID=="1.720"], type= 'l', col= 'darkorange3', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.1080"], pop.fin$depth[pop.fin$ID=="1.1080"], type= 'l', col= 'darkorange3', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID== "1.1440"], pop.fin$depth[pop.fin$ID== "1.1440"], type= 'l', col= 'goldenrod4', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID== "1.1800"], pop.fin$depth[pop.fin$ID=="1.1800"], type= 'l', col= 'chartreuse3', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID== "1.2160"], pop.fin$depth[pop.fin$ID== "1.2160"], type= 'l', col= 'forestgreen', lwd= 3)
lines(pop.fin$timestep_day[pop.fin$ID=="1.2520"], pop.fin$depth[pop.fin$ID=="1.2520"], type= 'l', col= 'darkslategray4', lwd= 3)
# lines(pop.fin$timestep_day[pop.fin$ID==8], pop.fin$depth[pop.fin$ID==8], type= 'l', col= 'cadetblue3', lwd= 3)
# lines(pop.fin$timestep_day[pop.fin$ID==9], pop.fin$depth[pop.fin$ID==9], type= 'l', col= 'darkorchid3', lwd= 3)
# lines(pop.fin$timestep_day[pop.fin$ID==10], pop.fin$depth[pop.fin$ID==10], type= 'l', col= 'brown', lwd= 3)
# lines(pop.fin$timestep_day[pop.fin$ID==11], pop.fin$depth[pop.fin$ID==11], type= 'l', col= 'gray', lwd= 3)

legend(x= "bottomleft", y=NULL, legend= seq(1,4,1), lty= 1, col = c('black', "darkred", 'darkorange3', 'goldenrod4', 'chartreuse3', 'forestgreen', 'darkslategray4', 'cadetblue3', 'darkorchid3', 'brown', 'gray'), title= 'Particle IDs', lwd= 3)

plot(pop.fin$timestep_day, pop.fin$depth, type= 'l', lwd= 3, xlab= 'timestep (days)', ylab= "depth (m)", main= "Time vs. Depth")
plot(pop.fin$timestep_day, pop.fin$lat, type= 'l', lwd= 3,  xlab= 'timestep (days)', ylab= expression(paste("Latitude ("^o,'N)')), main= 'Time vs. lat')


##########################################

### Plot entire particles trajectory

plot(OR_bathy, image = T, land = T, axes = T, lwd=0.1, 	deep=-3000,shallow=0, step=100, bpal = list(c(0, max(OR_bathy), "grey"), c(min(OR_bathy),0,blues)),xlab=expression(paste("Longitude ("^o,'W)')), ylab= expression(paste("Latitude ("^o,'N)'))) # , main= paste("released on 05/01/2016 for 31 days", r+(dim(u_avg)[4]*(p-1)), p, sep="_"))  

	
# plot outline of continent
plot(OR_bathy, n = 1, lwd = 0.5, add = TRUE) # plot outline of OR coast
	
points(pop1$lon, pop1$lat, pch=19)
points(pop21$lon, pop21$lat, col= 'red', pch= 19)
points(pop41$lon, pop41$lat, col= 'green', pch= 19)	
points(pop61$lon, pop61$lat, col= 'blue', pch= 19)	
points(pop81$lon, pop81$lat, col= 'purple', pch= 19)	
points(pop101$lon, pop101$lat, col= 'pink', pch= 19)	
points(pop121$lon, pop121$lat, col= 'brown', pch= 19)	
points(pop141$lon, pop141$lat, col= 'gray', pch= 19)	


plot(pop1$lat, pop1$lon)

plot(pop.fin$lon, pop.fin$lat, type="l")


### Pau :)

