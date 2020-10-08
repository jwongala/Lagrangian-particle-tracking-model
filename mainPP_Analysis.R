### Jennifer Wong-Ala 
### Date last worked on: 10/08/2020
### Purpose: The purpose of this script is to be a one stop shop for all data post processing, analysis, and creating figures for the analyses
### - this will be used to "clean" the raw LPT model output, and analyze the LPT model data and create figures from the Post_Process_LPT_Model_Output, Dispersal_and_Retention_Analyses, SR_test, ConnectivityMatrices, and Data_comparison_insitu_ROMS


####################################################################
### clear the workspace

rm(list=ls()) # clear the workspace


####################################################################
### Load libraries

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
### source functions 

setwd("/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model")

source('distance.function.R')

### I will add in more source functions later


####################################################################
### load in files needed to plot maps  

### *** NOTE TO SELF: need to move these shape files over to here so I can source it from .proj folder

### Marine Reserve Outlines

MR<-readOGR(dsn=path.expand("/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"), layer= 'MPA_MR_COMP_Boundaries_UTM10N')

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
### initial variables or constants

# Radius of the Earth  
R<-6378137 # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method


## Lat and lon boundaries

lat1<-42 
lat2<-45.24775

lon1<-(-125.10) # -126.0909
lon2<-(-123.9263)  

####################################################################
### unlist the raw output and put into a data frame where rows are particles at each time step and columns are the variables being tracked for each particles over time (e.g., lat, lon, velocity, etc.)

# input to create mod function would only be pop.roms (I)

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

# l=1
# tmp<-pop.roms[[l]]
# 
# for(j in 1:length(tmp)){
#   tmp[[j]]$del.t<-j+length(tmp)*(l-1)
#   pop.fish<-tmp[[j]]
#   pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
# }
# 
# pop.roms3<-pop.roms2 # rename tide data
# 
# rm(pop.roms2)
# 
# pop.fish<-as.data.frame(matrix(1:(out.col*out.row),ncol=out.col,nrow=out.row)*NA)
# pop.roms2<-list()
# pop.roms2[[1]]<-pop.fish
# 
# 
# # tides data 
# for(l in 1:(length(pop.roms4))){ 
#   tmp<-pop.roms4[[l]]
#   for(j in 1:length(tmp)){
#     tmp[[j]]$del.t<-j+length(tmp)*(l-1)
#     pop.fish<-tmp[[j]]
#     pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
#   }}

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

# output for source function would be only pop.mat (I think?)


####################################################################
### clean raw data set (i.e., remove initial positions and portion of a particles trajectory that left the model domain)

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











