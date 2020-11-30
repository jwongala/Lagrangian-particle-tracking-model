### Jennifer Wong-Ala 
### Date last worked on: 10/08/2020
### Purpose: The purpose of this is to be used to post process and analyze all of the LPT model output script 

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
library(lattice)
library(grDevices)
library(sgeostat)
library(maps)
library(fields)


####################################################################
### initial variables or constants

# Radius of the Earth  
R<-6378137 # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method


####################################################################
### source functions 

setwd("/Volumes/TOSHIBA EXT/LPT_model_output_analysis")

source('distance.function.R')
source('ppAnalysis.R')


####################################################################
### The below function does the following (ppAnalysis.R): 

# - unlist and convert raw data to data.frame
# - remove particle trajectories that leave ROMS spatial domain)

## Input variables and directories

MR_directory<-"/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"
transform_MR<-"+proj=longlat +lat_1=43 +lat_2=48 +lat_0=41 +lon_0=-117 +x_0=700000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 "

bathy_file<-"etopo1_bedrock.xyz"

rawfile_directory<-'/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/model_output_chapter1&2/used_low_res_ROMS/2016'

file_name<-'LPT_data_pp_06.10.2020_JWA_6hour_2016_low.RData'

lat1<-42
lat2<-45.24775

lon1<-(-125.10) # -126.0909
lon2<-(-123.9263) 


ppAnalysis_list<-ppAnalysis(MR_directory=MR_directory,
                            transform_MR=transform_MR,
                            bathy_file=bathy_file,
                            rawfile_directory=rawfile_directory,
                            file_name=file_name,
                            lat1=lat1, 
                            lat2=lat2, 
                            lon1=lon1, 
                            lon2=lon2)

    ## unlist returned variables from ppAnalysis_list
    MR_transformed<-ppAnalysis_list[[1]]
    bathy_dat<-ppAnalysis_list[[2]]
    data_use<-ppAnalysis_list[[3]]
    pop.mat<-ppAnalysis_list[[4]]


####################################################################
### The below function calculates (disretAnalysis.R):

# - frequencey of occurrence (FO)
# - time spent in model domain based on initial depth layer
# - average distance traveled in the latitudinal and longitudinal direction

## Input variables and directories
lat1<-42 
lat2<-45.24775
lon1<-(-125.10) # -126.0909
lon2<-(-123.9263)  
Rad<-6378137 
greatest_delt<-1344 # only contains particle trajectories for 14 days or less (1344 del.t)
init_delt<-15 
    
surfIDs<-seq(1,250,1)
midIDs<-seq(251,500,1)
deepIDs<-seq(501,750,1)
nlat<-360/10 # dimensions of maps in lat direction
nlon<-130/10 # dimensions of maps in lon direction
    
empty<-3.24774 # conversion factors? (can't remember what these numbers are)
empty2<-1.6 
mpmile<-111195
CI<-0.95
numdelts<-96
endday<-2208    

distretAnalysis_list<-distretAnalysis(data_use=data_use,
                                      lat1=lat1,
                                      lat2=lat2,
                                      lon1=lon1,
                                      lon2=lon2,
                                      Rad=Rad, 
                                      greatest_delt=greatest_delt,
                                      init_delt=init_delt,
                                      surfIDs=surfIDs,
                                      midIDs=midIDs,
                                      deepIDs=deepIDs,
                                      nlat=nlat,
                                      nlon=nlon,
                                      empty=empty,
                                      empty2=empty2,
                                      mpmile=mpmile,
                                      CI=CI, 
                                      numdelts=numdelts,
                                      endday=endday)    
    

    ## unlist returned variables from disretAnalysis_list
    init_par<-distretAnalysis_list[[1]]
    surf_par<-distretAnalysis_list[[2]]
    mid_par<-distretAnalysis_list[[3]]
    deep_par<-distretAnalysis_list[[4]]
    surf_end<-distretAnalysis_list[[5]]
    mid_end<-distretAnalysis_list[[6]]
    deep_end<-distretAnalysis_list[[7]]
    fo_mat<-distretAnalysis_list[[8]] 
    avg_lat<-distretAnalysis_list[[9]]
    avg_lon<-distretAnalysis_list[[10]] 
    sd_lat<-distretAnalysis_list[[11]]
    sd_lon<-distretAnalysis_list[[12]]
    err_lat<-distretAnalysis_list[[13]]
    err_lon<-distretAnalysis_list[[14]]
    tmp_tot<-distretAnalysis_list[[15]] # how long particles travel in model domain (all particles)
    tmp2<-distretAnalysis_list[[16]] # how long particles travel in model domain (surface)
    tmp3<-distretAnalysis_list[[17]] # how long particles travel in model domain (middle)
    tmp4<-distretAnalysis_list[[18]] # how long particles travel in model domain (deep)
    
    
####################################################################
### plot frequency of occurrence map

foo_title<-"2018 mean FO" # title for frequency of occurrene maps
    
blues<-c("lightsteelblue4",  "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")  
    
pal<-colorRampPalette(c("white","cyan", "cyan4", "darkblue"))
    
# setwd("/Users/wongalaj/Desktop/Chapter2_Figures")
    
dev.new(width= 6, height=10, res=200)
par(mai=c(1.1,1.1,0.4,0.3))
    
image.plot(z.lon,z.lat, fo18_hi,xlab=expression(paste("Longitude ("^o,'W)')),ylab=expression(paste("Latitude ("^o,'N)')), main= "fo18_hi", col=pal(500), cex.lab=1.5,cex.axis=1.5, cex.main=1.2, zlim=c(0,7)) 
    
map("worldHires",fill=T,col="grey", add=T)
    
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(50,1000,2000),labcex=0.7,add=T,col='black', lwd= 0.1) 
    
lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons
    
text(-124.46,42.85,"CB", cex=1.5)
text(-123.97, 44.6368, 'N', cex=1.5)    
    

####################################################################
### The below function calculates (secondaryCalcs.R):

# - depth range (initial locations and for the entire simulations)
# - calculate the percent of particles advected out of the ROMS model domain
# - calculate the percent of particles kept inside of the ROMS model domain

## Input variables and directories
days<-23 # days for calculating how many particles left in model
delt_days<-960 # 10 days in del.t
tot_len<-250 # total number of particles subset is already contains (each layer has 250 at beginning of simulation)

surfrange<-1:250
midrange<-251:500
deeprange<-501:750


secondaryCalcs_list<-secondaryCalcs(pop.mat=pop.mat, 
                                    surf_par=surf_par,
                                    mid_par=mid_par,
                                    deep_par=deep_par,
                                    tmp2=tmp2, 
                                    tmp3=tmp3,
                                    tmp4=tmp4,
                                    days=days,
                                    delt_days=delt_days,
                                    tot_len=tot_len,
                                    surfrange=surfrange,
                                    midrange=midrange,
                                    deeprange=deeprange)

    ## unlist returned variables from secondaryCalcs_list
    sinitrange<-secondaryCalcs_list[[1]]
    minitrange<-secondaryCalcs_list[[2]]
    dinitrange<-secondaryCalcs_list[[3]]
    srange<-secondaryCalcs_list[[4]]
    mrange<-secondaryCalcs_list[[5]]
    drange<-secondaryCalcs_list[[6]]
    per_out2<-secondaryCalcs_list[[7]]
    per_out3<-secondaryCalcs_list[[8]]
    per_out4<-secondaryCalcs_list[[9]]
    a2<-secondaryCalcs_list[[10]]
    a3<-secondaryCalcs_list[[11]]
    a4<-secondaryCalcs_list[[12]]
    c2<-secondaryCalcs_list[[13]]
    c3<-secondaryCalcs_list[[14]]
    c4<-secondaryCalcs_list[[15]]


####################################################################
### The below function calculates (ConnectivityMatrices.R):

# - uses the LPT model to create connectivity matrices for each simulation

## Input variables and directories
delt_days2<-1344 # 10 days in del.t


ConnectivityMatrices_list<-ConnectivityMatrices(data_use=data_use,
                                                delt_days2=delt_days2)


    ## unlist returned variables from secondaryCalcs_list
    conn_matt<-ConnectivityMatrices_list[[1]]
    test<-ConnectivityMatrices_list[[2]]

    
####################################################################
### plot connectivity matrix

cols<-colorRampPalette(c('gray','blue','green', 'yellow', 'orange', 'red'))

x<-1:35

dev.new(height=8, width=8, units='inch', res=200)

image.plot(conn_mat[,nrow(conn_mat):1], xlab='', ylab='source', xaxt='n', yaxt='n', 
           col=cols(150), zlim=c(0, 0.45)) 
# need to plot matrix using df[,nrow(df):1] so that it plots in the same orientation

axis(side=2, at=seq(1,0, length.out=35), labels=x, cex.axis=0.8) # figure out how to relabel yaxis
axis(side=3, at=seq(0,1, length.out=35), labels=x, cex.axis=0.8) # figure out how to relabel xaxis

mtext('sink', side=3, line=2.5)
mtext('250 m 2018', side=1, line=1)











