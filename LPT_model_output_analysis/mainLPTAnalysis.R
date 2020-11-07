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
### Load in data and clean in post processing script (ppAnalysis.R)

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

    ## unlist ppAnalysis_list
    MR_transformed<-ppAnalysis_list[[1]]
    bathy_dat<-ppAnalysis_list[[2]]
    data_use<-ppAnalysis_list[[3]]


####################################################################
### Load in data and clean in post processing script (ppAnalysis.R)

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

distretAnalysis_list<-distretAnalysis(data_use,
                                      lat1,
                                      lat2,
                                      lon1,
                                      lon2,
                                      Rad, 
                                      greatest_delt,
                                      init_delt,
                                      surfIDs,
                                      midIDs,
                                      deepIDs,
                                      nlat,
                                      nlon,
                                      empty,
                                      empty2,
                                      mpmile,
                                      CI, 
                                      numdelts,
                                      endday)    
    
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
tmp_tot<-distretAnalysis_list[[15]]
tmp2<-distretAnalysis_list[[16]]
tmp3<-distretAnalysis_list[[17]]
tmp4<-distretAnalysis_list[[18]]
    

####################################################################

foo_title<-"2018 mean FO" # title for frequency of occurrene maps

timeleft_domaindelt<-2208

blues<-c("lightsteelblue4",  "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")  

surfrange<-1:250
midrange<-251:500
deeprange<-501:750

b<-250 

bathypic_name<-"large_domain.tiff"

