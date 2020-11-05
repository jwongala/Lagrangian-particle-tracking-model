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








