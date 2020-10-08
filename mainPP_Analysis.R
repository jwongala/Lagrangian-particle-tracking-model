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
### initial variables or constants

# Radius of the Earth  
R<-6378137 # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method


####################################################################
### post processing of raw LPT model output


####################################################################
### create dataframe to plot with 

# input to create mod function would only be pop.roms (I think?)

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

















