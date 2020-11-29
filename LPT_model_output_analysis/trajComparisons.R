### The purpose of this script is to create maps to compare and contrast the 
### depth of particles over time and space, and then vertical velocity over time and space

## Jennifer Wong-Ala
## 20201110

####################################################################
### clear the workspace

rm(list=ls())


####################################################################
### load libraries

library(marmap)
library(raster)
# library(rgdal)
library(maps)
library(rgeos)
library(plotly)
library(fields)
library(plot3D)
library(lattice)

####################################################################
### load in data

## load in LPT model output

setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/model_output_chapter1&2/used_high_res_ROMS/2018')

load("LPT_data_pp_07.03.2020_JWA_6hour_2018_high.RData") 

	# LPT_data_pp_06.13.2020_JWA_6hour_2016_high.RData
	# LPT_data_pp_06.27.2020_JWA_6hour_2017_high.RData
	# LPT_data_pp_07.03.2020_JWA_6hour_2018_high.RData 


## load in ROMS bathy and lat and lon

setwd('/Volumes/TOSHIBA EXT/Desktop_backup/bathy_data')

load('highres_ROMS_bathy.RData')
load('lat_hires_bathy.RData')
load('lon_hires_bathy.RData')


## load in MR polygons for mapping

MR_directory<-"/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"
transform_MR<-"+proj=longlat +lat_1=43 +lat_2=48 +lat_0=41 +lon_0=-117 +x_0=700000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 "

MR<-readOGR(dsn=path.expand("/Volumes/TOSHIBA EXT/Desktop_backup/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"))
# class(MR) # SpatialPolygonsDataFrame
# crs(MR)
# extent(MR)
crsmerc=CRS(transform_MR) # this transforms the model 
MR_transformed<-spTransform(MR, CRS=crsmerc)


### Bathymetry

setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/other_data_for_plotting')

bathy.dat<-read.table('etopo1_bedrock.xyz', sep='')
names(bathy.dat)<-c('lon', "lat", 'depth')

bathy.dat$depth[bathy.dat$depth>0]<-NA # Avoid points above water
# head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]


title<-'2018 ebb/spring tide'


####################################################################
### load files 

### Bathymetry (this will be used if the NOAA bathy server is down)

setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/other_data_for_plotting')
str_name<-"250m_ROMS.tiff"
OR_bathy<-raster(str_name)

library(graphics)
OR_bathy2<-as.bathy(OR_bathy)

OR_bathy2<-getNOAA.bathy(lon1=-125, lon2=-123.5, lat1=42, lat2=46) # make sure to put longitude as negative


####################################################################
### separate data based depth (shallow pars <= 200 m; deep pars >= 250 m)

dim(data_use)
# head(data_use)

# range(data_use$depth)

deep_pars2<-subset(data_use, depth<=-250)
dim(deep_pars2)

deep_IDs<-unique(deep_pars2$ID) # how many deep particles? and what are the IDs?

	# 2016 = 31 deep particles
	# 2017 = 13
	# 2018 = 17

deep_pars<-subset(data_use, ID %in% deep_IDs)
dim(deep_pars) # 103192     16

deep_init2<-subset(deep_pars, timestep==1)
dim(deep_init2)

deep_init<-deep_init2[!duplicated(deep_init2$ID), ] # remove duplicated rows taht have multiple ID numbers
dim(deep_init)

####################################################################
### subset out every nth 

deep_pars2<-deep_pars[order(deep_pars$ID),] # order df by the ID column
dim(deep_pars2)

deep_pars3<-deep_pars2[seq(1,nrow(deep_pars2), 40), ] # select out and plot only every 8th particle location 
dim(deep_pars3)


####################################################################
### plot of particle initial locations

blues<-c("lightsteelblue4",  "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")  


plot(OR_bathy2, image = T, land = T, axes = T, lwd=1, deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy2), "grey"), c(min(OR_bathy2),0,blues)),xlab=expression(paste("Longitude ("^o,'W)')), ylab= expression(paste("Latitude ("^o,'N)')), cex.lab=1.3, cex.axis=1.2, xlim= c(-125.1,-123.5))  #, xlim= c(-125.5,-123.5), ylim=c(40.2, 47) 
points(deep_pars$lon[1:31], deep_pars$lat[1:31], col='black', pch=21, bg='green3', cex=0.9) # plot initial location of deep particles 

# plot specific isobaths
plot(OR_bathy2, deep=-50,shallow=-50, step=50, lwd = 0.3, add = TRUE, col='black', lty=2) # 50 m isobath
plot(OR_bathy2, deep=-100,shallow=-100, step=100, lwd = 0.3, add = TRUE, col='black', lty=2) # 100 m isobath
plot(OR_bathy2, deep=-500,shallow=-500, step=500, lwd = 0.3, add = TRUE, col='black', lty=2) # 500 m isobath
plot(OR_bathy2, deep=-1000,shallow=-1000, step=1000, lwd = 0.3, add = TRUE, col='black', lty=2) # 1000 m isobath

text(-124.4, 42.85, "CB", cex=1.2)
text(-123.95, 44.6368, 'N', cex=1.2)

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons


####################################################################
### plot particle trajectory where z changes color with depth 
# options(digits=3)
a<-seq(-750,0, length.out=10) # for depth cut 

	# seq(-0.01,0.01,by=0.002) # for w_vel cut

b<-length(a) # how many deep particles

# a2<-format(seq(-1.0,1.0,by=0.2), digits=3)

plot(OR_bathy2, image = T, land = T, axes = T, lwd=1, deep=-10,shallow=10, step=10, bpal = list(c(0, max(OR_bathy2), "grey"), c(min(OR_bathy2),0,'white')),xlab=expression(paste("Longitude ("^o,'W)')), ylab= expression(paste("Latitude ("^o,'N)')), cex.lab=1, cex.axis=1, xlim= c(-125.1,-123.5), main = title)  #, xlim= c(-125.5,-123.5), ylim=c(40.2, 47) 

rbPal<-colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue')) # creating a continuous color
	
	# 2016 = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue')) # creating a continuous color palette
	# 2017 = colorRampPalette(c('orange', 'yellow', 'green', 'blue')) # creating a continuous color palette
	# 2018 = colorRampPalette(c('orange', 'yellow', 'green', 'blue')) # creating a continuous color palette

deep_pars3$Col<-rbPal(b)[as.numeric(cut(deep_pars3$depth, breaks=a))] # add a column of color values based on the depth values

# leg<-c('(-734 to -660)', "(-660 to -587)", '(-587 to -513)', '(-440 to -367)', '(-367 to -293)', '(-293 to -220)', '(-220 to -147)', '(-147 to -37.3)', '(-73.3 to 0.73)') 

plot(deep_pars3$lon, deep_pars3$lat, col= deep_pars3$ID, pch=19, cex=0.4) # plot initial location of deep particles 

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(0,50,100,500,1000), add=T,labcex=0.7,col='black', lwd= 0.5) 

legend("bottomright",title='Depth range (m)',legend=levels(cut(deep_pars3$depth, breaks=a)), col=rbPal(b),pch=20, cex=0.6, bty='n') # [2:9] Vertical velocity (cm/s) deep_pars3$w_vel*100, breaks=a*100

# plot specific isobaths
plot(OR_bathy2, deep=-50,shallow=-50, step=50, lwd = 0.4, add = TRUE, col='black', lty=2) # 50 m isobath
plot(OR_bathy2, deep=-100,shallow=-100, step=100, lwd = 0.4, add = TRUE, col='black', lty=2) # 100 m isobath
plot(OR_bathy2, deep=-500,shallow=-500, step=500, lwd = 0.4, add = TRUE, col='black', lty=2) # 500 m isobath
plot(OR_bathy2, deep=-1000,shallow=-1000, step=1000, lwd = 0.4, add = TRUE, col='black', lty=2) # 1000 m isobath


text(-124.4, 42.85, "CB", cex=1.2)
text(-123.95, 44.6368, 'N', cex=1.2)

lines(MR_transformed, col= 'black', lwd= 1.5, lty= 1) # plot MR polygons


points(deep_init$lon, deep_init$lat, col='black', pch=21, bg='gray80', cex=0.8) # plot initial location of deep particles 


####################################################################
### plot particle trajectory where z changes color with depth 


h_hi2<-(h_hi)*-1 # 

h_hi2[h_hi2>0]<-NA

# wireframe(h_hi2, shade=T, aspect = c(1,0.2), ylab='Latitude', xlab='Longitude', zlab='Depth')

IDs2<-unique(deep_pars3$ID)

IDs<-sample(IDs2, size=10, replace=F)

y<-lon_hi
x<-lat_hi
z<-h_hi2


-125.1000 -124.8190 -124.5379 -124.2569 -123.9759

p<-plot_ly(x=x,y=y, z=z, type='surface') %>%

	layout(title= "2018 250 m ROMS - ebb/spring tide", scene=list(xaxis=list(title='Latitude'), yaxis=list(title='Longitude'), zaxis=list(title='Depth (m)'))) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[1]], y=deep_pars3$lon[deep_pars3$ID==IDs[1]], z=deep_pars3$depth[deep_pars3$ID==IDs[1]], mode='markers', type='scatter3d', marker=list(size=3, color='red', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[2]], y=deep_pars3$lon[deep_pars3$ID==IDs[2]], z=deep_pars3$depth[deep_pars3$ID==IDs[2]], mode='markers', type='scatter3d', marker=list(size=3, color='orange', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[3]], y=deep_pars3$lon[deep_pars3$ID==IDs[3]], z=deep_pars3$depth[deep_pars3$ID==IDs[3]], mode='markers', type='scatter3d', marker=list(size=3, color='yellow', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[4]], y=deep_pars3$lon[deep_pars3$ID==IDs[4]], z=deep_pars3$depth[deep_pars3$ID==IDs[4]], mode='markers', type='scatter3d', marker=list(size=3, color='green', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[5]], y=deep_pars3$lon[deep_pars3$ID==IDs[5]], z=deep_pars3$depth[deep_pars3$ID==IDs[5]], mode='markers', type='scatter3d', marker=list(size=3, color='blue', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[6]], y=deep_pars3$lon[deep_pars3$ID==IDs[6]], z=deep_pars3$depth[deep_pars3$ID==IDs[6]], mode='markers', type='scatter3d', marker=list(size=3, color='purple', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[7]], y=deep_pars3$lon[deep_pars3$ID==IDs[7]], z=deep_pars3$depth[deep_pars3$ID==IDs[7]], mode='markers', type='scatter3d', marker=list(size=3, color='pink', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[8]], y=deep_pars3$lon[deep_pars3$ID==IDs[8]], z=deep_pars3$depth[deep_pars3$ID==IDs[8]], mode='markers', type='scatter3d', marker=list(size=3, color='gray', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[9]], y=deep_pars3$lon[deep_pars3$ID==IDs[9]], z=deep_pars3$depth[deep_pars3$ID==IDs[9]], mode='markers', type='scatter3d', marker=list(size=3, color='brown', symbol=104)) %>%
	add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[10]], y=deep_pars3$lon[deep_pars3$ID==IDs[10]], z=deep_pars3$depth[deep_pars3$ID==IDs[10]], mode='markers', type='scatter3d', marker=list(size=3, color='white', symbol=104)) # %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[11]], y=deep_pars3$lon[deep_pars3$ID==IDs[11]], z=deep_pars3$depth[deep_pars3$ID==IDs[11]], mode='markers', type='scatter3d', marker=list(size=3, color='darkgolden', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[12]], y=deep_pars3$lon[deep_pars3$ID==IDs[12]], z=deep_pars3$depth[deep_pars3$ID==IDs[12]], mode='markers', type='scatter3d', marker=list(size=3, color='cadetblue', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[13]], y=deep_pars3$lon[deep_pars3$ID==IDs[13]], z=deep_pars3$depth[deep_pars3$ID==IDs[13]], mode='markers', type='scatter3d', marker=list(size=3, color='aquamarine', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[14]], y=deep_pars3$lon[deep_pars3$ID==IDs[14]], z=deep_pars3$depth[deep_pars3$ID==IDs[14]], mode='markers', type='scatter3d', marker=list(size=3, color='bisque3', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[15]], y=deep_pars3$lon[deep_pars3$ID==IDs[15]], z=deep_pars3$depth[deep_pars3$ID==IDs[15]], mode='markers', type='scatter3d', marker=list(size=3, color='burlywood3', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[16]], y=deep_pars3$lon[deep_pars3$ID==IDs[16]], z=deep_pars3$depth[deep_pars3$ID==IDs[16]], mode='markers', type='scatter3d', marker=list(size=3, color='darkgoldenrod', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[17]], y=deep_pars3$lon[deep_pars3$ID==IDs[17]], z=deep_pars3$depth[deep_pars3$ID==IDs[17]], mode='markers', type='scatter3d', marker=list(size=3, color='cornflowerblue', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[18]], y=deep_pars3$lon[deep_pars3$ID==IDs[18]], z=deep_pars3$depth[deep_pars3$ID==IDs[18]], mode='markers', type='scatter3d', marker=list(size=3, color='darkcyan', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[19]], y=deep_pars3$lon[deep_pars3$ID==IDs[19]], z=deep_pars3$depth[deep_pars3$ID==IDs[19]], mode='markers', type='scatter3d', marker=list(size=3, color='darkgreen', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[20]], y=deep_pars3$lon[deep_pars3$ID==IDs[20]], z=deep_pars3$depth[deep_pars3$ID==IDs[20]], mode='markers', type='scatter3d', marker=list(size=3, color='darkorchid4', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[21]], y=deep_pars3$lon[deep_pars3$ID==IDs[21]], z=deep_pars3$depth[deep_pars3$ID==IDs[21]], mode='markers', type='scatter3d', marker=list(size=3, color='darkslategray', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[22]], y=deep_pars3$lon[deep_pars3$ID==IDs[22]], z=deep_pars3$depth[deep_pars3$ID==IDs[22]], mode='markers', type='scatter3d', marker=list(size=3, color='gold', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[23]], y=deep_pars3$lon[deep_pars3$ID==IDs[23]], z=deep_pars3$depth[deep_pars3$ID==IDs[23]], mode='markers', type='scatter3d', marker=list(size=3, color='firebrick', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[24]], y=deep_pars3$lon[deep_pars3$ID==IDs[24]], z=deep_pars3$depth[deep_pars3$ID==IDs[24]], mode='markers', type='scatter3d', marker=list(size=3, color='darkolivegreen', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[25]], y=deep_pars3$lon[deep_pars3$ID==IDs[25]], z=deep_pars3$depth[deep_pars3$ID==IDs[25]], mode='markers', type='scatter3d', marker=list(size=3, color='darkviolet', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[26]], y=deep_pars3$lon[deep_pars3$ID==IDs[26]], z=deep_pars3$depth[deep_pars3$ID==IDs[26]], mode='markers', type='scatter3d', marker=list(size=3, color='deeppink2', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[27]], y=deep_pars3$lon[deep_pars3$ID==IDs[27]], z=deep_pars3$depth[deep_pars3$ID==IDs[27]], mode='markers', type='scatter3d', marker=list(size=3, color='darkseagreen1', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[28]], y=deep_pars3$lon[deep_pars3$ID==IDs[28]], z=deep_pars3$depth[deep_pars3$ID==IDs[28]], mode='markers', type='scatter3d', marker=list(size=3, color='darkorange', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[29]], y=deep_pars3$lon[deep_pars3$ID==IDs[29]], z=deep_pars3$depth[deep_pars3$ID==IDs[29]], mode='markers', type='scatter3d', marker=list(size=3, color='deepskyblue', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[30]], y=deep_pars3$lon[deep_pars3$ID==IDs[30]], z=deep_pars3$depth[deep_pars3$ID==IDs[30]], mode='markers', type='scatter3d', marker=list(size=3, color='darksalmon', symbol=104)) %>%
	# add_trace(x=deep_pars3$lat[deep_pars3$ID==IDs[31]], y=deep_pars3$lon[deep_pars3$ID==IDs[31]], z=deep_pars3$depth[deep_pars3$ID==IDs[31]], mode='markers', type='scatter3d', marker=list(size=3, color='chartreuse', symbol=104))

setwd('/Users/wongalaj/Desktop')

saveWidget(p, '2018_LPT.html', selfcontained=F, libdir='lib')

htmlwidgets::saveWidget(as_widget(p), "2018_index.html", selfcontained=F, libdir='lib')




















