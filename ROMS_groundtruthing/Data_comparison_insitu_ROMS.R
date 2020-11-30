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
library(ggplot2)
library(zoo)


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

setwd('/Volumes/TOSHIBA EXT/ROMS_HiRes/2016') # high res ROMS directory
	
	# setwd('/Volumes/TOSHIBA EXT/ROMS_HiRes/2016') # high res ROMS directory
	# setwd('/Volumes/TOSHIBA EXT/April2016_LowRes') # low res ROMS directory	

files<-list.files('/Volumes/TOSHIBA EXT/ROMS_HiRes/2016', pattern= '*.nc', full.names=T) # list of the ROMS files names to open adn loop through

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


MR<-readOGR(dsn=path.expand("/Volumes/TOSHIBA EXT/Desktop_backup/OregonMarineReserves_LPK_ArcGIS/commondata/gis_files"), layer= 'MPA_MR_COMP_Boundaries_UTM10N')

class(MR) # SpatialPolygonsDataFrame
crs(MR)

extent(MR)

crsmerc=CRS("+proj=longlat +lat_1=43 +lat_2=48 +lat_0=41 +lon_0=-117 +x_0=700000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ") # this transforms the model 

MR_transformed<-spTransform(MR, CRS=crsmerc)

##########################################################
### plot map to select regions 

blues <- c("lightsteelblue4", "lightsteelblue3","lightsteelblue2", "lightsteelblue1")

plot.bathy(OR_bathy2, deep=-1000, shallow=-10, step=50,image=T, land=T, lwd=0.1, bpal=list(c(0, max(OR_bathy2), 'grey'), c(min(OR_bathy2), 0, blues)), xlab=expression(paste("Longitude ("^o,'W)')), ylab=expression(paste("Latitude ("^o,'N)')))

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(50, 100), labcex=0.7,add=T,col='grey', lwd= 0.5)

lines(MR_transformed, col= 'black', lwd= 1, lty= 1) # plot MR polygons

points(-124.86, 44.1, pch=19, cex=1, col='darkred') # shelf location (HB)

points(-124.25, 44.1, pch=19, cex=1, col='darkred') # inshore location (HB)

points(-124.86, 42.82, pch=19, cex=1, col='darkred') # shelf location (CB)

points(-124.64, 42.82, pch=19, cex=1, col='darkred') # inshore location (CB)

plot(OR_bathy2, n = 1, lwd = 0.5, add = TRUE) # plot outline of OR coast		

text(-124.15, 42.836294, "Cape Blanco")
text(-123.81, 44.621745, "Newport")

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
	datehour2<-as.POSIXct(ROMS_time, origin= "2005-01-01 00:00:00", tz= "GMT") # 2005 for low res and 0001 
# 2005-01-01 00:00:00  0001-01-01 00:00:00
	b<-cbind(ROMS_w, ROMS_time)
	
	# ROMS_2018time<-rbind(ROMS_2018time, b) 
	tmp<-rbind(tmp, b)
	
} 

setwd('/Users/wongalaj/Desktop')
inshoreCB_hires18<-tmp
save(inshoreCB_hires18, file='vert_vel_hires18_inshoreCB.RData')


##########################################################
### load in post processed data and plot figures 

setwd('/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/ROMS_groundtruthing/vertical_velocity_comparisons/')

## low res data
# 2016 = 420
# 2017 = 420
# 2018 = 444

low16<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(inshoreHB_lowres16))
length(low16) # 420
hi16<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-02 00:00:00"), length.out=nrow(inshoreHB_hires16))
length(hi16)

low17<-seq(as.POSIXct("2017-03-27 02:00:00"), as.POSIXct("2017-05-02 00:00:00"), length.out=nrow(inshoreHB_lowres17))
length(low17)
hi17<-seq(as.POSIXct("2017-03-26 00:00:00"), as.POSIXct("2017-04-28 11:00:00"), length.out=nrow(inshoreHB_hires17))
length(hi17)

low18<-seq(as.POSIXct("2018-03-26 02:00:00"), as.POSIXct("2018-05-02 00:00:00"), length.out=nrow(inshoreHB_lowres18))
length(low18)
hi18<-seq(as.POSIXct("2018-03-26 00:00:00"), as.POSIXct("2018-05-01 00:00:00"), length.out=nrow(inshoreHB_hires18))
length(hi18)

##########################################################
### decompose data to look at time series components

## frequency of del.t data
freq<-12 # low res = 12
freq2<-24 # high res = 24
	
## time string 
t1<-low18 # low res
t2<-hi18 # hi res


## example of decomposing timeseries signal 
ts_dat<-ts(data=inshoreHB_lowres16[,1], frequency=freq, start= low16[1]) # convert to timeseries (ts)
summary(ts_dat) # summary of ts data
str(ts_dat) # structure of ts
dc<-decompose(ts_dat) # decompose data to show observed, trend, seasonal, and random data 
plot(dc, yax.flip=T) # plot the decomposed data
plot(dc$trend*100, lwd=5) # only plot the trend data 


## assign data to be decomposed and then plotted 
l1<-shelfHB_lowres18[,1]*100
l2<-inshoreHB_lowres18[,1]*100
l3<-shelfCB_lowres18[,1]*100
l4<-inshoreCB_lowres18[,1]*100

h1<-shelfHB_hires18[,1]*100
h2<-inshoreHB_hires18[,1]*100
h3<-shelfCB_hires18[,1]*100
h4<-inshoreCB_hires18[,1]*100


## low res
HB_shelf<-ts(data=l1, frequency=freq, start= t1[1])
length(HB_shelf)
tmp1<-decompose(HB_shelf)
HB_shelf<-tmp1$trend 

HB_inshore<-ts(data=l2, frequency=freq, start= t1[1]) # inshoreHB_lowres18*100
length(HB_inshore)
tmp2<-decompose(HB_inshore)
HB_inshore<-tmp2$trend

CB_shelf<-ts(data=l3, frequency=freq, start= t1[1])
length(CB_shelf)
tmp3<-decompose(CB_shelf)
CB_shelf<-tmp3$trend

CB_inshore<-ts(data=l4, frequency=freq, start= t1[1])
length(CB_inshore)
tmp4<-decompose(CB_inshore)
CB_inshore<-tmp4$trend


## high res
HB_shelf2<-ts(data=h1, frequency=freq, start= t2[1])
length(HB_shelf2)
tmp5<-decompose(HB_shelf2)
HB_shelf2<-tmp5$trend

HB_inshore2<-ts(data=h2, frequency=freq, start= t2[1])
length(HB_inshore2)
tmp6<-decompose(HB_inshore2)
HB_inshore2<-tmp6$trend

CB_shelf2<-ts(data=h3, frequency=freq, start= t2[1])
length(CB_shelf2)
tmp7<-decompose(CB_shelf2)
CB_shelf2<-tmp7$trend

CB_inshore2<-ts(data=h4, frequency=freq, start= t2[1])
length(CB_inshore2)
tmp8<-decompose(CB_inshore2)
CB_inshore2<-tmp8$trend


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

## ranges for all data (so they can be comparable)

# -0.0182, 0.0137 (cm/s)


## 2 km 
# 2016 = 13:length(yr_plt)
# 2017 = 1:389 
# 2018 = 1:432

## 250 m 
# 2016 = 25:length(yr_plt2)-1
# 2017 = 27:length(yr_plt2)
# 2018 = 3:length(yr_plt2)

## Time plot
a<-1:432 # 2 km
b<-3:length(yr_plt2)


quartz(width=9, height=7)
dev.copy(jpeg, paste(title_use, '.jpg', sep=''), height=7, width=9, res=200, units='in')
par(mfrow=c(2,2), mai=c(0.8, 0.8, 0.6, 0.6))


plot(yr_plt[a], HB_inshore[a], type='l', xlab='', ylab='', main='Heceta Bank', ylim=c(-0.0182, 0.0137)) # , ylim=c(-2e-03, 6.6e-04)
lines(yr_plt[a], HB_shelf[a], type='l', xlab='', ylab='', lty=2, col='red')
legend(x="topright", y=NULL, legend=c("inshore","shelf"), lty=c(1,2), col=c('black', 'red'), cex=0.7, bty='n')
mtext(text_title, side=3, line=0, adj=0)


plot(yr_plt[a], CB_inshore[a], type='l', xlab='', ylab='', main='Cape Blanco', ylim=c(-0.0182, 0.0137)) # , ylim=c(-2e-03, 6.6e-04)
lines(yr_plt[a], CB_shelf[a], type='l', xlab='', ylab='',  lty=2, col='red')


plot(yr_plt2[b], HB_inshore2[b], type='l', xlab='Time (hourly)', ylab='', main='', ylim=c(-0.0182, 0.0137)) # , ylim=c(-2e-03, 6.6e-04) 
lines(yr_plt2[b], HB_shelf2[b], type='l', xlab='', ylab='', lty=2, col='red')
mtext(text_title2, side=3, line=0, adj=0)


plot(yr_plt2[b], CB_inshore2[b], type='l', xlab='Time (hourly)', ylab='', main='', ylim=c(-0.0182, 0.0137)) # , ylim=c(-2e-03, 6.6e-04)
lines(yr_plt2[b], CB_shelf2[b], type='l', xlab='', ylab='', lty=2, col='red')
mtext('Vertical velocity (cm/s)', side=2, line=30, adj=4)


dev.off()


##########################################################
### calculate average, sd, CI, and sum of vertical velocity for the entire simulation period 

# calculations where I just change the data input

dat_tmp<-shelfCB_hires18[,1] # change data input for the analyses

avg_vert<-mean(dat_tmp)*100 # units = cm/s
sum_vert<-sum(dat_tmp)*100 # units = cm/s

sd_vert<-sd(dat_tmp)*100 # units = cm/s
se_vert<-sd_vert/sqrt(length(dat_tmp))

err_vert<-qnorm(0.95)*sd_vert/sqrt(length(dat_tmp)) # 95% CI err

l_vert<-avg_vert-err_vert # left and right side of the error
r_vert<-avg_vert+err_vert

avghi18_shCB<-avg_vert
sumhi18_shCB<-sum_vert
sehi18_shCB<-se_vert
errhi18_shCB<-err_vert


### combine data sets based on year and then spatial resolution

## low res averages, sum, sd, and error
lowres_avg<-c(avglow16_inHB, avglow16_shHB,avglow16_inCB, avglow16_shCB, avglow17_inHB, avglow17_shHB,avglow17_inCB, avglow17_shCB, avglow18_inHB, avglow18_shHB,avglow18_inCB, avglow18_shCB)
lowres_sum<-c(sumlow16_inHB, sumlow16_shHB,sumlow16_inCB, sumlow16_shCB, sumlow17_inHB, sumlow17_shHB,sumlow17_inCB, sumlow17_shCB, sumlow18_inHB, sumlow18_shHB,sumlow18_inCB, sumlow18_shCB)
lowres_se<-c(selow16_inHB, selow16_shHB,selow16_inCB, selow16_shCB, selow17_inHB, selow17_shHB,selow17_inCB, selow17_shCB, selow18_inHB, selow18_shHB,selow18_inCB, selow18_shCB)
lowres_err<-c(errlow16_inHB, errlow16_shHB,errlow16_inCB, errlow16_shCB, errlow17_inHB, errlow17_shHB,errlow17_inCB, errlow17_shCB, errlow18_inHB, errlow18_shHB,errlow18_inCB, errlow18_shCB)
lab_low<-c('inshoreHB_2016', 'shelfHB_2016', 'inshoreCB_2016', 'shelfCB_2016', 'inshoreHB_2017', 'shelfHB_2017', 'inshoreCB_2017', 'shelfCB_2017','inshoreHB_2018', 'shelfHB_2018', 'inshoreCB_2018', 'shelfCB_2018')
low_dat<-data.frame(name=lab_low, average=lowres_avg, se=lowres_se, sum= lowres_sum, error=lowres_err) # dataframe

## high res averages, sum, sd, and error
hires_avg<-c(avghi16_inHB, avghi16_shHB,avghi16_inCB, avghi16_shCB, avghi17_inHB, avghi17_shHB,avghi17_inCB, avghi17_shCB, avghi18_inHB, avghi18_shHB,avghi18_inCB, avghi18_shCB)
hires_sum<-c(sumhi16_inHB, sumhi16_shHB,sumhi16_inCB, sumhi16_shCB, sumhi17_inHB, sumhi17_shHB,sumhi17_inCB, sumhi17_shCB, sumhi18_inHB, sumhi18_shHB,sumhi18_inCB, sumhi18_shCB)
hires_se<-c(sehi16_inHB, sehi16_shHB,sehi16_inCB, sehi16_shCB, sehi17_inHB, sehi17_shHB,sehi17_inCB, sehi17_shCB, sehi18_inHB, sehi18_shHB,sehi18_inCB, sehi18_shCB)
hires_err<-c(errhi16_inHB, errhi16_shHB,errhi16_inCB, errhi16_shCB, errhi17_inHB, errhi17_shHB,errhi17_inCB, errhi17_shCB, errhi18_inHB, errhi18_shHB,errhi18_inCB, errhi18_shCB)
lab_hi<-c('inshoreHB_2016', 'shelfHB_2016', 'inshoreCB_2016', 'shelfCB_2016', 'inshoreHB_2017', 'shelfHB_2017', 'inshoreCB_2017', 'shelfCB_2017','inshoreHB_2018', 'shelfHB_2018', 'inshoreCB_2018', 'shelfCB_2018')
hi_dat<-data.frame(name=lab_hi, average=hires_avg, se=hires_se, sum= hires_sum, error=hires_err ) # datafrane 

## plots each data set 
ggplot(low_dat) + 
  geom_bar( aes(x=name, y=average), stat='identity', fill='darkblue', alpha=0.7) + 
  ylim(-0.00103, 0.00065) +
  xlab("Location") +
  ylab("average vertical velocity (cm/s)") + 
  ggtitle("2 km ROMS") + 
  theme(text=element_text(size=20, color='black'), axis.text=element_text(size=10, color='black')) +
  geom_errorbar( aes(x=name, ymin=average-se, ymax=average+se), width=0.4, colour='gray4', alpha=0.9, size=0.5)

ggplot(hi_dat) + 
  geom_bar( aes(x=name, y=average), stat='identity', fill='darkblue', alpha=0.7) + 
  ylim(-0.00103, 0.00065) + 
  xlab("Location") +
  ylab("average vertical velocity (cm/s)") + 
  ggtitle("250 m ROMS") + 
  theme(text=element_text(size=20, color='black'), axis.text=element_text(size=10, color='black')) +
  geom_errorbar( aes(x=name, ymin=average-se, ymax=average+se), width=0.4, colour='gray4', alpha=0.9, size=0.5)

