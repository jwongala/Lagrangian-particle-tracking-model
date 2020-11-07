### The purpose of this script is to create maps of the LPT model output using bathymetry maps.
### Jennifer Wong-Ala
### 20201106

####################################################################
### load files 

### Bathymetry (this will be used if the NOAA bathy server is down)

setwd('/Users/jennifer/Desktop/')
str_name<-"large_domain.tiff"
OR_bathy<-raster(str_name)

library(graphics)
OR_bathy2<-as.bathy(OR_bathy)

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

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(50,1000,2000),labcex=0.7,add=T,col='black', lwd= 0.07) 

points(init_par$lon, init_par$lat, col='coral3', cex=0.5, pch=20)

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

