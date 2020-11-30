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

