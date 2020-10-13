### Jennifer Wong-Ala
### Purpose: The purpose of this script is to compare and contrast the velocities particles experience for the entire simulation period

##########################################################
### clear workspace 

rm(list=ls())

##########################################################
### load libraries


##########################################################
### load data 
setwd("/Volumes/TOSHIBA EXT/Lagrangian-particle-tracking-model/model_output_chapter1&2/used_high_res_ROMS/2016")
load('LPT_data_pp_06.13.2020_JWA_6hour_2016_high.RData')


# simulation used is 250 m ROMS ebb tide
# LPT_06.13.2020_JWA_6hour_2016_high.RData

head(data_use)
dim(data_use)

range(data_use$ID) # 1 750 

# data_use$ID[which(data_use$depth<=-700)] # figure out which particle is the really deep one

# plot(data_use$del.t, data_use$depth)

# subset out two particles 
## one that primarily stays in the surface and the second particle is the one that goes all the way to 600 m depths. 

par_shal<-subset(data_use, ID==50)
dim(par_shal) # 3360   16
head(par_shal)

par_deep2<-subset(data_use, ID==268)
dim(par_deep2) # 3432   16
head(par_deep2)

par_deep<-par_deep2[1:nrow(par_shal),]
dim(par_deep)

date<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(par_deep))
# date1<-seq(as.POSIXct("2016-03-27 00:00:00"), as.POSIXct("2016-05-01 00:00:00"), length.out=nrow(par_shal))


sum_deep<-sum(par_deep$w_vel, na.rm=T)*100 # units = m/s (-0.4186448)

sum_shal<-sum(par_shal$w_vel, na.rm=T)*100 # units = m/s (-0.07111202)

avg_shal<-mean(par_shal$w_vel, na.rm=T)*100

avg_deep<-mean(par_deep$w_vel, na.rm=T)*100

sum_shal2<-paste('shallow w_vel sum = ', round(sum_shal, digits=2), 'cm/s' ,sep= ' ')
sum_deep2<-paste('deep w_vel sum = ', round(sum_deep, digits=2), 'cm/s',sep= ' ')

avg_shal2<-paste('shallow w_vel avg = ', round(avg_shal, digits=3),'cm/s' , sep= " ")
avg_deep2<-paste('deep w_vel avg = ', round(avg_deep, digits=3),'cm/s' , sep= " ")


dev.new()
plot(date, par_deep$depth, pch=19, cex=0.5, xlab= "Time", ylab='Depth (m)')
points(date, par_shal$depth, col='red', pch=19, cex=0.5)

text(x=as.POSIXct("2016-04-05 00:00:00"), y=-500, sum_shal2, cex=0.7)
text(x=as.POSIXct("2016-04-05 00:00:00"), y=-530, sum_deep2, cex=0.7)

text(x=as.POSIXct("2016-04-05 00:00:00"), y=-600, avg_shal2, cex=0.7)
text(x=as.POSIXct("2016-04-05 00:00:00"), y=-630, avg_deep2, cex=0.7)

