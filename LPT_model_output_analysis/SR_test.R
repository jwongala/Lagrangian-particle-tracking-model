### SR metric

init_set<-init_par
	# init_par
	# init_surf
	# init_mid
	# init_deep

dat_set<-fin_end
	# fin_end
	# surf_end
	# mid_end
	# deep_end 

title_use<-'all par- 250m ROMS'

setwd('/Users/jennifer/Desktop')
bathy.dat<-read.table('etopo1_bedrock.xyz',sep='') # This is a big file, so we need to subsample it 
names(bathy.dat)<-c('lon','lat','depth')

bathy.loess<-loess(depth~lon*lat,span=0.02,degree=1,data=bathy.dat)
                  
summary(lm(bathy.loess$fitted~bathy.dat$depth))
                  
dat_set$bot.depth<-predict(bathy.loess,newdata=dat_set)



lat_seq<-seq(45.5,42, -0.5)

lat_list_end<-as.list(1:(length(lat_seq)-1))

lat_list_start<-as.list(1:(length(lat_seq)-1))

self_rec<-NA*1:(length(lat_seq)-1)



for(i in 1:(length(lat_seq)-1)){

	lat_list_end[[i]]<-dat_set$ID[dat_set$lat<=lat_seq[i]&dat_set$lat>lat_seq[i+1]&dat_set$bot.depth>(-100)]
	
	lat_list_start[[i]]<-init_set$ID[init_set$lat<=lat_seq[i]&init_set$lat>lat_seq[i+1]]
	
	self_rec[i]<-sum(lat_list_end[[i]] %in% lat_list_start[[i]])/length(lat_list_start[[i]])	
	
	
}

# plot(self_rec, ylim=c(0,1), type='h', lwd=10, xlab='Regions', ylab='% of particles self-recruiting', lty=1)


####################################################################
## save and load data 
setwd('/Volumes/TOSHIBA EXT/Desktop_backup/model_output_chapter1&2/SR_lowres')

# setwd('/Users/jennifer/Desktop/model_output_chapter1&2/SR_lowres')

load("SR_data_pp_06.17.2020_JWA_init_2016_low.RData")

ls()

# SR_low16_02<-SR_low16_01

# save(SR_low16_02, file='SR_data_pp_06.11.2020_JWA_init_2016_low.RData') # redo SR_low18_04



####################################################################
## take average for all model simulations

low_res_mean<-rowMeans(cbind(SR_low16_01, SR_low16_02, SR_low16_03, SR_low16_04, SR_low17_01, SR_low17_02, SR_low17_03, SR_low17_04, SR_low18_01, SR_low18_02, SR_low18_03, SR_low18_04))

hi_res_mean<-rowMeans(cbind(SR_hi16_01,SR_hi16_02,SR_hi16_03,SR_hi16_04,SR_hi17_01,SR_hi17_02,SR_hi17_03,SR_hi17_04,SR_hi18_01,SR_hi18_02,SR_hi18_03,SR_hi18_04))

hi16<-rowMeans(cbind(SR_hi16_01,SR_hi16_02,SR_hi16_03,SR_hi16_04))

hi17<-rowMeans(cbind(SR_hi17_01,SR_hi17_02,SR_hi17_03,SR_hi17_04))

hi18<-rowMeans(cbind(SR_hi18_01,SR_hi18_02,SR_hi18_03,SR_hi18_04))

low16<-rowMeans(cbind(SR_low16_01, SR_low16_02, SR_low16_03, SR_low16_04))

low17<-rowMeans(cbind(SR_low17_01, SR_low17_02, SR_low17_03, SR_low17_04))

low18<-rowMeans(cbind(SR_low18_01, SR_low18_02, SR_low18_03, SR_low18_04))

# mean2016<-rowMeans(cbind(SR_low16_01, SR_low16_02, SR_low16_03, SR_low16_04, SR_hi16_01, SR_hi16_02, SR_hi16_03, SR_hi16_04))

# mean2017<-rowMeans(cbind(SR_low17_01, SR_low17_02, SR_low17_03, SR_low17_04, SR_hi17_01, SR_hi17_02, SR_hi17_03, SR_hi17_04))

# mean2018<-rowMeans(cbind(SR_low18_01, SR_low18_02, SR_low18_03, SR_low18_04, SR_hi18_01, SR_hi18_02, SR_hi18_03, SR_hi18_04))

####################################################################
# plot
# setwd('/Users/jennifer/Desktop/Chapter2_Figures/ALL_PARS')

title<-"2018 hi res mean"

# png(paste(title_use, '.png', sep=''), height=500, width=700, res=200)
# setwd('/Users/wongalaj/Desktop')
# dev.copy(jpeg,title,height=7,width=5,res=200,units='in')
quartz(height=5,width=6)
par(mai=c(1,1,0.5,0.5))

barplot(hi18,names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), ylim=c(0,0.4), las=3, main=title, ) # #,ylim=c(0,1), xlab='Regions', ylab='% of particles self-recruiting',names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), main= title_use, cex=2, cex.lab=2, cex.main=1.5, cex.axis=2)   

# dev.off()



## panel plot

par(mfrow=c(3,2), mai=c(1,0.8,0.5,0.5))

barplot(low16,ylim=c(0,1), xlab='', ylab='',names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), main= "2 km", cex.lab=2, cex.main=2.5, cex.axis=2, cex=2, las=3)  
mtext('2016', side=3, line=-1.7, adj=0.01, cex=1.2)

barplot(hi16,ylim=c(0,1), xlab='', ylab='',names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), main= "250 m", cex.lab=2, cex.main=2.5, cex.axis=2, cex=2)  

barplot(low17,ylim=c(0,1), xlab='', ylab='% of particles self-recruiting',names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), main= "", cex.lab=2, cex.main=2.5, cex.axis=2, cex=2)  
mtext('2017', side=3, line=-1.7, adj=0.01, cex=1.2)

barplot(hi17,ylim=c(0,1), xlab='', ylab='',names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), main= "", cex.lab=2, cex.main=2.5, cex.axis=2, cex=2)  

barplot(low18,ylim=c(0,1), xlab='Regions', ylab='',names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), main= "", cex.lab=2, cex.main=2.5, cex.axis=2, cex=2)  
mtext('2018', side=3, line=-1.7, adj=0.01, cex=1.2)

barplot(hi18,ylim=c(0,1), xlab='Regions', ylab='',names.arg=c("A" ,"B" ,"C" ,"D","E" ,"F" ,"G" ), main= "", cex.lab=2, cex.main=2.5, cex.axis=2, cex=2)  

# check
# surf_end[surf_end$ID==217,]







