### The purpose of this script is to calculate: [insert here]
### Jennifer Wong-Ala
### 20201106

distretAnalysis<-function(data_use,lat1,lat2,lon1,lon2,Rad,greatest_delt,
                          init_delt,surfIDs,midIDs,deepIDs,nlat,nlon,empty,
                          empty2,mpmile,CI,numdelts,endday){
    
  
    ####################################################################
    ### Convert km size of MR to degrees so that I can divide up coastal ocean 
    ### into grid cells that are no bigger than the largest size of the Otter Rock MR
    
    lat1<-lat1 
    lat2<-lat2
    
    lon1<-lon1 # -126.0909
    lon2<-lon2
    
    Rad<-Rad # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method
    
    deg.lat<-lat2-lat1 # length of coastline measured in degrees 
    
    deg.lon<-lon2-lon1
    
    ## conversion of coastline length in degrees to km
    
    lat.len<-Rad*pi*(deg.lat/180)*0.001 # length of coastline (latitude in km)
    
    lon.len<-(Rad*pi*(cos(pi*(lat.len/180)))*(deg.lon/180))*0.001 # width of coastline (longitude in km)
    
    
    ####################################################################
    ### Remove part of trajectory that exceeds most least common del.t for all 
    ### particles (23 days=2208 del.t)
    
    fin_dat<-subset(data_use, del.t<=greatest_delt) # only contains particle trajectories for 14 days or less (1344 del.t)
    
    # dim(fin_dat)
    
    fin_end<-subset(data_use, del.t==greatest_delt) # only contains end location of particles 
    
    ####################################################################
    ### subset particles based on depth layer and initial location 
    
    ## init_par - init locations for all particles
    
    init_par<-subset(fin_dat, del.t==init_delt)
    
    ## surface 
    surf_par<-subset(fin_dat, ID %in% surfIDs)
    # dim(surf_par) 
    surf_end<-subset(fin_end, ID %in% surfIDs)
    
    ## middle
    mid_par<-subset(fin_dat, ID %in% midIDs)
    # dim(mid_par)
    mid_end<-subset(fin_end, ID %in% midIDs)
    
    deep_par<-subset(fin_dat, ID %in% deepIDs)
    # dim(deep_par)
    deep_end<-subset(fin_end, ID %in% deepIDs)
    
    
    ####################################################################
    ### Making Frequency of Occurrence Heat maps
    
    # change what data set using to create the plots
    
    freq_dat<-fin_dat
    	# fin_dat - all par trajectories
    	# surf_par
    	# mid_par
    	# deep_par
    	
    title_foo<-"2018 mean FO" # title to be used for figures 
    	# "2km ROMS, all particles, initial, 23 day - FoO"
    	# "2km ROMS, surface, initial, 23 day - FoO"
    	# "2km ROMS, middle, initial, 23 day - FoO"
    	
    	
    # Part 1: make a regular grid and count stations within each grid cell. Plot results
       
        nlat= nlat # lat.len/13 # 698/10 # determine resolution of grid
        nlon= nlon # lon.len/13 # 238/10
        latd=seq(lat1, lat2,length.out=nlat)
        lond=seq(lon1,lon2,length.out=nlon)
        
        
      grid.lon=data.frame(
        lon1=rep(lond[-length(lond)],(nlat-1)),
        lon2=rep(lond[-1],(nlat-1)),
        lon3=rep(lond[-1],(nlat-1)),
        lon4=rep(lond[-length(lond)],(nlat-1)))#make dataframe of just longitude
      
      grid.lat=data.frame(
        lat1=sort(rep(latd[-length(latd)],(nlon-1))),
        lat2=sort(rep(latd[-length(latd)],(nlon-1))),
        lat3=sort(rep(latd[-1],(nlon-1))),
        lat4=sort(rep(latd[-1],(nlon-1)))) # lat dataframe
       
      ### sort init_pars data
      dev.new(width=6,height=10)
      plot(freq_dat$lon, freq_dat$lat,pch='.',ylim=c(lat1, lat2),xlim=c(lon1,lon2))
      # points(maynew2$lon, maynew2$lat, col='red', pch= '.')
      n.stations=NA*(1:nrow(grid.lon))
       
      for(i in 1:length(n.stations)){
        print(i)
        tmp=in.chull(freq_dat$lon, freq_dat$lat, grid.lon[i,], grid.lat[i,])
        n.stations[i]=sum(tmp)#This decides what goes into each grid pixel (may16_3$ID*tmp)
        points(freq_dat$lon[tmp], freq_dat$lat[tmp],col=i,pch=16)
        polygon(grid.lon[i,],grid.lat[i,])
      }
      map("worldHires",fill=T,col="grey",add=T)
      
    num_ts<-nrow(freq_dat)
    
    
    ## creating grid to plot for heatmap 
    z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2 # 
    z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2 # 
    z.mat<-((matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F))/num_ts)*100 # /uniID_nov16 # units are percentage (%)
    fo_mat<-z.mat
    # range(z.mat, na.rm=T)
    
    
    ####################################################################
    ### Averaged distance travelled
    
    fin_trav<-fin_dat # entire data set to be looking at 
    init_trav<-init_par # only initial locations of data set
    
    end_trav<-fin_end # only end locations of data set
    
    empty<-NA*1:length(unique(end_trav$ID))
    empty2<-NA*1:length(unique(end_trav$ID))
    par_ids<-unique(end_trav$ID)
    
    for(i in 1:length(empty)){
    	begin<-subset(init_trav, ID %in% par_ids[i])
    	begin_lat<-begin[1,1]
    	begin_lon<-begin[1,2]
    	
    	end<-subset(end_trav, ID %in% par_ids[i]) 
    	end_lat<-end[1,1] 
    	end_lon<-end[1,2] 
    	
    	dlat<-begin_lat-end_lat
    	dlon<-begin_lon-end_lon
    	
    	empty[i]<-dlat # distance in degrees
    	empty2[i]<-dlon 
    }
    
    empty<-empty 
    empty2<-empty2 
    
    # convert degrees to km (negative= travel north, positive=travel south) 
    
    test_lat<-(pi*empty/180)
    
    shrink<-cos(test_lat)
    
    
    test_lon<-mpmile*sqrt((empty2*shrink)^2 + (test_lat)^2) # meters
    
    dist_lon<-test_lon/1000 # lon dist travelled in km 
    dist_lat<-Rad*pi*(empty/180)*0.001 # convert to m 
    
    tot_dist<-as.data.frame(cbind(par_ids, dist_lat, dist_lon))
    
    ##
    avg_lat<-mean(tot_dist$dist_lat)
    sd_lat<-sd(tot_dist$dist_lat)
    
    err_lat<-qnorm(CI)*sd_lat/sqrt(length(par_ids))
    err_lat
    
    l_lat<-avg_lat-err_lat
    r_lat<-avg_lat+err_lat
    
    ##
    avg_lon<-mean(tot_dist$dist_lon)
    sd_lon<-sd(tot_dist$dist_lon)
    
    err_lon<-qnorm(CI)*sd_lon/sqrt(length(par_ids))
    err_lon
    
    l_lon<-avg_lon-err_lon
    r_lon<-avg_lon + err_lon
    
    
    ####################################################################
    ### find out what is the max del.t to determine the time at which the particle left the model at
    data_use2<-subset(data_use, del.t<=endday)
    dim(data_use2)
    
    # all particles in model simulation
    par_ID<-unique(data_use2$ID)
    tmp<-NA*1:length(par_ID)
    
    # only surface released particles
    surf_par<-subset(data_use2, ID %in% surfIDs)
    dim(surf_par) # 179386     16
    par_ID2<-unique(surf_par$ID)
    tmp2<-NA*1:length(par_ID2)
    
    # only middle released particles
    mid_par<-subset(data_use2, ID %in% midIDs)
    dim(mid_par) # 177727     16
    par_ID3<-unique(mid_par$ID)
    tmp3<-NA*1:length(par_ID3)
    
    # only deep released particles
    deep_par<-subset(data_use2, ID %in% deepIDs)
    dim(deep_par) # 180791     16
    par_ID4<-unique(deep_par$ID)
    tmp4<-NA*1:length(par_ID4)
    
    ## for loop to find which is max del.t for each particle and then put into 
    ## empty vector for each group
    
    # all particles
    for(i in 1:length(par_ID)){
    	print(i)
    	tmp[i]<-max(data_use2$del.t[which(data_use2$ID %in% par_ID[i])])
    }
    
    # surface particles
    for(i in 1:length(par_ID2)){
    	print(i)
    	tmp2[i]<-max(surf_par$del.t[which(surf_par$ID %in% par_ID2[i])])
    }
    
    # middle particles
    for(i in 1:length(par_ID3)){
    	print(i)
    	tmp3[i]<-max(mid_par$del.t[which(mid_par$ID %in% par_ID3[i])])
    }
    
    # deep particles
    for(i in 1:length(par_ID4)){
    	print(i)
    	tmp4[i]<-max(deep_par$del.t[which(deep_par$ID %in% par_ID4[i])])
    }
    
    # convert del.t into days
    tmp_tot<-tmp/96
    # max(tmp_tot)
    tmp2<-tmp2/96
    tmp3<-tmp3/96
    tmp4<-tmp4/96
    
    # add NA at the end if lengths of tmp's are different
    tmp2<-c(tmp2,NA)
    # tmp3<-c(tmp3,NA)
    tmp4<-c(tmp4,NA)
    

    ####################################################################
    ### save and return lists
    
    distretAnalysis_list<-list(init_par, surf_par, mid_par, deep_par, surf_end, 
                               mid_end, deep_end, fo_mat, avg_lat, avg_lon, 
                               sd_lat, sd_lon, err_lat, err_lon, tmp_tot, tmp2,
                               tmp3, tmp4)
    
    return(distretAnalysis_list)
    
    # Pau 
}
    


