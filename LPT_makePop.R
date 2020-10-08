### This function is to create the initial dataframe for the model 
### Jennifer Wong-Ala 
### 20191204

makePop<-function(lon_release, lat_release, layer_release){	

	#### empty list to store results for each ROMS files inside 
	pop.roms<-list()

	#### empty list to store results for each time-step
	pop<-list() 
	pop[[1]]<-data.frame(lat=NA*(1:(length(lon_release)*length(layer_release))),lon=NA,sim_length=NA,layer=NA,depth=NA,mask=NA,u_vel=NA,v_vel=NA,w_vel=NA,ROMS_file=NA,datehour=NA,ID=(1:(length(lon_release)*length(layer_release))),lat.index=NA*(1:(length(lon_release)*length(layer_release))),lon.index=NA*(1:(length(lon_release)*length(layer_release)))) # (1:(length(lon_release)*length(layer_release)))

	# Translate from layer to depth based on location and initiaize list
	pop[[1]]$lon<-rep(lon_release,length(layer_release))
	pop[[1]]$lat<-rep(lat_release,length(layer_release))
	pop[[1]]$layer<-rep(layer_release,each=length(lat_release))

		for(k in 1:length(pop[[1]]$lon)){
			dist<-distance.function(pop[[1]]$lat[k],pop[[1]]$lon[k],lat_rho[,,1],lon_rho[,,1])
			tmp<-1*(dist==min(dist))
			lat.index<-(1:dim(lat_rho)[1])[apply(tmp,1,sum)==1]
			lon.index<-(1:dim(lat_rho)[2])[apply(tmp,2,sum)==1]
			pop[[1]]$depth[k]<-z_rho[lat.index,lon.index,pop[[1]]$layer[k]]
			pop[[1]]$mask[k]<-mask.rho[lat.index,lon.index]
			pop[[1]]$lat.index[k]<-lat.index
			pop[[1]]$lon.index[k]<-lon.index
			pop[[1]]$lon[k]<-ifelse(mask.rho[lat.index,lon.index]==1,lon_release[1],lon_rho[boundary[lon.index],1,40]) # lon_release[k] (keep lon_release[1] if releasing particles in same lat and lon locations), change back to k index when not
			}

		for(k in 1:length(pop[[1]]$lon)){
			dist<-distance.function(pop[[1]]$lat[k],pop[[1]]$lon[k],lat_rho[,,1],lon_rho[,,1])
			tmp<-1*(dist==min(dist))
			lat.index<-(1:dim(lat_rho)[1])[apply(tmp,1,sum)==1]
			lon.index<-(1:dim(lat_rho)[2])[apply(tmp,2,sum)==1]
			pop[[1]]$mask[k]<-mask.rho[lat.index,lon.index]
			pop[[1]]$lat.index[k]<-lat.index
			pop[[1]]$lon.index[k]<-lon.index
			}

	pop_list<-list(pop.roms,pop)
	
	return(pop_list)
	
} # end of function 