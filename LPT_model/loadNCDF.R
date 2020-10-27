### This function loads in the ROMS files at time p and downloads the various forcings to the workspace
### Jennifer Wong-Ala
### 2020-01-03

loadNCDF<-function(p, pop.roms, pop, files, lon_exd, lat_exd, latd, lond, depthd, depth_exd){

	setwd("/Volumes/TOSHIBA EXT/ROMS_HiRes/2016") # set directory to where ROMS is
		
		# setwd(/Volumes/TOSHIBA EXT/ROMS_HiRes/2016)
		# setwd("/Users/jennifer/Desktop/April2016_LowRes")
		# setwd("/Volumes/TOSHIBA EXT/April2016_LowRes") 
		# setwd("/Users/wongalaj/Desktop/April2016_LowRes")
		# setwd("/home/jwongala/OR_coast_AK/April2016_LowRes")
		# setwd("/Users/Jenn/Desktop/May2016_ROMSlowres")
	
	if(p>1){pop[[1]]<-pop.roms[[p-1]][[sim_time]] }
	
	# files_2<-files[]
	
	### Open netcdf files
	nc.files<-nc_open(files[p])
	# print(nc.files)
	# attributes(nc.files$var)$names # look at the names of the variables for each variable 

	### Extract s_rho
	s_rho<-ncvar_get(nc.files, varid= 's_rho')
	srho.dims<-(1:length(s_rho))

	### Extract v
	v_dat<-ncvar_get(nc.files, varid= 'v') # save v variable as object
	# dim(v_dat) # 87 355  40  12 (each object is 5 dimensions that has: v, lon, lat, depth and time (hour))
	fillvalue_v<-ncatt_get(nc.files, 'v', '_FillValue') # replace the NAs with 'NA'
	v_dat[v_dat==fillvalue_v]<-NA
	v_avg<-0.5*(v_dat[2:lat_exd,,,] + v_dat[1:latd,,,]) # only surface layer
	# dim(v_avg) # 86 355  40  12
	# head(v_avg)
	v_avg1<-(v_avg/Rad)*(180/pi) # converting from m/s to degrees latitude
	v_avg1[is.na(v_avg1)]<-min(v_avg1,na.rm=T)

	### Extract u
	u_dat<-ncvar_get(nc.files, varid= 'u') # save u variable as object
	# dim(u_dat) #  86 356  40  12 (each object is 5 dimensions that has: u, lon, lat, depth and time (hour))
	fillvalue_u<-ncatt_get(nc.files, 'u', "_FillValue") # replace the NAs with 'NA'
	u_dat[u_dat==fillvalue_u]<-NA 
	u_avg<-0.5*(u_dat[,2:lon_exd,,] + u_dat[,1:lond,,]) # only surface layer
	# dim(u_avg) # 86 355  40  12
	# range(u_avg,na.rm=T)
	u_avg1<-(u_avg/(Rad*cos(pi*v_avg/180)))*(180/pi) # converting from m/s to degrees
	u_avg1[is.na(u_avg1)]<-min(u_avg1,na.rm=T)

	### Extract w
	w_dat<-ncvar_get(nc.files, varid= 'w')
	fillvalue_w<-ncatt_get(nc.files, 'w', "_FillValue") # replace the NAs with 'NA'
	w_dat[w_dat==fillvalue_w]<-NA 
	w_avg2<-0.5*(w_dat[,,2:depth_exd,] + w_dat[,,1:depthd,])
	w_avg1<-w_avg2[1:latd,1:lond,,]
	w_avg1[is.na(w_avg1)]<-min(w_avg1,na.rm=T)

	### Extract salinity vertical diffusion coefficient
	AKs2<-ncvar_get(nc.files, varid= 'AKt') # units = m2/s (AKt for high res, AKs for lowres)
	dim(AKs2)
	AKs1<-0.5*(AKs2[,,2:depth_exd,] + AKs2[,,1:depthd,])
	AKs<-AKs1[1:latd,1:lond,,]
	dim(AKs) # 86 355  40  12
	minaks<-min(AKs[which(AKs>0)]) # how to find the minimum value that doesn't equal zero
	del_tAKs<-1 #  # 1/del_tAKs changed from 200 sec to 1 based on north et al paper

	### Extract time/dates
	
	time_dat<-ncvar_get(nc.files, varid= 'ocean_time') # save time variable as an object; units= seconds since 01-01-2005 (gregorian calendar)
	
	# close nc.files
	nc_close(nc.files)
	
	datehour2<-as.POSIXct(time_dat, origin= "2005-01-01 00:00:00", tz= "GMT") # 2005 for low res and 0001 
	datehour<-as.character(datehour2)


	#### constants needed to interpolate u,v,w to run model on 15 minute timestep. (below is new calculations for time_step of 15 mins, but sim_time of 96)
	b<-rep(seq(0,1,length.out=t_d),sim_time/(t_d))
	a<-1-b
	dd<-sort(rep(1:dim(time_dat),t_d))

	loadNCDF_list<-list(s_rho, srho.dims, v_avg1, u_avg1, w_avg1, AKs, minaks, del_tAKs, b, a, dd, datehour, p, pop, pop.roms) # variables need to export from function
	
	
	return(loadNCDF_list)
	
} # end of function