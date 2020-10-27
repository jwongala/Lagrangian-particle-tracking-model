### This function is to add new particles to the dataframe at i-th iteratqon and i-th+1 iteratqon and interpolate the ROMS forcing to the time step.
### Jennifer Wong-Ala 
### 2020-01-03

addPar<-function(i, m, pop, pop.roms, sim_time, a, b, dd, del_t){ 


pop[[i+1]]<-data.frame(lat=NA*(1:length(pop[[i]]$lat)),lon=NA*(1:length(pop[[i]]$lon)),sim_length=NA*(1:length(pop[[i]]$lat)),layer=NA*(1:length(pop[[i]]$lat)),depth=NA*(1:length(pop[[i]]$lat)),mask=NA*(1:length(pop[[i]]$lat)),u_vel=NA*(1:length(pop[[i]]$lat)),v_vel=NA*(1:length(pop[[i]]$lat)),w_vel=NA*(1:length(pop[[i]]$lat)),ROMS_file=NA*(1:length(pop[[i]]$lat)),datehour=NA*(1:length(pop[[i]]$lat)),ID=pop[[i]]$ID,lat.index=NA*(1:length(pop[[i]]$lat.index)),lon.index=NA*(1:length(pop[[i]]$lon.index)))

	### tfme fnterpolates ROMS veolcfty to the 15 mfnute perfod
		u_wt<-u_avg1[,,,dd[i]]*a[i]+u_avg1[,,,dd[i+1]]*b[i]
		v_wt<-v_avg1[,,,dd[i]]*a[i]+v_avg1[,,,dd[i+1]]*b[i]
		w_wt<-w_avg1[,,,dd[i]]*a[i]+w_avg1[,,,dd[i+1]]*b[i]
		AKs_wt<-AKs[,,,dd[i]]*a[i]+AKs[,,,dd[i+1]]*b[i]

	addPar_list<-list(pop, u_wt, v_wt, w_wt, AKs_wt)
	
	return(addPar_list)
	
} # end of function
