### This function is to move each particle with the ROMS forcing and then add the new locations and other information back into the dataframe
### Jennifer Wong-Ala 
### 2020-01-03

movePar<-function(p, i, pop, del_t, del_tAKs, sim_step, files, datehour, u_wt, v_wt, w_wt, AKs_wt, w){
	

	
		## Find and calculate lat.index and lon.index
		layer<-pop[[i]]$layer[w]
		dist<-distance.function(pop[[i]]$lat[w], pop[[i]]$lon[w],lat_rho[,,1],lon_rho[,,1])
		tmp<-1*(dist==min(dist))
		lat.index<-(1:dim(lat_rho)[1])[apply(tmp,1,sum)==1]
		lon.index<-(1:dim(lat_rho)[2])[apply(tmp,2,sum)==1]

		pop[[i]]$lat.index[w]<-lat.index
		pop[[i]]$lon.index[w]<-lon.index

		# extract currents at particle positions - the kick
		f1_u<-u_wt[lat.index,lon.index,layer]
		f1_v<-v_wt[lat.index,lon.index,layer]
		f1_w<-w_wt[lat.index,lon.index,layer]
		f1_AKs<-AKs_wt[lat.index,lon.index,layer]

		## Euler integration - the push: predicts new postion for each particle based on predicted currents at particle positions

		## Euler integration (longitude)
		D2_lon<-sim_step*f1_u  
		new.lon.tmp<-max(D2_lon+pop[[i]]$lon[w],min.lon)

		## Coastal boundary conditions
		new.lon1<-ifelse(mask.rho[lat.index,lon.index]==1,new.lon.tmp,lon_rho[boundary[lon.index],1,40])

		### Euler integration (latitude)
		D2_lat<-sim_step*f1_v  

		## Latitudinal boundary conditions
		new.lat.tmp3<-D2_lat+pop[[i]]$lat[w]
		new.lat.tmp2<-max(new.lat.tmp3,min.lat)
		new.lat.tmp<-min(new.lat.tmp2,max.lat)

		## Calculate 1st derivative to calculate AKs
		derivative<- deriv(~ x^2, "x") 
		x<-f1_AKs
		AKs_datP<-eval(derivative)
		AKs_datP<-attr(AKs_datP, "gradient")[1,]
		rnorm2<-function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
		R<-rnorm2(100,0,1)[1,]
		mean(R)  ## 4
		sd(R)    # R<-runif(seq(-1,1), min= -1, max= 1)[1]
		d<-1

		## Euler integration (depth)
		D2_w<-sim_step*f1_w
		new.z3<-D2_w+pop[[i]]$depth[w]
		
		### random displacement model
		rand_disp1<-((AKs_datP*del_tAKs) + R*(((2*(d^-1))+ (0.5*AKs_datP*del_tAKs))^0.5))
		new.z_rand<-new.z3 + rand_disp1
		
		new.z2<-min(min.depth, new.z_rand)
		deepest.depth<-z_rho[lat.index,lon.index,1]
		new.z<-max(new.z2,deepest.depth)

		## Return back to layer from depth
		delta_layer<-(z_rho[lat.index,lon.index,srho.dims[2:40]]-z_rho[lat.index,lon.index,srho.dims[1:39]])/2
		upper_layer<-c(z_rho[lat.index,lon.index,srho.dims[1:39]]+delta_layer,0)
		layer<-(srho.dims)[new.z<upper_layer][1]

		pop[[i+1]]$lat[w]<-new.lat.tmp
		pop[[i+1]]$lon[w]<-new.lon1
		pop[[i+1]]$layer[w]<-layer
		pop[[i+1]]$depth[w]<-new.z
		pop[[i+1]]$mask[w]<-mask.rho[lat.index,lon.index]
		pop[[i+1]]$sim_length[w]<-ifelse(is.na(pop[[i]]$sim_length[w])==T, del_t, (pop[[i]]$sim_length[w] + del_t)) # ifelse(p==1, p*i*del_t, pop[[i]]$sim_length[w]+del_t) # p*i*del_t # (pop[[i]]$sim_length[w] + del_t) #  # units= minutes always make this 15 when adding a new particle 
		pop[[i]]$datehour<-datehour[6] # is this the best way to do this? 
		pop[[i]]$ROMS_file<-seq_along(files)[p]
		pop[[i]]$u_vel[w]<-f1_u
		pop[[i]]$v_vel[w]<-f1_v
		pop[[i]]$w_vel[w]<-f1_w
		
	
	movePar_list<-list(pop)	
		

	
	return(movePar_list) # return list from function
	
} # end of function

