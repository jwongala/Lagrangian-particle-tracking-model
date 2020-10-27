### This function is to move each particle with the ROMS forcing and then add the new locations and other information back into the dataframe
### Jennifer Wong-Ala 
### 20191205

makeGrid<-function(grid_dir, latd, lond, depthd, lat_exd, lon_exd, grid_name, mat_name){
	
	grid_file<-nc_open(grid_name) # open file with grid locations MarRes_grid_3s.nc

	lon_rho1<-ncvar_get(grid_file, varid= "lon_rho") # get lon rho data
	# dim(lon_rho1) # 87 356

	lon_rho2<-array(rep(lon_rho1, depthd),dim=c(lat_exd,lon_exd, depthd))
	# dim(lon_rho2) # 87 356  40 

	lon_rho<-lon_rho2[1:latd,1:lond,]
	# dim(lon_rho) # 86 355 40

	min.lon<-min(lon_rho) # -126.0909

	lat_rho1<-ncvar_get(grid_file, varid= "lat_rho") # get lat rho data
	# dim(lat_rho1) # 87 356

	lat_rho2<-array(rep(lat_rho1,depthd),dim=c(lat_exd, lon_exd,depthd))
	# dim(lat_rho2) # 87 356  40 

	lat_rho<-lat_rho2[1:latd,1:lond,]
	# dim(lat_rho) # 86 355 40

	max.lat<-max(lat_rho) # 46.9973
	min.lat<-min(lat_rho) # 40.65895

	#### Load Depth layer info via .mat files
	ak_grid<-readMat(mat_name)
	
	# low res
	z_rho1<-array(unlist(ak_grid$Grid[31,,]),dim=c(latd,lond,depthd))
	
	# high res
	# z_rho1<-ak_grid[[4]] 
	
	# range(z_rho1)
	# dim(z_rho1) # 86 355  40
	
	z_rho<-z_rho1[1:latd,1:lond,]
	# range(z_rho)
	# dim(z_rho)

	#### masks that will be used for coastal boundary conditions 
	# 1= water, 0= land

	mask.u<-ncvar_get(grid_file, varid= "mask_u")
	# dim(mask.u) # 86  356 

	mask.v<-ncvar_get(grid_file, varid= "mask_v")
	# dim(mask.v) # 87 355

	mask_rho2<-ncvar_get(grid_file, varid= "mask_rho")
	# dim(mask_rho2) # 87 356
	mask.rho<-mask_rho2[1:latd,1:lond]
	# range(mask.rho)
	# dim(mask.rho)

	# Identify grid points which are in the boundary 
	boundary<-(1:dim(mask.rho)[2])*NA
		for(q in 1:dim(mask.rho)[2]){
			# This finds out where zero (land) is and then moves two space to the left to create it as a boundary point so that it doesn't actually touch land
			boundary[q]<-(((1:dim(mask.rho)[1])[mask.rho[,q]==0])[1]-2)
		}	
		
		boundary[is.na(boundary)]<-dim(mask.rho)[1] # wherever there are NAs it replaces it with the first dimension of mask.rho (87)


	nc_close(grid_file) # close grid file

	makeGrid_list<-list(lon_rho, lat_rho, min.lon, min.lat, max.lat, z_rho, mask.u, mask.v, mask.rho, boundary)

return(makeGrid_list)
	
} # end of function

