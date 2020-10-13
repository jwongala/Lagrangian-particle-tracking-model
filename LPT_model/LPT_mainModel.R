### Jennnifer A.T.K. Wong-Ala
## MS Thesis
## Particle Tracking Model Code for Oregon Coast (including Marine Reserves, Yaquina Bay and Estuary)
## Last worked on: 2020-01-03 (year-month-day format)

####################################################################
##### NOTES: Read before continuing to use script!

# 1. Sections noted with (***) contain placeholders or are notes that Jenn needs to attend to

####################################################################
### Prep the workspace
rm(list=ls()) # clear the workspace
options(warn=0) # warn= 0 will let me know there is an error, but won't stop the code from running; use warn=2 to stop code from running when error occurs


####################################################################
### Load libraires
# library(marmap) # to plot bathymetry
library(ncdf4) # to load netcdf advection files
library(data.table) # to use fread
library(mgcv) # to do GAMs
library(sp)
library(raster)
# library(ncdf.tools)
# library(rgeos)
# library(maptools)
# library(grDevices)
library(fields)
library(pracma)
# library(OceanView)
# library(animation)
library(R.matlab)
# library(plot3D)
# library(lattice)
# library(plotly)
library(Deriv)
library(foreach)
library(doParallel)
library(iterators)
# library(itertools)
# library(mapplots)
library(sgeostat)


####################################################################
### Load functions
setwd("/Users/jennifer/Desktop/Modularized Model") 
	# /Users/Jenn/Desktop/Model
	# /Users/jennifer/Desktop/Modularized Model
	# /Users/wongalaj/Desktop/Model/Modularized Model
	# /home/jwongala/OR_coast_AK/Modularized Model
	# /data2/modularized_model
	
source("distance.function.R") # LC function to calculate distances between two points
source("makeGrid.R")
source("makePop.R") # make initial population
source("loadNCDF.R") # load netcdf files (ROMS)
source("addPar.R") # add particles to population
source("movePar.R") # move particles with ROMS forcing


####################################################################
### Load files for lat and lon var

## Grid variables

latd<-86 # low res= 86; high res= 371
lond<-355 # low res= 355; high res= 1446
depthd<-40 # both 40
depth_exd<-41
lat_exd<-87 # low res= 87; high res= 372
lon_exd<-356 # low res= 356; high res= 1447
grid_name<-"OR_coast_AK_subdomain_grid.nc" # low res= 'OR_coast_AK_subdomain_grid.nc'; high res= MarRes_grid_3s.nc
mat_name<-"OR_coast_AK_grid.mat" # low res= "OR_coast_AK_grid.mat"; high res= "MarRes_llzh_3.mat"
grid_dir<-setwd("/Users/jennifer/Desktop/OR_coast_AK_grid")
	# setwd("/Users/jennifer/Desktop/OR_coast_AK_grid")
	# setwd("/Users/wongalaj/Desktop/Model/OR_coast_AK_grid") # directory of grid files
	# setwd("/home/jwongala/OR_coast_AK/OR_coast_AK_grid") # directory for grid file
	# setwd("/Users/Jenn/Desktop/Model/grid")
	# setwd("/data2/OR_coast_AK_grid") 

makeGrid_list<-makeGrid(grid_dir=grid_dir,
						latd=latd,
						lond=lond,
						depthd=depthd,
						lat_exd=lat_exd,
						lon_exd=lon_exd,
						grid_name=grid_name,
						mat_name=mat_name)

	## unlist makeGrid_list 
	lon_rho<-makeGrid_list[[1]] # lon coordinates at rho points
	lat_rho<-makeGrid_list[[2]] # lat coordinates at rho points
	min.lon<-makeGrid_list[[3]] # minimum lon coordinate
	min.lat<-makeGrid_list[[4]] # minimum lat values
	max.lat<-makeGrid_list[[5]]
	z_rho<-makeGrid_list[[6]]
	mask.u<-makeGrid_list[[7]]
	mask.v<-makeGrid_list[[8]]
	mask.rho<-makeGrid_list[[9]]
	boundary<-makeGrid_list[[10]]
	

####################################################################
#### INITIAL CONDITIONS AND CONSTANTS

# Radius of the Earth  
Rad<-6378137 # measured in meters, Earth’s radius, sphere, to convert m to decimal degrees via Haversine method

lon_release<-rep(-123.9263, 600) # rep(-123.9263, length.out= 100)  # (-123.9263)  #100 where to release particles (longitude coordinates)
lat_release2<-rep(44, 300)
lat_release1<-rep(45,300)
lat_release<-c(lat_release1,lat_release2)
# lat_release<-c(44,45) # seq(from= 42, to= 45.24775, length.out= 100)  # 45; where to release particles (latitude coordinates)
layer_release<-c(40,20,1) # seq(40,1,-4) # c(40, 30, 20, 10, 1) # (40) #  #  # depth layer to release  

min.depth<-max(z_rho) # minimum depth (use max because it is negative)

del_t<-15 # (units= mins; delta t for each particle that it is interpolated to) 

# total_mod_run<-0 # 44640 # (units= minutes; 31 days= this is so I can create a sequence of release par)
# hour6<-0 # 360 # 6 hours = 360 mins
# stop_par_release<-0 # 1440 # 1.21e+6 # stop releasing particles after 2 weeks; units=minutes
# release_par<-seq(from=hour6, to=stop_par_release, by= hour6 )  # sequence fed to model determining how ofter particles are released from initial conditions 

roms_step<-7200 # time_dat[2]-time_dat[1] # 7200 secs (2 hour time step from ROMS)
t_d<-8
sim_step<-roms_step/t_d # in seconds (changed it froms 7200 which was 2 hours in seconds to 15 mins; 900 in seconds)
sim_time<-(12*t_d) # 96 number of iterations 

files<-list.files("/Volumes/TOSHIBA EXT/ROMS_HiRes/2016", pattern= '*.nc', full.names=T)
	# /Volumes/TOSHIBA EXT/ROMS_HiRes/2016
	# /Volumes/TOSHIBA EXT/April2016_LowRes 
	# /Volumes/TOSHIBA EXT/OR_coast_AK_ROMS 
	# /Users/wongalaj/Desktop/April2016_LowRes 
	# /home/jwongala/OR_coast_AK/April2016_LowRes 
	# /Users/Jenn/Desktop/May2016_ROMSlowres
	# /Users/jennifer/Desktop/April2016_LowRes 
	# /data2/ROMS_output/ORcoast_HiRes/2016
	
nc.file<-NULL
dim.files<-length(files)
# files2<-floor(dim.files/12)

####################################################################
### MAKE INITIAL POPULATION

pop_list<-makePop(lon_release=lon_release,
				  lat_release=lat_release,
				  layer_release=layer_release)
		
	## unlist pop.roms and pop				  
	pop.roms<-(pop_list[[1]]) # save first list as pop.roms
	pop<-pop_list[[2]] # save second df as pop
	

####################################################################
### MODEL STARTS HERE
# seq_along(files)

system.time(

for(p in seq_along(files)) {
	
		print(paste0("p= ", p)) # print what ROMS file using
	
		### load ROMS netcdf file
	
		loadNCDF_list<-loadNCDF(p=p,
								pop.roms=pop.roms,
								pop=pop,
								files=files,
								lon_exd=lon_exd,
								lat_exd=lat_exd,
								latd=latd,
								lond=lond,
								depthd=depthd,
								depth_exd=depth_exd)
		
			## unlist loadNCDF_list
			s_rho<-loadNCDF_list[[1]] # s-coordinate at rho points (unitless) 
			srho.dims<-loadNCDF_list[[2]] # dimension of s-rho (unitless)
			v_avg1<-loadNCDF_list[[3]] # v momentum component (units= degree/s)
			u_avg1<-loadNCDF_list[[4]] # u momentum compenent (units= degree/s)
			w_avg1<-loadNCDF_list[[5]] # w momentum component (units= m/s)
			AKs<-loadNCDF_list[[6]] # salinity vertical diffusion coefficient (m2/s)
			minaks<-loadNCDF_list[[7]] # minimum AKs value
			del_tAKs<-loadNCDF_list[[8]] # time step used to in random diffusivity equation
			b<-loadNCDF_list[[9]] # coefficient needed to time interpolate ROMS forcings
			a<-loadNCDF_list[[10]] # coefficient needed to time interpolate ROMS forcings
			dd<-loadNCDF_list[[11]] # coefficient needed to time interpolate ROMS forcings
			datehour<-loadNCDF_list[[12]]
			p<-loadNCDF_list[[13]]
			pop<-loadNCDF_list[[14]]
			pop.roms<-loadNCDF_list[[15]]
	
	for(i in 1:(sim_time-1)) {
				
		print(paste0("i= ", i)) # print iteration model is on
		
			### add particles to initial population and time interpolate ROMS forcings
			addPar_list<-addPar(i=i,
							 	m=m,
							 	pop=pop,
							 	pop.roms=pop.roms,
							 	sim_time=sim_time,
							 	a=a, 
							 	b=b, 
							 	dd=dd,
							 	del_t=del_t)	
			
				## unlist addPar_list
				pop<-addPar_list[[1]]
				u_wt<-addPar_list[[2]]
				v_wt<-addPar_list[[3]]
				w_wt<-addPar_list[[4]]
				AKs_wt<-addPar_list[[5]]
			
		
		for(w in 1:nrow(pop[[i]])) {
				
		## move particles to initial population and time interpolate ROMS forcings
			movePar_list<-movePar(p=p,
								  i=i,
								  pop=pop,
								  del_t=del_t,
								  del_tAKs=del_tAKs,
								  sim_step=sim_step,
								  files=files,
								  datehour=datehour,
								  u_wt=u_wt,
								  v_wt=v_wt,
								  w_wt=w_wt,
								  AKs_wt=AKs_wt,
								  w=w)
								  
				## unlist movePar_list
				pop<-movePar_list[[1]]
		
		}  # end of w-th for loop		

		   ### save pop to pop.roms
	       pop.roms[[p]]<-pop # save each population matrix at every p-th ROMS file

	} # end of i-th for loop
		

} # end of p-th for loop


) # end of system timer



# save model run
save(pop.roms, file= "LPT_06.02.2020_JWA_init_2016_high_diffusivityTEST_5days.RData")
print("model is pau!")


### Pau :) 

