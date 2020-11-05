### The purpose of this script is to be a one stop shop for all data post processing, analysis, and creating figures for the analyses
### Jennifer Wong-Ala
### 20201105

ppAnalysis<-function(MR_directory, transform_MR, bathy_file, rawfile_directory, file_name, lat1, lat2, lon1, lon2){
  
    ### load in MR polygons for mapping
    MR<-readOGR(dsn=path.expand(MR_directory))
    # class(MR) # SpatialPolygonsDataFrame
    # crs(MR)
    # extent(MR)
    crsmerc=CRS(transform_MR) # this transforms the model 
    MR_transformed<-spTransform(MR, CRS=crsmerc)

    ### Bathymetry
    bathy.dat<-read.table(bathy_file, sep='')
    names(bathy.dat)<-c('lon', "lat", 'depth')

    bathy.dat$depth[bathy.dat$depth>0]<-NA # Avoid points above water
    # head(bathy.dat)
    bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]


    ####################################################################
    ### initial variables or constants
    
    ## Lat and lon boundaries

    lat1<-lat1 
    lat2<-lat2

    lon1<-lon1 # -126.0909
    lon2<-lon2 

    
    ####################################################################
    ### load in raw data file
    
    setwd(rawfile_directory)
    load(file_name) 
    
    ####################################################################
    ### unlist the raw output and put into a data frame where rows are particles at each time step and columns are the variables being tracked for each particles over time (e.g., lat, lon, velocity, etc.)

    # input to create mod function would only be pop.roms (I)

    out.row<-nrow(pop.roms[[1]][[1]])
    out.col<-1+ncol(pop.roms[[1]][[1]])
    pop.fish<-as.data.frame(matrix(1:(out.col*out.row),ncol=out.col,nrow=out.row)*NA)
    pop.roms2<-list()
    pop.roms2[[1]]<-pop.fish
    pop.mat<-NULL
    # pop.mat2<-NULL # TIDES 
    # pop.roms4<-pop.roms[-1] # TIDES: use this for data after 1st day 

   
    for(l in 1:(length(pop.roms))){
      tmp<-pop.roms[[l]]
      for(j in 1:length(tmp)){
        tmp[[j]]$del.t<-j+length(tmp)*(l-1)
        pop.fish<-tmp[[j]]
        pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
    }}

    
    ## TIDES
    # l=1
    # tmp<-pop.roms[[l]]
    # 
    # for(j in 1:length(tmp)){
    #   tmp[[j]]$del.t<-j+length(tmp)*(l-1)
    #   pop.fish<-tmp[[j]]
    #   pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
    # }
    # 
    # pop.roms3<-pop.roms2 # rename tide data
    # 
    # rm(pop.roms2)
    # 
    # pop.fish<-as.data.frame(matrix(1:(out.col*out.row),ncol=out.col,nrow=out.row)*NA)
    # pop.roms2<-list()
    # pop.roms2[[1]]<-pop.fish
    # 
    # 
    # # tides data 
    # for(l in 1:(length(pop.roms4))){ 
    #   tmp<-pop.roms4[[l]]
    #   for(j in 1:length(tmp)){
    #     tmp[[j]]$del.t<-j+length(tmp)*(l-1)
    #     pop.fish<-tmp[[j]]
    #     pop.roms2[[(length(tmp)*(l-1))+j]]<-tmp[[j]]
    #   }}

    
    for(m in 1:(length(pop.roms3))){
      print(m)
      nw<-pop.roms3[[m]]
      nw$timestep<-m
      pop.mat2<-rbind(pop.mat2, nw)
    }
    
    for(m in 1:(length(pop.roms2))){
      print(m)
      nw<-pop.roms2[[m]]
      nw$timestep<-m
      pop.mat<-rbind(pop.mat, nw)
    }
    
    
    # pop.mat<-rbind(pop.mat, pop.mat2) # TIDES: combine day 1 ROMS and rest of days of data
    
    # dim(pop.mat)
    # head(pop.mat)
    # tail(pop.mat)
    

    ####################################################################
    ### clean raw data set (i.e., remove initial positions and portion of a particles trajectory that left the model domain)
    
    ## rename data to then clean
    pop.mat1<-pop.mat 
    # dim(pop.mat1) 
    
    ## unique ID numbers that have lon values less than lon1 
    lon.uni<-sort(unique(pop.mat1$ID[which(pop.mat1$lon<lon1)]))
    # length(lon.uni) 
    
    lat.uni2<-sort(unique(pop.mat1$ID[which(pop.mat1$lat>=lat2)])) 
    # length(lat.uni2) # 68
    lat.uni1<-sort(unique(pop.mat1$ID[which(pop.mat1$lat<=lat1)])) 
    # length(lat.uni1) # 700
    
    lat.uni<-sort(c(lat.uni1,lat.uni2)) # levels(lat.uni1) # only included lat.uni beceause both include all of the particles released, particles are the same for both 
    length(lat.uni) # 23958
    
    #### now I need to subset out each df for all of the individual ID numbers from uni_ids and remove the trajectory once lat or lon is equal to a boundary condition value
    uni_ids<-sort(c(lat.uni, lon.uni))
    # length(uni_ids) # 24,024
    uni_ids<-as.numeric(uni_ids)
    
    ### new dataframe without particle IDs that need to be edited to remove trajectory that touched domain border
    
    pop.mat5<-pop.mat1[!pop.mat1$ID %in% uni_ids,] 
    # dim(pop.mat5)
    
    pop.mat1$ID<-as.numeric(pop.mat1$ID)
    
    par_fill<-NULL
    
    for(i in 1:length(uni_ids)){
      
      print(i)
      
      par<-subset(pop.mat1, ID==uni_ids[i]) # subset out particles with specific particle ID
      # dim(par)
      
      num_remove<-which(par$lat<=lat1 | par$lat>=lat2 | par$lon<lon1)[1] # need to find the row number where the particle touched the boundary
      
      par_new<-par[-(num_remove:nrow(par)),] # go to df and remove the part of the trajectory that touched the boundary
      # dim(par_new)
      
      # rbind new dataset to par_fill	
      par_fill<-rbind(par_fill, par_new)	
      
    }
    
    data_use<-rbind(pop.mat5, par_fill) # final post processed data set to send back to main directory
    # dim(data_use)

    ppAnalysis_list<-list(MR_transformed, bathy_dat, data_use) # save data as a list and put into return() to return back to main working directory

return(ppAnalysis_list) # what do I want to return to list?
    
    } # end of function

