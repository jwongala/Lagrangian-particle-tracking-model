### The purpose of this script is to do the secondary calculations such as:

# - create connectivity matrices for the LPT model output 

### Jennifer Wong-Ala
### 20201119

ConnectivityMatrices<-function(data_use,delt_days2){

  
    ##########################################################
    ### post process data so there are not duplicates and it is in order based on ID2 column
    
    ## first subset only particles released at first time step and then add new column to init_par and end_par to show which post-processed grid cell the particle began and ended in
    
    data_use<-hi18 # assign dataset that I want to use to here
    
    init_par<-subset(data_use, del.t==1)
    dim(init_par) # 1293   17
    
    pix_id_init<-rep(NA, nrow(init_par))
    init_par2<-cbind(init_par, pix_id_init)
    
    end_par<-subset(data_use, del.t==delt_days2) # doing for two weeks (1344)
    dim(end_par) # 82 17
    
    pix_id_end<-rep(NA,nrow(end_par)) # empty column to put end grid cell for each particle 
    end_par2<-cbind(end_par, pix_id_end)
    
    ## first remove duplicated rows that have the same ID2 number
    
    tmp_init<-init_par2[!duplicated(init_par2[, c("ID2")]),]
    dim(tmp_init) # 2891  18
    head(tmp_init)
    
    tmp_end<-end_par2[!duplicated(end_par2[, c("ID2")]),] 
    dim(tmp_end) # 957  18
    head(tmp_end) 
    
    ## combine init_par2 and end_par2 into one df that will be used to make connectivity matrix 
    
    init_par3<-tmp_init[order(tmp_init$ID2), c(1:18)]
    dim(init_par3) 
    # head(init_par3)
    # tail(init_par3)
    
    end_par3<-tmp_end[order(tmp_end$ID2), c(1:18)]
    dim(end_par3)
    # head(end_par3)
    # tail(end_par3)
    
    ##########################################################
    ### find number of particles initialized at each of the 35 locations 
    
    lat_seq<-seq(45.24775, 42.00000, length.out=36) # changed it from 35 to 20 (it needs to be one more than how many grids because it is making grid cells for inbetween the lat and lon locations I give it)
    
    ## create empty list of grid cells for particle initial and end post-processed grid cell location
    lat_list_start<-as.list(1:(length(lat_seq)-1))
    lat_list_end<-as.list(1:(length(lat_seq)-1))
    
    
    ## determine what post-processed grid cell a particle began in and ended in
    for(i in 1:(length(lat_seq)-1)){
    	
    	lat_list_start[[i]]<-init_par3$ID2[init_par3 $lat<=lat_seq[i] & init_par3 $lat>lat_seq[i+1]]
    	lat_list_end[[i]]<-end_par3 $ID2[end_par3 $lat<=lat_seq[i] & end_par3 $lat>lat_seq[i+1]]
    	
    }		
    
    
    ## count how many particles begin in each of the 35 grid cells
    
    a<-NULL
    
    for (i in 1:length(lat_list_start)){
    
    	print(paste('i=', i))
    
    	a[i]<-length(lat_list_start[[i]])
    	
    }
    
    ## replace NAs in pixel_id_init column to post-processed grid cell that the particle was found in. 
     
    for(i in 1:nrow(init_par3)){
    	print(paste('i=',i))
    	
    	pp_gridcell<-NULL # empty vector to save values in; where 1 is, is the grid cell the particle is located. 
    	
    	ID<-init_par3 $ID2[i] # what is the ID number of the particle? 
    	
    	
    	for(j in 1:length(lat_list_start)){
    		
    		print(paste('j=',j))
    		
    		pp_gridcell[j]<-sum(which(lat_list_start[[j]] %in% ID))
    		
    		index<-which(pp_gridcell>=1) # what grid cell did the particle originate in? 
    		
    	}
    	
    	init_par3 $pix_id_init[i]<-index # assign the index to the pix_id_init column for the inidividual particle
    	
    }
    
    # head(init_par2)
    # tail(init_par2)
    
    
    ## replace NAs in pixel_id_end column to post-processed grid cell that the particle was found in. 
    
    for(i in 1:nrow(end_par3)){
    	
    	print(paste('i=', i))
    	
    	ppgridcell<-NULL
    	
    	ID<-end_par3$ID2[i]
    
    	for(j in 1:length(lat_list_end)){
    		
    		print(paste('j=', j))
    		
    		pp_gridcell[j]<-sum(which(lat_list_end[[j]] %in% ID))
    		
    		index<-which(pp_gridcell>=1)
    		
    	}	
    	
    	end_par3 $pix_id_end[i]<-index
    	
    }
    
    head(end_par3)
    tail(end_par3)
    
    
    ##########################################################
    ### combine init and end df to get one df 
    
    tmp_end_id<-end_par3$ID2 # create vector of particle IDs from end df
    # length(unique(tmp_end_id)) # same number of rows as end_par3! :)
    
    
    init_par4<-subset(init_par3, ID2 %in% tmp_end_id) # subset out only the particles that are found in the end df from the init df
    dim(init_par4)
    # head(init_par4)
    # tail(init_par4)
    
    pix_id_init2<-init_par4[,c(18)] # create vector of only the init grid cell locations from the init df
    length(pix_id_init2)
    
    end_par4<-cbind(end_par3, pix_id_init2) # use cbind to combine the init gridcell locations to the end_par df 
    dim(end_par4)
    head(end_par4)
    
    
    ##########################################################
    ### create connectivity matrix 
    
    ## loop through list to see if count the number of particles that begin in grid cell N_i and the number of particles that end in grid cell N_j
    
    conn_mat<-matrix(NA, length(lat_list_start), length(lat_list_end))
    dim(conn_mat) # 35  35
    
    
    for(i in 1:nrow(conn_mat)){
    	
    	print(paste('i=', i))
    	
    	# begin_par<-length(which(end_par4$pix_id_init2==i)) # how many particles are initialized with at the i-th grid cell
    	
    	# begin_seq<-rep(begin_par, length.out=nrow(conn_mat)) # create a sequence of the number of particles initialized at that grid cell 
    	
    	end_par<-rep(NA, length.out=nrow(conn_mat)) # create an empty sequence of NAs
    	
    	for(j in 1:ncol(conn_mat)){
    		
    		print(paste('j=', j)) 
    		
    		end_par[j]<-length(which(end_par4$pix_id_end==j & end_par4$pix_id_init2==i)) 
    		
    	}
    	
    	conn_vals<-end_par/a[i]
    	length(conn_vals)	
    	
    	conn_mat[i,]<-conn_vals
    	
    }
    
    
    ##########################################################
    ### test to see if each row is less than or eqaul to 1 
    
    test<-apply(conn_mat, 1, sum)
    length(test)
    head(test)
    range(test)
    
    
    ##########################################################
    ## list to return to mainLPTAnalysis.R script
  
  ConnectivityMatrices_list<-list(conn_mat, test)

return(ConnectivityMatrices_list)

}

