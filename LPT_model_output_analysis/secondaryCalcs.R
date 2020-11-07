####################################################################
#### Depth range (release location)
sinitrange<-range(pop.mat$depth[1:250])
minitrange<-range(pop.mat$depth[251:500])
dinitrange<-range(pop.mat$depth[501:750])


####################################################################
#### Depth range (entire simulation)
srange<-range(surf_par$depth, na.rm=T)
mrange<-range(mid_par$depth, na.rm=T)
drange<-range(deep_par$depth, na.rm=T)


####################################################################
#### % of particles advected out of domain after 10 days
d<-deep_par

tot_len<-250 # total number of particles subset is already contains (each layer has 250 at beginning of simulation)

tenday<-subset(d, del.t==960) # why did I choose 960? 
len<-length(unique(tenday$ID)) # number of particles left at 10 days
per_left<-(len/tot_len)*100 # percent of particles left at 10 days
per_out<-100-per_left  # subtract percent left in model from 100 to get % of particles advected out after 10 days
per_out


####################################################################                                                    
#### % of particles with simulation lengths >=35 days
tmp_use<-tmp4
# par_use<-mid_par

a<-length(which(tmp_use>=23)) # 23 days 

c<-a/b*100 


