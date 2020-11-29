### The purpose of this script is to do the secondary calculations such as:

# - depth range (initial locations and for the entire simulations)
# - calculate the percent of particles advected out of the ROMS model domain
# - calculate the percent of particles kept inside of the ROMS model domain

### Jennifer Wong-Ala
### 20201108

secondaryCalcs<-function(pop.mat,surf_par,mid_par,deep_par,tmp2,tmp3,tmp4,days,
                         delt_days,tot_len,surfrange,midrange,deeprange){
    
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
    
    tenday2<-subset(surf_par, del.t==960) # why did I choose 960? 
    len2<-length(unique(tenda2y$ID)) # number of particles left at 10 days
    per_left2<-(len2/tot_len)*100 # percent of particles left at 10 days
    per_out2<-100-per_left2  # subtract percent left in model from 100 to get % of particles advected out after 10 days
    
    tenday3<-subset(mid_par, del.t==960) # why did I choose 960? 
    len3<-length(unique(tenday3$ID)) # number of particles left at 10 days
    per_left3<-(len3/tot_len)*100 # percent of particles left at 10 days
    per_out3<-100-per_left3 # subtract percent left in model from 100 to get % of particles advected out after 10 days
    
    
    tenday4<-subset(deep_par, del.t==960) # why did I choose 960? 
    len4<-length(unique(tenday$ID)) # number of particles left at 10 days
    per_left4<-(len4/tot_len)*100 # percent of particles left at 10 days
    per_out4<-100-per_left4  # subtract percent left in model from 100 to get % of particles advected out after 10 days
    
    
    ####################################################################                                                    
    #### % of particles with simulation lengths >=35 days
    
    a2<-length(which(tmp2>=23)) # 23 days 
    c2<-a/tot_len*100 
    
    a3<-length(which(tmp3>=23)) # 23 days 
    c3<-a/tot_len*100
    
    a4<-length(which(tmp4>=23)) # 23 days 
    c4<-a/tot_len*100

    ####################################################################                                                    
    
    secondaryCalcs_list<-list(sinitrange,minitrange,dinitrange,srange,mrange,
                              drange,per_out2,per_out3,perout4,a2,a3,a4,c2,c3,c4) 
    
    return(secondaryCalcs_list)
}