#define function to make plot output.
selectionWithSoy <- function(years = c(0:30), Nsoil0=3, pLsoil0 = 0.5, # initial conditions
                             logNnod = 5, s.nodL =2, overwinterR = 0.80, s.soilH=0.17,# paramters for population growth and selection 
                                                                                       #in and out of nodules
                             freqDependentSanctions=FALSE,a = 30, c=0.99,Smax=1.5){ # parameters for frequency-dependent sanctions
                             
  #### Nnod = of rhizobia released from nodules, RnodLoverH = # of offspring 
  #per nodule for strain L/# of offspring per nodule of srain H., winterdecline.log = about how much the rhizobia
  # population goes down between harvest and planting. survdiff = difference in # of descendants per cell at next planting season (H-L)
  ####Make vectors for output:
  pLsoil <- numeric(length(years)) # frequency of low-PHB allele
  Nsoil <- numeric(length(years)) # total population size of rhizobia in soil at planting
  Nharvest <- numeric(length(years)) # population size in the soil at harvest (for looking at oscillation)
  pLchange <- c(NA,numeric(length(years)-1)) # change in frequency of low-PHB allele
  WL <- c(NA,numeric(length(years)-1)) # fitness of low-PHB strain from year t to t+1 (measured at planting)
  WH <- c(NA,numeric(length(years)-1)) ## fitness of high-PHB strain
  s.tot = c(NA,numeric(length(years)-1)) ## cumulative selection coefficent against high-PHB strain (can be neg.)
  pLsoil[1] <- c(pLsoil0) # proportion of low-PHB strain in soil before the first host crop
  Nsoil[1] <- c(10^Nsoil0) # rhizobia population size in soil before the first host crop
  Nharvest[1] <- c(10^Nsoil0)
  Nnod = (10^logNnod)
  for(i in 2:length(years)){
    ## phase 1: selection in symbiosis:
    nodoccupancy_L = pLsoil[i-1] # nodule occupancy matches proportion in soil
    if(freqDependentSanctions==FALSE){
      nodreproductionratio = s.nodL + 1
    }else{## change selection differential in nodules based on frequency in soil (i.e., nod occupancy)
      snodLHF = Smax*c*exp(-a*(pLsoil[i-1]-0.5))/(c*exp(-a*(pLsoil[i-1]-0.5))+(1-c))
      nodreproductionratio = snodLHF + 1
    }
    # Nnod = number of rhizobia released from nodules (constant)
    # proportion of rhizobia released from nodules that are made of of strain L.
    pLnod = nodoccupancy_L/(nodreproductionratio^(-1) * (1-nodoccupancy_L) + nodoccupancy_L) 
    NLs <- pLsoil[i-1]*Nsoil[i-1] + pLnod*Nnod  ## population size of low-PHB strain after nodule senescence
    Ns <- Nsoil[i-1] + Nnod ##total population size in soil
    Nharvest[i] <- Ns
    NHs <- Ns - NLs # population size of high-PHB strain 
    
    ##phase 2: overwintering: soil population declines slightly between release from nodules and next planting season.
    ##let's say populations typically decline at about at the rate pressented by Revellin (log10 goes down by 0.18 each year.)
    # Since the time period I'm looking at is shorter, I'll say it goes down a third as much (log10 goes down by 0.06 between planting seasons)
    ## note: populations with a larger fraction of high-PHB cells will decline more slowly.
    Rwinter = overwinterR  # number of descendants per cell for low-PHB strain
    #The high-PHB strain has higher survival than the 
    NL.spring <- NLs * (Rwinter) # population size of low-PHB strain left the following winter
    NH.spring <- NHs*(Rwinter*(1+s.soilH)) # pop size of high-PHB strain left the following winter, plus fitness advantage.
    Ns.spring <- NL.spring + NH.spring
    Nsoil[i] <- Ns.spring
    pLsoil[i] <-  NL.spring/Ns.spring
    pLchange[i] <- pLsoil[i] - pLsoil[i-1]
    ##calculate net selection coefficient against the high-PHB strain (can be negative)
    WL[i] = NL.spring / (pLsoil[i-1] * Nsoil[i-1]) #fitness (descendants per parent from the previous year).
    WH[i] = NH.spring / ((1-pLsoil[i-1]) * Nsoil[i-1]) 
    s.tot[i] = (WL[i] - WH[i])/WH[i] # cumulative selection differential between growing seasons 
  }
  return(list(pL = pLsoil, Nsoil = Nsoil,Nharvest=Nharvest,pLchange = pLchange, WL=WL, WH=WH, s.tot = s.tot))
}


selectionWithSoy_stochastic <- function(years = c(0:30), Nsoil0=3, pLsoil0 = 0.5, # initial conditions
                             logNnod = 5, s.nodL =2, overwinterR = 0.80, s.soilH=0.17,# paramters for population growth and selection 
                             #in and out of nodules
                             freqDependentSanctions=FALSE,a = 30, c=0.99,Smax=1.5){ # parameters for frequency-dependent sanctions
  
  #### Nnod = of rhizobia released from nodules, RnodLoverH = # of offspring 
  #per nodule for strain L/# of offspring per nodule of srain H., winterdecline.log = about how much the rhizobia
  # population goes down between harvest and planting. survdiff = difference in # of descendants per cell at next planting season (H-L)
  ####Make vectors for output:
  pLsoil <- numeric(length(years)) # frequency of low-PHB allele
  Nsoil <- numeric(length(years)) # total population size of rhizobia in soil at planting
  Nharvest <- numeric(length(years)) # population size in the soil at harvest (for looking at oscillation)
  pLchange <- c(NA,numeric(length(years)-1)) # change in frequency of low-PHB allele
  WL <- c(NA,numeric(length(years)-1)) # fitness of low-PHB strain from year t to t+1 (measured at planting)
  WH <- c(NA,numeric(length(years)-1)) ## fitness of high-PHB strain
  s.tot = c(NA,numeric(length(years)-1)) ## cumulative selection coefficent against high-PHB strain (can be neg.)
  pLsoil[1] <- c(pLsoil0) # proportion of low-PHB strain in soil before the first host crop
  Nsoil[1] <- c(10^Nsoil0) # rhizobia population size in soil before the first host crop
  Nharvest[1] <- c(10^Nsoil0)
  Nnod = (10^logNnod)
  for(i in 2:length(years)){
    ## phase 1: selection in symbiosis:
    nodoccupancy_L = pLsoil[i-1] # nodule occupancy matches proportion in soil
    if(freqDependentSanctions==FALSE){
      nodreproductionratio = s.nodL + 1
    }else{## change selection differential in nodules based on frequency in soil (i.e., nod occupancy)
      snodLHF = Smax*c*exp(-a*(pLsoil[i-1]-0.5))/(c*exp(-a*(pLsoil[i-1]-0.5))+(1-c))
      nodreproductionratio = snodLHF + 1
    }
    # Nnod = number of rhizobia released from nodules (constant)
    # proportion of rhizobia released from nodules that are made of of strain L.
    pLnod = nodoccupancy_L/(nodreproductionratio^(-1) * (1-nodoccupancy_L) + nodoccupancy_L) 
    NLs <- pLsoil[i-1]*Nsoil[i-1] + pLnod*Nnod  ## population size of low-PHB strain after nodule senescence
    Ns <- Nsoil[i-1] + Nnod ##total population size in soil
    Nharvest[i] <- Ns
    NHs <- Ns - NLs # population size of high-PHB strain 
    
    ##phase 2: overwintering: soil population declines slightly between release from nodules and next planting season.
    ##let's say populations typically decline at about at the rate pressented by Revellin (log10 goes down by 0.18 each year.)
    # Since the time period I'm looking at is shorter, I'll say it goes down a third as much (log10 goes down by 0.06 between planting seasons)
    ## note: populations with a larger fraction of high-PHB cells will decline more slowly.
    lambda = 100* overwinterR
    ### overwinter survival (# of offspring per 100 cells is poisson distributed)
    Rwinter = rpois(n=1,lambda)/100  # number of descendants per cell for low-PHB strain
    #The high-PHB strain has higher survival than the 
    NL.spring <- NLs * (Rwinter) # population size of low-PHB strain left the following winter
    NH.spring <- NHs*(Rwinter*(1+s.soilH)) # pop size of high-PHB strain left the following winter, plus fitness advantage.
    Ns.spring <- NL.spring + NH.spring
    Nsoil[i] <- Ns.spring
    pLsoil[i] <-  NL.spring/Ns.spring
    pLchange[i] <- pLsoil[i] - pLsoil[i-1]
    ##calculate net selection coefficient against the high-PHB strain (can be negative)
    WL[i] = NL.spring / (pLsoil[i-1] * Nsoil[i-1]) #fitness (descendants per parent from the previous year).
    WH[i] = NH.spring / ((1-pLsoil[i-1]) * Nsoil[i-1]) 
    s.tot[i] = (WL[i] - WH[i])/WH[i] # cumulative selection differential between growing seasons 
  }
  return(list(pL = pLsoil, Nsoil = Nsoil,Nharvest=Nharvest,pLchange = pLchange, WL=WL, WH=WH, s.tot = s.tot))
}



# this is similar to the density-dependent sanctions setting above, except it just switches after the 
#first year (assuming differential reproduction in nodules is greater the first year due to lower population 
# density in soil.
selectionWithSoy_2phase <- function(years = c(0:30), Nsoil0=3, pLsoil0 = 0.5, # initial conditions
                             logNnod = 5, s.nodL_LD =2, s.nodL_HD = 1, overwinterR = 0.80, s.soilH = 0.17) { # constants
  #### Nnod = of rhizobia released from nodules, RnodLoverH = # of offspring 
  #per nodule for strain L/# of offspring per nodule of srain H., winterdecline.log = about how much the rhizobia
  # population goes down between harvest and planting. survdiff = difference in # of descendants per cell at next planting season (H-L)
  ####Make vectors for output:
  pLsoil <- numeric(length(years)) # frequency of low-PHB allele
  Nsoil <- numeric(length(years)) # total population size of rhizobia in soil
  pLchange <- c(NA,numeric(length(years)-1)) # change in frequency of low-PHB allele
  WL <- c(NA,numeric(length(years)-1)) # fitness of low-PHB strain from year t to t+1 (measured at planting)
  WH <- c(NA,numeric(length(years)-1)) ## fitness of high-PHB strain
  s.tot = c(NA,numeric(length(years)-1)) ## cumulative selection coefficent against high-PHB strain (can be neg.)
  pLsoil[1] <- c(pLsoil0) # proportion of low-PHB strain in soil before the first host crop
  Nsoil[1] <- c(10^Nsoil0) # rhizobia population size in soil before the first host crop
  Nnod = (10^logNnod)
  for(i in 2:length(years)){
    ## phase 1: selection in symbiosis:
    nodoccupancy_L = pLsoil[i-1] # nodule occupancy matches proportion in soil
    if(i==2){
      nodreproductionratio = s.nodL_LD + 1
    }else{## change selection differential in nodules based on soil population size
      nodreproductionratio = s.nodL_HD + 1
    }
    # Nnod = number of rhizobia released from nodules (constant)
    # proportion of rhizobia released from nodules that are made of of strain L.
    pLnod = nodoccupancy_L/(nodreproductionratio^(-1) * (1-nodoccupancy_L) + nodoccupancy_L) 
    NLs <- pLsoil[i-1]*Nsoil[i-1] + pLnod*Nnod  ## population size of low-PHB strain after nodule senescence
    Ns <- Nsoil[i-1] + Nnod ##total population size in soil 
    NHs <- Ns - NLs # population size of high-PHB strain 
    
    ##phase 2: overwintering: soil population declines slightly between release from nodules and next planting season.
    ##let's say populations typically decline at about at the rate pressented by Revellin (log10 goes down by 0.18 each year.)
    # Since the time period I'm looking at is shorter, I'll say it goes down a third as much (log10 goes down by 0.06 between planting seasons)
    ## note: populations with a larger fraction of high-PHB cells will decline more slowly.
    Rwinter = overwinterR  # number of descendants per cell for low-PHB strain
    #The high-PHB strain has higher survival than the 
    NL.spring <- NLs * (Rwinter) # population size of low-PHB strain left the following winter
    NH.spring <- NHs*(Rwinter*(1+s.soilH)) # pop size of high-PHB strain left the following winter, plus fitness advantage.
    Ns.spring <- NL.spring + NH.spring
    Nsoil[i] <- Ns.spring
    pLsoil[i] <-  NL.spring/Ns.spring
    pLchange[i] <- pLsoil[i] - pLsoil[i-1]
    ##calculate net selection coefficient against the high-PHB strain (can be negative)
    WL[i] = NL.spring / (pLsoil[i-1] * Nsoil[i-1]) #fitness (descendants per parent from the previous year).
    WH[i] = NH.spring / ((1-pLsoil[i-1]) * Nsoil[i-1]) 
    s.tot[i] = (WL[i] - WH[i])/WH[i] # cumulative selection differential between growing seasons
  }
  return(list(pL = pLsoil, Nsoil = Nsoil,pLchange = pLchange, WL=WL, WH=WH, s.tot = s.tot))
}


variablepars <- function(s.soilmax=0.25, dimsoil=100, dimnod=100,s.nodmax=2, N0=3, pL0=0.5,Nnod=5,overwinterR=0.5){
  ## scenario 1: selection remains constant.
  s.soiltest <- seq(0,s.soilmax,length.out=dimsoil) # selection differential in soil (H relative to L), bounded so strain H has zero net population growth.
  s.nodtest <- seq(0,s.nodmax,length.out=dimnod) # selection differential due to reproduction in nodules (L relative to H)
  ##make 3D array for combined selection differential (L relative to H), frequency of low-PHB allele, and change 
  # in frequency of low-PHB allele
  s.out1<- array(NA,dim=c(length(s.nodtest),length(s.soiltest),31)) 
  pL.out1 <- s.out1
  Nsoil.out1 <- s.out1 ## need the soil population size and allele frequency to figure out cumulative selection differentials over multiple years
  pLchange.out1 <- s.out1
  pL.out1[,,1] <- pL0 # starting allele frequency is always 0.5
  Nsoil.out1[,,1] <- 10^N0 # starting soil population size is always 10^3 rhizobia/g soil

## Run model for every combination of s.soil and s.nod from the vectors above.
for(i in 1:length(s.nodtest)){
  for(j in 1:length(s.soiltest)){
    m.out <- selectionWithSoy(s.soilH=s.soiltest[j],s.nodL=s.nodtest[i],pLsoil0 = pL0,
                              logNnod = Nnod,Nsoil0 = N0, overwinterR = overwinterR) # keep defaults for other parameters
    for(k in 2:31){
      Nsoil.out1[i,j,k] <- m.out$Nsoil[k]
      s.out1[i,j,k] <- m.out$s.tot[k]
      pL.out1[i,j,k] <- m.out$pL[k]
      pLchange.out1[i,j,k] <- m.out$pLchange[k]
    }}}
  pLchange1to5 <- pL.out1[,,6] - pL.out1[,,2]
  pLchange5to30 <- pL.out1[,,31] - pL.out1[,,6]
  pLchange0to1 <- pLchange.out1[,,2]
  return(list(pL.out=pL.out1,Nsoil.out=Nsoil.out1,s.out=s.out1,pLchange0to1=pLchange0to1, pLchange1to5=pLchange1to5, pLchange5to30=pLchange5to30))
}

# this adds a selection differential in the rhizosphere during the growing season (deviation from 1)
# (selection coefficients remain constant)
selectionWithSoyPlusRhizosphere <- function(years = c(0:30), Nsoil0=3, pLsoil0 = 0.5, # initial conditions
                             logNnod = 5, s.nodL =2, s.rhizH=0.1,overwinterR = 0.80, s.soilH=0.17){ # constants
  #### Nnod = of rhizobia released from nodules, RnodLoverH = # of offspring 
  #per nodule for strain L/# of offspring per nodule of srain H., winterdecline.log = about how much the rhizobia
  # population goes down between harvest and planting. survdiff = difference in # of descendants per cell at next planting season (H-L)
  ####Make vectors for output:
  pLsoil <- numeric(length(years)) # frequency of low-PHB allele
  Nsoil <- numeric(length(years)) # total population size of rhizobia in soil
  pLchange <- c(NA,numeric(length(years)-1)) # change in frequency of low-PHB allele
  WL <- c(NA,numeric(length(years)-1)) # fitness of low-PHB strain from year t to t+1 (measured at planting)
  WH <- c(NA,numeric(length(years)-1)) ## fitness of high-PHB strain
  s.tot = c(NA,numeric(length(years)-1)) ## cumulative selection coefficent against high-PHB strain (can be neg.)
  pLsoil[1] <- c(pLsoil0) # proportion of low-PHB strain in soil before the first host crop
  Nsoil[1] <- c(10^Nsoil0) # rhizobia population size in soil before the first host crop
  Nnod = (10^logNnod)
  for(i in 2:length(years)){
    ## phase 1: selection in symbiosis:
    nodoccupancy_L = pLsoil[i-1] # nodule occupancy matches proportion in soil
      nodreproductionratio = s.nodL + 1
      # Nnod = number of rhizobia released from nodules (constant)
    # proportion of rhizobia released from nodules that are made of of strain L.
    pLnod = nodoccupancy_L/(nodreproductionratio^(-1) * (1-nodoccupancy_L) + nodoccupancy_L) 
    ##meanwhile, there is a round of selection in the rhizosphere:
    NLs_planting <- pLsoil[i-1]*Nsoil[i-1] # number of low-PHB genotype present in the soil at planting.
    NHs_planting <- Nsoil[i-1]- NLs_planting # number of low-PHB genotype present in the soil at planting.
    NLs_harvest <- NLs_planting*(1- s.rhizH/2)
    NHs_harvest <- NHs_planting*(1+ s.rhizH/2)
    NLs <- NLs_harvest + pLnod * Nnod  ## population size of low-PHB strain after nodule senescence
    NHs <- NHs_harvest + (1-pLnod)*Nnod
    Ns <- NLs + NHs  ##total population size in soil 
  
    ##phase 2: overwintering: soil population declines slightly between release from nodules and next planting season.
    ##let's say populations typically decline at about at the rate pressented by Revellin (log10 goes down by 0.18 each year.)
    # Since the time period I'm looking at is shorter, I'll say it goes down a third as much (log10 goes down by 0.06 between planting seasons)
    ## note: populations with a larger fraction of high-PHB cells will decline more slowly.
    Rwinter = overwinterR  # number of descendants per cell for low-PHB strain
    #The high-PHB strain has higher survival than the 
    NL.spring <- NLs * (Rwinter) # population size of low-PHB strain left the following winter
    NH.spring <- NHs*(Rwinter*(1+s.soilH)) # pop size of high-PHB strain left the following winter, plus fitness advantage.
    Ns.spring <- NL.spring + NH.spring
    Nsoil[i] <- Ns.spring
    pLsoil[i] <-  NL.spring/Ns.spring
    pLchange[i] <- pLsoil[i] - pLsoil[i-1]
    ##calculate net selection coefficient against the high-PHB strain (can be negative)
    WL[i] = NL.spring / (pLsoil[i-1] * Nsoil[i-1]) #fitness (descendants per parent from the previous year).
    WH[i] = NH.spring / ((1-pLsoil[i-1]) * Nsoil[i-1]) 
    s.tot[i] = (WL[i] - WH[i])/WH[i] # cumulative selection differential between growing seasons 
  }
  return(list(pL = pLsoil, Nsoil = Nsoil,pLchange = pLchange, WL=WL, WH=WH, s.tot = s.tot))
}



selectionWithSoyPlusRhizosphere_carryingCapacity <- function(years = c(0:30), Nsoil0=3, pLsoil0 = 0.5, # initial conditions
                                            logNnod = 7, s.nodL =2, s.rhizH=0,overwinterR = 0.80, 
                                            s.soilH=0.17, logNcarryingCapacity=6){ # constants
  #### Nnod = of rhizobia released from nodules, RnodLoverH = # of offspring 
  #per nodule for strain L/# of offspring per nodule of srain H., winterdecline.log = about how much the rhizobia
  # population goes down between harvest and planting. survdiff = difference in # of descendants per cell at next planting season (H-L)
  ####Make vectors for output:
  pLsoil <- numeric(length(years)) # frequency of low-PHB allele
  Nsoil <- numeric(length(years)) # total population size of rhizobia in soil
  pLchange <- c(NA,numeric(length(years)-1)) # change in frequency of low-PHB allele
  WL <- c(NA,numeric(length(years)-1)) # fitness of low-PHB strain from year t to t+1 (measured at planting)
  WH <- c(NA,numeric(length(years)-1)) ## fitness of high-PHB strain
  s.tot = c(NA,numeric(length(years)-1)) ## cumulative selection coefficent against high-PHB strain (can be neg.)
  pLsoil[1] <- c(pLsoil0) # proportion of low-PHB strain in soil before the first host crop
  Nsoil[1] <- c(10^Nsoil0) # rhizobia population size in soil before the first host crop
  Nnod = (10^logNnod)
  for(i in 2:length(years)){
    ## phase 1: selection in symbiosis:
    nodoccupancy_L = pLsoil[i-1] # nodule occupancy matches proportion in soil
    nodreproductionratio = s.nodL + 1
    # Nnod = number of rhizobia released from nodules (constant)
    # proportion of rhizobia released from nodules that are made of of strain L.
    pLnod = nodoccupancy_L/(nodreproductionratio^(-1) * (1-nodoccupancy_L) + nodoccupancy_L) 
    ##meanwhile, there is a round of selection in the rhizosphere:
    NLs_planting <- pLsoil[i-1]*Nsoil[i-1] # number of low-PHB genotype present in the soil at planting.
    NHs_planting <- Nsoil[i-1]- NLs_planting # number of low-PHB genotype present in the soil at planting.
    NLs_harvest <- NLs_planting*(1- s.rhizH/2)
    NHs_harvest <- NHs_planting*(1+ s.rhizH/2)
    NLs <- NLs_harvest + pLnod * Nnod  ## population size of low-PHB strain after nodule senescence
    NHs <- NHs_harvest + (1-pLnod)*Nnod
    Ns <- NLs + NHs  ##total population size in soil 
    
    ##phase 2: overwintering: soil population declines slightly between release from nodules and next planting season.
    ##let's say populations typically decline at about at the rate pressented by Revellin (log10 goes down by 0.18 each year.)
    # Since the time period I'm looking at is shorter, I'll say it goes down a third as much (log10 goes down by 0.06 between planting seasons)
    ## note: populations with a larger fraction of high-PHB cells will decline more slowly.
    Rwinter = overwinterR  # number of descendants per cell for low-PHB strain
    #The high-PHB strain has higher survival than the 
    NL.spring <- NLs * (Rwinter) # population size of low-PHB strain left the following winter
    NH.spring <- NHs*(Rwinter*(1+s.soilH)) # pop size of high-PHB strain left the following winter, plus fitness advantage.
    Ns.spring <- NL.spring + NH.spring
    ##set an arbitrary carrying capacity--basically, keep the population size in the soil from 
    # growing after the first crop of soybean. 
    Nsoil[i] <- min(Ns.spring,10^logNcarryingCapacity)# doesn't let soil population size get above an arbitrary carrying capacity
    pLsoil[i] <-  NL.spring/Ns.spring
    pLchange[i] <- pLsoil[i] - pLsoil[i-1]
    ##calculate net selection coefficient against the high-PHB strain (can be negative)
    WL[i] = NL.spring / (pLsoil[i-1] * Nsoil[i-1]) #fitness (descendants per parent from the previous year).
    WH[i] = NH.spring / ((1-pLsoil[i-1]) * Nsoil[i-1]) 
    s.tot[i] = (WL[i] - WH[i])/WH[i] # cumulative selection differential between growing seasons 
  }
  return(list(pL = pLsoil, Nsoil = Nsoil,pLchange = pLchange, WL=WL, WH=WH, s.tot = s.tot))
}


##make a color key:
##allele frequency change:
colorkey_deltapL <- function(margins=c(5,4,4,2)+0.1,axiscex= 1.3, labcex=1.4, line=1,ylim=c(-0.3,0.3)){
  key <- matrix(seq(ylim[1],ylim[2],by=0.01),ncol=1)
  par(mar=margins)
  image(t(key),col=NA,xaxt="n",yaxt="n")
  axis(2,at=seq(1,length(key),length.out=6)/length(key), round(seq(-0.2,0.3,by=0.1),1),cex.axis=axiscex,las=1)
  
  image(t(key), zlim=c(-0.7,0),col=negpallet(50),add=T)
  image(t(key), zlim=c(0,0.3),col=pospallet(50),add=T)
}


##plot so I can see the intersection of positive and negative values at each year
negpallet <- colorRampPalette(c("blue","white"))
pospallet <- colorRampPalette(c("white","red"))
#highpallet <- colorRampPalette(c("red","yellow")) # set a range of colors for high positive values (selection diff.) 

##functions for making plots:
plotimagewithneg.deltaPL <- function(plmat,X=s.nodtest, Y=s.soiltest,plottitle="test" ,snodlim=c(0,2),ssoillim=c(0,0.5),
                                     ylim=c(-0.7,0.5),margins=c(5,4,4,2)+0.1,
                                     axiscex= 1.4, labcex=1.4,line=2.5,yaxt="s",ytext=expression("soil S"["H"]),
                                     xtext=expression("nodule S"["L"]), xat= seq(0,5,by=0.25)){
  par(mar=margins)
  image(x=X,y=Y,z=plmat,col=grey(0.5),main= plottitle,xlab="",ylab="",
        cex.lab=labcex, cex.axis=axiscex,xlim=snodlim,ylim=ssoillim,yaxt=yaxt,xaxt="s")
  #axis(1,at=seq(0,(max(ssoillim)),by=0.25),labels=xat,cex.axis=axiscex)
  mtext(side=2,text = ytext, line=line,cex=labcex)
  mtext(side=1,text = xtext, line=3.5,cex=labcex)
  image(x=X,y=Y,z=plmat,add=T,zlim=c(0,ylim[2]),col=pospallet(50))
  image(x=X,y=Y,z=plmat,add=T,zlim=c(ylim[1],0),col=negpallet(50))
  par(mar=c(5,4,4,2)+0.1)
}






