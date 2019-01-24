##This script runs the Bayesian models used to compare PHB #data between experimental treatments (soil sources used for #inoculation). It is called with the source function in other #scripts. 

source("DBDA2E-utilities.R")
##convert data columns to xy variable names for model
y = as.numeric(datFrm[,yName])
x1 = as.numeric(datFrm[,x1Name])
x2 = as.numeric(datFrm[,x2Name])
x3 = as.numeric(datFrm[,x3Name])
x1levels <- levels(datFrm[,x1Name])
x2levels <- levels(factor(datFrm[,x2Name]))
x3levels <- levels(factor(datFrm[,x3Name]))

#Get number of observations and levels
Ntotal = length(y)
Nx1Lvl = length(unique(x1levels))
Nx2Lvl = length(unique(x2levels))
Nx3Lvl = length(unique(x3levels))
x1Inx2 = tapply(x1,x2,unique)
x2Inx3 = tapply(x2,x3,unique)
x1Inx3 = tapply(x1,x3,unique)
Nx3Inx2 <- tapply(x3,x2,function(x)length(unique(x)))
Nx2Inx1 <- tapply(x2,x1, function(x)length(unique(x)))
Nx3Inx1 <- tapply(x3,x1,function(x)length(unique(x)))



#calculate parameters for hyper-priors on sd:
agammaShRa <- unlist(gammaShRaFromModeSD(mode = sd(y)/2, sd = 2*sd(y)))

##for setting lower limit on within group SD:
cellSDs = tapply(y, x1 ,sd)
medianCellSD = median(cellSDs, na.rm=TRUE)

##for priors on group-level deviance (plant, plot)
# I want sd to be larger to prevent shrinkage:


# specify data in a form that's compatible with the model:
dataList = list(
  y = y, # comment out for priors
  x1 = x1,
  x2 = x2,
  x3 = x3,
  Ntotal = Ntotal,
  Nx1Lvl = Nx1Lvl,
  Nx2Lvl = Nx2Lvl,
  Nx3Lvl = Nx3Lvl,
  x1Inx2 = x1Inx2,
  x2Inx3 = x2Inx3,
  x1Inx3 = x1Inx3,
  Nx3Inx2 = Nx3Inx2,
  Nx2Inx1 = Nx2Inx1,
  Nx3Inx1 = Nx3Inx1,
  # data properties for scaling the prior:
  yMean = mean(y),
  ySD = sd(y),
  agammaShRa = agammaShRa,
  medianCellSD = medianCellSD
)

dataList_priors = list(
  #y = y, # comment out for priors
  x1 = x1,
  x2 = x2,
  x3 = x3,
  Ntotal = Ntotal,
  Nx1Lvl = Nx1Lvl,
  Nx2Lvl = Nx2Lvl,
  Nx3Lvl = Nx3Lvl,
  x1Inx2 = x1Inx2,
  x2Inx3 = x2Inx3,
  x1Inx3 = x1Inx3,
  Nx3Inx2 = Nx3Inx2,
  Nx2Inx1 = Nx2Inx1,
  Nx3Inx1 = Nx3Inx1,
  # data properties for scaling the prior:
  yMean = mean(y),
  ySD = sd(y),
  agammaShRa = agammaShRa,
  medianCellSD = medianCellSD
)


##Calculate initial values:
x1means <- as.numeric(tapply(y,x1,mean))
x1sds <- as.numeric(tapply(y,x1,sd))
x2means <- as.numeric(tapply(y,x2,mean))
x2deflection <- x2means*0
for(i in 1:Nx2Lvl){
  x2deflection[i] <- x2means[i] - x1means[x1Inx2[i]]
}
x3means <- as.numeric(tapply(y,x3,mean))
x3deflection <- x3means*0
for(i in 1:Nx3Lvl){
  x3deflection[i] <- x3means[i] - x2means[x2Inx3[i]]
}
initsList = list(
  a0 = mean(y),
  a1 = c(x1means - mean(y)),
  a2 = x2deflection,
  a3 = x3deflection
  #ySigma = mean(c(xsds)), Model doesn't work if I specify initial values for ySigma
  #a1Sigma = sd(x1means - mean(y)),
  #a2Sigma = sd(x2deflection),
  #a3Sigma = sd(x3deflection)
)

# prior for SD of deflections (broad gamma distribution)
#a1Sigma ~ dgamma(agammaShRa[1], agammaShRa[2]) # 
#a2Sigma ~ dgamma(agammaShRa[1], agammaShRa[2])
#a3Sigma ~ dgamma(agammaShRa[1], agammaShRa[2])
####Figure out what priors to use for SD of deflections.
#plotpoints <- qgamma(seq(0.025,0.975,length.out=1000),shape=agammaShRa[1],rate=agammaShRa[2])
#plot(plotpoints, dgamma(plotpoints, shape=agammaShRa[1], rate=agammaShRa[2]))
#abline(v=c(initsList$a1Sigma,initsList$a2Sigma,initsList$a3Sigma), col=c("red","blue","green"))
# the actual between-group dev. is below maximum probability on gamma dist.
# two options are to make it a large constant
#plotpoints <- qunif(seq(0.025,0.975,length.out=1000),0,(10*sd(y)))
#plot(plotpoints, dunif(plotpoints, 0,10*sd(y)))

#plotpoints <- qnorm(seq(0.025,0.975,length.out=1000),mean=0,sd=(100*sd(y)))
#plot(plotpoints, dnorm(plotpoints, mean=0,sd=100*sd(y)),xlim=c(-3,3))
#abline(v=c(x1means-mean(y)))
#abline(v=c(x2deflection), col="red")
#abline(v=c(x3deflection), col="green")

###Define the model:
modelstring = "
model {
##Likelihood (in JAGS, mu could be specified in or out of the loop.)
  for( i in 1:Ntotal){
    y[i] ~ dt( mu[i], 1/ySigma[x1[i]]^2, nu)
    mu[i] <- a0 + a1[x1[i]] + a2[x2[i]] + a3[x3[i]]
 }
  
# prior for grand mean:
a0 ~ dnorm(yMean, 1/(ySD*10)^2) # set non-committal prior (data sd *10)
## Priors for deflections (deflection due to plant includes effects of rotation and plot)
for(j1 in 1:Nx1Lvl){  # deflection due to rotation normally distributed around 0.
  a1[j1] ~ dnorm(0, 1/aSigma^2)}
for(j2 in 1:Nx2Lvl){
  a2[j2] ~ dnorm(0, 1/aSigma^2)}
for(j3 in 1:Nx3Lvl){
  a3[j3] ~ dnorm(0, 1/aSigma^2)}
# prior for SD of deflections (broad gamma distribution)
# setting these as gamma leads to high shrinkage. Try 
# setting sd of deflections as a
#a1Sigma ~ dgamma(agammaShRa[1], agammaShRa[2]) # 
#a2Sigma ~ dgamma(agammaShRa[1], agammaShRa[2])
#a3Sigma ~ dgamma(agammaShRa[1], agammaShRa[2])
aSigma = 100* ySD # set sd of deflections as large constant

## Prior for Nu
  nu <- nuMinus1 + 1 # nu must be >=1
  nuMinus1 ~ dexp(1/29) # prior on nu-1, arbitray, so it's broad.

  ##allow variance to differ among rotation treatments:
for(j in 1:Nx1Lvl){
  ySigma[j] <- max(sigma[j], medianCellSD/1000) # prevent zero ySigma
  sigma[j] ~ dgamma( ySigmaSh , ySigmaRa ) # estimate group's scale from hierarchical distribution. 
  #sigma[j] ~ dgamma(agammaShRa[1], agammaShRa[2]) # estimates each group's scale separately 
}
# priors for mode and sd of ySigma, converted to rate and shape
ySigmaSh <- 1 + ySigmaMode * ySigmaRa 
ySigmaRa <- ( ( ySigmaMode + sqrt( ySigmaMode^2 + 4 * ySigmaSD^2))/ (2*ySigmaSD^2))
ySigmaMode ~ dgamma(agammaShRa[1], agammaShRa[2])
ySigmaSD ~ dgamma(agammaShRa[1], agammaShRa[2])
#

##Convert a0, a[] to sum-to-zero b0, b[]
# get cell means (organize in 3D array)
### plant within plot:
## Get means for each nested group:
for(j3 in 1:Nx3Lvl){ # mean for plantID
  m3[j3] <- a0 + a1[x1Inx3[j3]]  + a2[x2Inx3[j3]] + a3[j3]
}
# arrange into a 2D array to get means for plot:
for(j3 in 1:Nx3Lvl){for(j2 in 1:Nx2Lvl){
 m3v2[j3,j2] <- ifelse(x2Inx3[j3]==j2 ,
                       m3[j3], 0) 
}}
for(j2 in 1:Nx2Lvl){
  m2[j2] <- sum(m3v2[1:Nx3Lvl, j2])/Nx3Inx2[j2]
}
## arrange into a 2D array to get means for rotation:
for(j2 in 1:Nx2Lvl){for(j1 in 1:Nx1Lvl){
  m2v1[j2,j1] <- ifelse(x1Inx2[j2] == j1, 
                 m2[j2] , 0)
}}
## mean estimates for rotation:
for(j1 in 1:Nx1Lvl){ 
  m1[j1] <- sum(m2v1[1:Nx2Lvl, j1])/Nx2Inx1[j1]
}
## get sum-to-zero parameters:
b0 <- mean(m3) # grand mean
for(j1 in 1:Nx1Lvl){
  b1[j1] <- m1[j1] - b0 # deflection due to rotation.
}
for(j2 in 1:Nx2Lvl){
  b2[j2] <- m2[j2] - m1[x1Inx2[j2]] # deflection due to plot (within rotation)
}
for(j3 in 1:Nx3Lvl){
  b3[j3] <- m3[j3] - m2[x2Inx3[j3]]
} 
}" # close quotes for modelstring
writeLines(modelstring, con="TEMPmodel.txt") # write model to a file


set.seed(seednum)
parameters = c("b0","b1","b2","b3","m1","ySigma","nu","ySigmaMode","ySigmaSD")
#parameters = c("b0","b1","b2","b3","m1","ySigma","nu")
require(runjags)
adaptSteps = 1000
burnInSteps = 2000
nChains = 4
thinSteps = 1
numSavedSteps = 50000
## Comment out to run the script without re-running the model
# runJagsOut <- run.jags(method="parallel",
#                        model = "TEMPmodel.txt",
#                        monitor = parameters,
#                        data = dataList,
#                        inits = initsList,
#                        n.chains = nChains,
#                        adapt = adaptSteps,
#                        burnin = burnInSteps,
#                        sample = ceiling(numSavedSteps/nChains),
#                        thin = thinSteps,
#                        summarise = FALSE,
#                        plots = FALSE)
# codaSamples3 = as.mcmc.list(runJagsOut)
# save(codaSamples3,datFrm, file = paste0(fileNameRoot,"Mcmc.Rdata")) # saved dataframe used
#to generate model alongside model output.

###For getting priors from model (comment out to avoid re-running model)
# runJagsOut_priors <- run.jags(method="parallel",
#                        model = "TEMPmodel.txt",
#                        monitor = parameters,
#                        data = dataList_priors,
#                        inits = initsList,
#                        n.chains = nChains,
#                        adapt = adaptSteps,
#                        burnin = burnInSteps,
#                        sample = ceiling(numSavedSteps/nChains),
#                        thin = thinSteps,
#                        summarise = FALSE,
#                        plots = FALSE)
#  
#   codaSamples_priors = as.mcmc.list(runJagsOut_priors)
#   save(codaSamples_priors, file = paste0(fileNameRoot,"priors-Mcmc.Rdata"))
#-------------------------------------------
# make a function to get summary info:
smryMCMC = function(codaSamples, datFrm = NULL, xName = NULL, 
                    contrasts = NULL,saveName=NULL) {
  # all single parameters:
  parameterNames = varnames(codaSamples)
  if ( !is.null(datFrm) & !is.null(xName) ) {
    xlevels = levels(as.factor(datFrm[,xName]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for(parName in parameterNames){
    summaryInfo = rbind(summaryInfo, summarizePost(mcmcMat[,parName]))
    thisRowName = parName
    if(!is.null(datFrm) & !is.null(xName)) {
      # for row name, extract numeric digits from parameter name
      levelVal = as.numeric(
        grep("^[1-9]", # grep substrings that begin with digits
             unlist(strsplit(parName , "\\[|,|\\]")),
             value=TRUE) )
      if(length(levelVal) > 0) {
        thisRowName = paste(thisRowName, xlevels[levelVal])
      }
    }
    rownames(summaryInfo)[nrow(summaryInfo)] = thisRowName
  }
  ## all contrasts:
  if( !is.null(contrasts)){
    if ( is.null(datFrm) | is.null(xName)){
      show(" *** You must specify the data file and factor names to do contrats. ***\n")
    } else{
      for( cIdx in 1:length(contrasts)){
        thisContrast = contrasts[[cIdx]]
        left = right = rep(FALSE, length(xlevels))
        for( nIdx in 1:length( thisContrast[[1]])) {
          left = left | xlevels==thisContrast[[1]][nIdx]
        }
        left = normalize(left)
        for( nIdx in 1:length(thisContrast[[2]] )) {
          right = right | xlevels == thisContrast[[2]][nIdx]
        }
      right = normalize(right)
      contrastCoef = matrix(left - right, ncol=1)
      postContrast = (mcmcMat[,paste("b1[",1:length(xlevels),"]", sep="")]
                      %*% contrastCoef )
      summaryInfo = rbind( summaryInfo, 
                           summarizePost( postContrast , 
                                          compVal = thisContrast$compVal , 
                                          ROPE = thisContrast$ROPE ))
      rownames(summaryInfo)[nrow(summaryInfo)] = (
        paste(paste(thisContrast[[1]], collapse=""),"v.",
              paste(thisContrast[[2]], collapse=""),sep="") )
      }# closes for( cIdx in...)
    }#closes else{...
  } # closes if(is.null(datFrm)...)
  # save results
  if( !is.null(saveName)) {
    write.csv(summaryInfo, file = paste(saveName, "SummaryInfo.csv",sep=""))
  }
  return( summaryInfo)
}# closes function

##save summary info (commented to avoid re-running)
#summaryInfo <- smryMCMC(codaSamples = codaSamples3, contrasts = NULL,saveName = fileNameRoot)


##----------------------------
