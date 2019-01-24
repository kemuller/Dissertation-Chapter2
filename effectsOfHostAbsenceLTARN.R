#############################8/2/18
# script organizes data and runs JAGS model with nested effects of rotation, plot and plant. 
# this is for the SROC experiment only 

load("phbplusplantdata_final.RDATA")
##for now, just take rotations for effect of host absence
hostab <- phbdata[phbdata$rotation %in% c("S15-C16","S14-C15-S16","S12-O13-C14-C15-A16"),
                  c("phbEst_zerobound","rotation","shootMassG","plantID","plot","plate","phbslope",
                    "plotrep","plantrep")]

str(hostab)
hostab$sqrt.phb <- sqrt(hostab$phbEst_zerobound) # square-root transformation for skew
hostab$rotation <- factor(hostab$rotation,levels=c("S15-C16","S14-C15-S16","S12-O13-C14-C15-A16"))
hostab$plot <- factor(hostab$plot)
hostab$yearsSinceSoy <- hostab$rotation
levels(hostab$yearsSinceSoy) <- c("0","1","3")


# I want plant ID to be a unique number so it's easy to put in order. 
hostab.pd <- pd[pd$rotation %in% unique(hostab$rotation),]
hostab.pd$rotation <- factor(hostab.pd$rotation, levels= levels(hostab$rotation) )
hostab.pd$yearsSinceSoy <- hostab.pd$rotation
levels(hostab.pd$yearsSinceSoy) <- c("0","1","3")
hostab.pd$plot <- factor(hostab.pd$plot)
hostab.pd <- hostab.pd[order(hostab.pd$rotation,hostab.pd$plot),]
hostab.pd$plotnum <- rep(1:length(unique(hostab.pd$plot)),each=4)
hostab.pd$plantnum <- c(1:nrow(hostab.pd))
hostab <- merge(hostab, hostab.pd[,c("plantID","plotnum","plantnum")])
# order the dataframe, just to help me wrap my head around stuff
hostab <- hostab[order(hostab$yearsSinceSoy,hostab$plotnum,hostab$plantnum),]

##format the data to get x,y varible names for the model.
fileNameRoot = "LTARN-"
graphFileType="png"
seednum  = 43
# specify column names relevant to analysis.
datFrm = hostab
yName = "sqrt.phb"
x1Name = "yearsSinceSoy"
x2Name = "plotnum"
x3Name = "plantnum"

#specify desired contrasts
xContrasts =  list(
  list(c("0"),c("1"),compVal=0.0,ROPE=c(-0.1,0.1)),
  list(c("0"),c("3"),compVal=0.0,ROPE=c(-0.1,0.1)),
  list(c("1"),c("3"),compVal=0.0,ROPE=c(-0.1,0.1))
)

## Run script with model:
source("bayesScript_rotationPlotPlant.R") # model is commented out so it doesn't re-run.
load("LTARN-Mcmc.Rdata") # MCMC output

#---------------------------------------------------
###Make diagnostic plots for convergence, shrinkage, and autocorrelation:
parameters <- varnames(codaSamples3)
length(parameters)
##plot code commented out to avoid re-running
# for(parname in parameters[1:35]){
#   diagMCMC(codaSamples3,parName=parname, saveName = fileNameRoot,saveType = "png")
# }
# for(parname in parameters[36:length(parameters)]){
#   diagMCMC(codaSamples3,parName=parname, saveName = fileNameRoot,saveType = "png")
# }

## Convergence: chains look well-mixed for all paramters.
## autocorrelation:
postpars <- read.csv("LTARN-SummaryInfo.csv",stringsAsFactors = F)
rownames(postpars) <- postpars$X
postpars[postpars$ESS<10000,]
## There's some autocorelation in within-group SD (ySigma). That's probably because 
## this data set has approx. equal variance among rotation treatments (unlike SROC data)
## nu also shows autocorrelation, but estimates are high enough to be close to a normal 
# distribution anyway. 
###########################------------------

####Do posterior predictive check:
b0 <- postpars["b0","Mode"]
b1 <- postpars[paste0("b1[",1:3,"]"),"Mode"]
b2 <- postpars[paste0("b2[",1:Nx2Lvl,"]"),"Mode"]
b3 <- postpars[paste0("b3[",1:Nx3Lvl,"]"),"Mode"]
sigma <- postpars[paste0("ySigma[",1:3,"]"),"Mode"]
nu <- postpars["nu","Mode"]
set.seed(735)
nsim <- 5e4 # make the same as chain length, just for consistency. 
ysim <- matrix(NA, nrow=Ntotal, ncol=nsim)

for(n in 1:nsim){
  for(i in 1:Ntotal){
    ## model from jags: y[i] ~ dt( mu[i], 1/ySigma[x1[i]]^2, nu)
    #mu[i] <- b0 + b1[x1[i]] + b2[x2[i]] + b3[x3[i]]
    mu <- b0 + b1[x1[i]] + b2[x2[i]] + b3[x3[i]]
    s <- sigma[x1[i]]
    pred <- mu + rt(n=1, df=nu)* s # random predicted values on a t distribution. 
    ysim[i,n] <- ifelse(pred>0, pred , 0 ) #since data were truncated at 0, simulation should be too.
  }}

##Then, take means within each rotation treatment and organize into columns
ysim.means <- apply(ysim, 1, mean)
ysim.HDI <- apply(ysim, 1, function(x)HDIofMCMC(x, credMass=0.95))

ysim.HDIlow.plant <- tapply(ysim.HDI[1,],x3,mean) # HDI of nodules within plant replicates. 
ysim.HDIhigh.plant <- tapply(ysim.HDI[2,],x3,mean) # (each nodule is a replicate, so I'll take the mean)

x3means.sim <- t(apply(ysim,2,function(idsim){
  return(as.numeric(tapply(idsim,x3,mean)))
}))
colnames(x3means.sim) <- unique(x3)

x1means.sim <- t(apply(ysim,2,function(idsim){
  return(as.numeric(tapply(idsim,x1,mean)))
}))
colnames(x1means.sim) <- levels(hostab$yearsSinceSoy)

x2means.sim <- t(apply(ysim,2,function(idsim){
  return(as.numeric(tapply(idsim,x2,mean)))
}))
colnames(x2means.sim) <- unique(hostab$plotnum)


##Get mean and 95% HDI of simulated mean estimates:
meanofx1means.sim <- apply(x1means.sim,2,mean)
HDI95.x1sim <- apply(x1means.sim, 2, function(x){HDIofMCMC(x, credMass=0.95)})
meanofx2means.sim <- apply(x2means.sim, 2,mean)
HDI95.x2sim <- apply(x2means.sim, 2, function(x){HDIofMCMC(x, credMass=0.95)})
meanofx3means.sim <- apply(x3means.sim, 2,mean)
HDI95.x3sim <- apply(x3means.sim, 2, function(x){HDIofMCMC(x, credMass=0.95)})

x3forplotting <- data.frame(meanest.plant = meanofx3means.sim,HDIlow.y =ysim.HDIlow.plant,
                            HDIhigh.y = ysim.HDIhigh.plant, 
                            meanfromdata.plant = tapply(hostab$sqrt.phb,hostab$plantnum,mean),
                            plantnum = 1:48,HDIlow = HDI95.x3sim[1,],HDIhigh=HDI95.x3sim[2,])
x3forplotting <- x3forplotting[order(x3forplotting$meanest.plant),]
x3forplotting$rank.plantmean <- c(1:Nx3Lvl)
hostab <- merge(hostab, x3forplotting[,c("plantnum","rank.plantmean")],by="plantnum")

# compare simulated to actual data for all factor levels

#tiff(file="LTARNPostPredCheck1.tiff", width=1000, height=550, res=100)
# plot plants in order of mean to help visualization
stretch <- 1
xlim=c(1,Nx3Lvl)
par(mar=c(5.1,4.1,4.1,2.1) + c(0,+1,-2,0) )
plot(hostab$rank.plantmean*stretch, hostab$sqrt.phb,xlim=xlim,col=grey(0,alpha=0.5),xlab="plantID (sorted by mean PHB)",
     ylab= expression(sqrt("PHB/cell (pg)")),cex.lab=1.4,cex.axis=1.3,pch=1,lwd=2,cex=1.5,xaxt="n")
points(meanest.plant ~I(rank.plantmean*stretch), data=x3forplotting,pch=17,col=grey(0,alpha=0.5),cex=2)
with(x3forplotting, segments(I(rank.plantmean*stretch),HDIlow.y,I(rank.plantmean*stretch),HDIhigh.y),
     col=grey(0,alpha=0.5))
#dev.off()
## there are a few outliers, but data are generally within the 95% HDI of the model. 

#tiff("LTARNPostPredCheck2.tiff", width=600, height=550, res=100)
yrange = c(min(x3forplotting$HDIlow)-0.01, max(x3forplotting$HDIhigh)+0.01)
plot(x3forplotting$meanfromdata.plant,x3forplotting$meanest.plant ,type="n",
     ylab = "simulated plant mean",xlab="observed plant mean",
     cex.axis=1.3, cex.lab=1.4,xlim=yrange,ylim=yrange)
with(x3forplotting, segments(meanfromdata.plant,HDIlow,meanfromdata.plant,HDIhigh, col=grey(0,alpha=0.5),lwd=2))
points(x3forplotting$meanfromdata.plant,x3forplotting$meanest.plant, col=grey(0,alpha=0.5),pch=21,lwd=2,cex=2)
abline(a=0,b=1,lty=2,lwd=2)
#dev.off()

#-----------------------------
# compare observed deflections:
x3plusx2deflection <- x3deflection*0
for(i in 1:length(x3deflection)){
  x3plusx2deflection[i] <- x3deflection[i] + x2deflection[x2Inx3[i]]
}
plot(1:3, x1means-mean(y), ylim=c(-0.3,0.3),xlim=c(0.5,3.5),pch=17,cex=2)
points(x1Inx2 -0.25, x2deflection,col="red")
points(x1Inx3 +0.25, x3deflection,col="blue")
points(x1Inx3 , x3plusx2deflection,col="green")
plot(b2~x1Inx2)
plot(b3~x1Inx3)
plot(b1)
## there doesn't seem to be a big difference in variation among plots in this dataset. 

##get residuals and model estiamtes
## Level 1 residuals (individual nodules)
for(i in 1:Ntotal){
  hostab$yhat[i] <- b0 + b1[x1[i]] + b2[x2[i]] + b3[x3[i]]
}
hostab$residual <- hostab$sqrt.phb - hostab$yhat




#---------------------------------------
##Look at the HDI of the difference of means with simulated data:
x1contrasts = list(
  list(c("0"),c("1")),
  list(c("0"),c("3")),
  list(c("1"),c("3"))
)
x1levels <- levels(hostab$yearsSinceSoy)
getcontrasts <- function(x1means.sim,x1contrasts,mnames=colnames(x1means.sim), snames=NULL){
  contmat1 <- matrix(NA, nrow = nrow(x1means.sim), ncol=length(x1contrasts))
  contmat2 <- matrix(NA, nrow = nrow(x1means.sim), ncol=length(x1contrasts))
  cnames <- rep(NA,length(x1contrasts))
  for(i in 1:length(x1contrasts)){
    left = which(x1levels==x1contrasts[[i]][1])
    right = which(x1levels==x1contrasts[[i]][2])
    contmatmeans <- x1means.sim[,mnames]
    contmat1[,i] <- apply(contmatmeans[,c(left,right)],1,function(x) x[2] - x[1] )
    cnames[i] <- paste0(x1contrasts[[i]][1],"to",x1contrasts[[i]][2])
    if(!is.null(snames)){
      contmatsds <- x1means.sim[,snames]
      contmat2[,i] <- apply(contmatsds[,c(left,right)],1,function(x)
        sqrt((x[1]^2 + x[2]^2)/2))
    }else{ contmat2[,i] <- rep(NA,nrow(contmat2))}
    
  }
  colnames(contmat1) <- cnames
  colnames(contmat2) <- cnames
  return(list(meandiff = contmat1,sdpooled = contmat2, effectsize = contmat1/contmat2))
}
##  get contrasts from MCMCchain:
mcmcMat <- as.matrix(codaSamples3,chains=T)
m1names <- paste0("m1[",c(1:3),"]")
sdnames <- paste0("ySigma[",c(1:3),"]")
contmatMCMC <- getcontrasts(mcmcMat,x1contrasts,mnames=m1names,snames=sdnames)


##Figure out reasonable range of equivalancy (based on effect size from sampling from population 
# with the global mean and pooled SD from posterior model)
m <- postpars[postpars$X=="b0","Mode"]
ysigma <- postpars[grep("ySigma\\[",postpars$X),"Mode"]
ysigmapooled <- sqrt(sum(ysigma^2)/3)
nu <- postpars[postpars$X=="nu","Mode"]
Ntotal = nrow(hostab)
nsim=5e4
ysimnull <- matrix(NA, nrow=Ntotal, ncol=nsim)
set.seed(30)
for(n in 1:nsim){
  for(i in 1:Ntotal){
    pred <- m + rt(n=1, df=nu)* ysigmapooled # random predicted values on a t distribution. 
    #since data were truncated at 0, simulation should be too.
    ysimnull[i,n] <- ifelse(pred>0, pred , 0 ) 
  }}
## get within rotation means and sample sd from null data:
x1means.nullsim <- t(apply(ysimnull,2,function(idsim){
  return(as.numeric(tapply(idsim,x1,mean)))
}))
x1sd.nullsim <- t(apply(ysimnull,2,function(idsim){
  return(as.numeric(tapply(idsim,x1,sd)))
}))
colnames(x1means.nullsim) <- paste0("mean",x1levels)
colnames(x1sd.nullsim) <- paste0("sd",x1levels)
x1meansd.nullsim <- cbind(x1means.nullsim,x1sd.nullsim )
nullcontrasts <- getcontrasts(x1meansd.nullsim,x1contrasts,mnames=colnames(x1means.nullsim),
                              snames=colnames(x1sd.nullsim))
##set ROPE based on the highest contrast for each set of pairwise comparisons (corrects for multiple tests)
nullcontrasts.max <- apply(nullcontrasts$effectsize,1,function(x){x[which.max(abs(x))]})
nullcontrasts.max.meandiff <- apply(nullcontrasts$meandiff,1,function(x){x[which.max(abs(x))]})

##get the HDI of null contrasts:
nullhdi <- HDIofMCMC(nullcontrasts.max,credMass=0.95)
nullhdi.meandiff <- HDIofMCMC(nullcontrasts.max.meandiff,credMass=0.95)
hist(nullcontrasts.max)
abline(v=nullhdi)
hist(nullcontrasts.max.meandiff)
rope <- round(max(abs(nullhdi)),3)* c(-1,1) # take the ROPE as the boundary of the null contrast HDI furthest from 0.
rope.meandiff <- round(max(abs(nullhdi.meandiff)),3)* c(-1,1)
##boundary for effect size is 0.295. Boundary for mean difference is 0.091.

##get summary stats for contrasts:
summarizecontrasts <- function(contmat,rope=c(-0.01,0.01)){
  contrastsummary <- data.frame(contrast = colnames(contmat), mean = NA, propInROPE = NA,
                                HDIlow=NA, HDIhigh=NA,rope1 = rope[1], rope2=rope[2])
  for(i in 1:length(x1contrasts)){
    contrastsummary$propInROPE[i] <-  sum(contmat[,i]>rope[1]& contmat[,i]<rope[2]) /nrow(contmat)
    contrastsummary$mean[i] <- mean(contmat[,i])
    hdis <- HDIofMCMC(contmat[,i],credMass=0.95)
    contrastsummary$HDIlow[i] <- hdis[1]
    contrastsummary$HDIhigh[i] <- hdis[2]
  }
  return(contrastsummary)
}

contsum.mcmc <- summarizecontrasts(contmatMCMC$meandiff, rope = rope.meandiff) # contrasts from MCMC chain
contsum.mcmceff <- summarizecontrasts(contmatMCMC$effectsize,rope=rope)

ltarndf <- hostab
ltarncontrasts <- list(contrasts=contmatMCMC, contrastsummary = list(effectsize = contsum.mcmceff, meandiff = contsum.mcmc))

ltarn.sim <- list(x1means = meanofx1means.sim,x1HDI95 = HDI95.x1sim,
                 x2means = meanofx2means.sim, x2HDI95 = HDI95.x2sim,
                 x3means = meanofx3means.sim, x3HDI95 = HDI95.x3sim,
                 withinplantHDI95= rbind(ysim.HDIlow.plant,ysim.HDIhigh.plant),simdata=ysim)
## save all the objects I need for making figure 2:
#save(ltarndf, ltarncontrasts, ltarn.sim, file="LTARNobjectsForFig2.RDATA")


