#############################8/2/18
# script organizes data and runs JAGS model with nested effects of rotation, plot and plant. 
# this is for the SROC experiment only 
load("phbplusplantdata_final.RDATA")
##for now, just take rotations for effect of host absence
hostab <- phbdata[phbdata$rotation %in% c("C5","C1","S5","Cn"),
                  c("phbEst_zerobound","rotation","shootMassG","nodCount","nodMassG","plantID","plot","plate","phbslope",
                    "plotrep","plantrep")]

str(hostab)
hostab$sqrt.phb <- sqrt(hostab$phbEst_zerobound) # square-root transformation for skew
hostab$rotation <- factor(hostab$rotation,levels=c("S5","C1","C5","Cn"))
hostab$plot <- factor(hostab$plot)
hostab$yearsSinceSoy <- hostab$rotation
levels(hostab$yearsSinceSoy) <- c("0","1","5","30")


# I want plant ID to be a unique number so it's easy to put in order. 
hostab.pd <- pd[pd$rotation %in% c(c("C5","C1","S5","Cn")),]
hostab.pd$rotation <- factor(hostab.pd$rotation, levels= c("S5","C1","C5","Cn") )
hostab.pd$yearsSinceSoy <- hostab.pd$rotation
levels(hostab.pd$yearsSinceSoy) <- c("0","1","5","30")
hostab.pd$plot <- factor(hostab.pd$plot)
hostab.pd <- hostab.pd[order(hostab.pd$rotation,hostab.pd$plot),]
hostab.pd$plotnum <- rep(1:16,each=4)
hostab.pd$plantnum <- c(1:64)
hostab <- merge(hostab, hostab.pd[,c("plantID","plotnum","plantnum")])
# order the dataframe, just to help me wrap my head around stuff
hostab <- hostab[order(hostab$yearsSinceSoy,hostab$plotnum,hostab$plantnum),]

##format the data to get x,y varible names for the model.
fileNameRoot = "SROChostab-"
graphFileType="png"
seednum  = 42
# specify column names relevant to analysis.
datFrm = hostab
yName = "sqrt.phb"
x1Name = "yearsSinceSoy"
x2Name = "plotnum"
x3Name = "plantnum"

#specify desired contrasts
xContrasts =  list(
  list(c("0"),c("1"),compVal=0.0,ROPE=c(-0.1,0.1)),
  list(c("0"),c("5"),compVal=0.0,ROPE=c(-0.1,0.1)),
  list(c("0"),c("30"),compVal=0.0,ROPE=c(-0.1,0.1)),
  list(c("1"),c("5"),compVal=0.0,ROPE=c(-0.1,0.1)),
  list(c("1"),c("30"),compVal=0.0,ROPE=c(-0.1,0.1)),
  list(c("5"),c("30"),compVal=0.0,ROPE=c(-0.1,0.1))
)

## Run script with model:
source("bayesScript_rotationPlotPlant.R",echo = T) # model is commented out to avoid re-running
load("SROChostab-Mcmc.Rdata")
load("SROChostab-priors-Mcmc.Rdata")

##check priors.
priormat <- as.matrix(codaSamples_priors)
mcmcMat <- as.matrix(codaSamples3, chains=T)

plotPost(priormat[,"b0"],main="b0",xpd=F)
zz <- seq(-50,50,by=0.01)
lines(zz, dnorm(zz, mean(hostab$sqrt.phb),sd(hostab$sqrt.phb)*10))

plotPost(priormat[,"ySigma[1]"],xlim=c(0,5),main="ySigma[1]",xpd=F)
zz <- seq(0,5,by=0.01)
lines(zz, dgamma(zz, agammaShRa[1],agammaShRa[2]))

###Make diagnostic plots for convergence, shrinkage, and autocorrelation:
parameters <- varnames(codaSamples3)
postpars <- read.csv("SROChostab-SummaryInfo.csv",stringsAsFactors = F)
postpars <- postpars[1:length(parameters),1:8]

# posterior variation in SD:
postpars[grep("ySigma",postpars$X), ]
# variation from data:
ysigmaObs <- tapply(y,x1,sd)
sd(ysigmaObs)
####Commented out functions to avoid re-drawing plots.
#  for(parname in parameters[1:30]){
#    diagMCMC(codaSamples3,parName=parname, saveName = fileNameRoot,saveType = "png")
#  }
# for(parname in parameters[31:60]){
#   diagMCMC(codaSamples3,parName=parname, saveName = fileNameRoot,saveType = "png")
# }
# 
 #for(parname in parameters[61:96]){
 #  diagMCMC(codaSamples3,parName=parname, saveName = fileNameRoot,saveType = "png")
 #}
#--------------------------------------
##Diagnostics:
##Chains all look well mixed in diagnostic plots. There's some 
# shrinkage in the within-group SD.

##Autocorrelation:
postpars[postpars$ESS<10000,] # 

#----------------------------------------------------------
##### Posterior predictive check
## Get credible parameter values from the model:
b0 <- postpars[grep("b0",postpars$X),"Mode"]
b1 <- postpars[grep("b1",postpars$X),"Mode"]
b2 <- postpars[grep("b2",postpars$X),"Mode"]
b3 <- postpars[grep("b3",postpars$X),"Mode"]
sigma <- postpars[grep("ySigma\\[",postpars$X),"Mode"]
nu <- postpars[grep("nu",postpars$X),"Mode"]

##Get simulated  y estimates from t distribution:
#predicted estimates will have the same structure as the data. 
# first, get a big matrix of simulated data from posterior estimates:
set.seed(731)
nsim <- 5e4 # make the same as chain length, just for consistency. 
ysim <- matrix(NA, nrow=Ntotal, ncol=nsim)

for(n in 1:nsim){
  for(i in 1:Ntotal){
    ## model from jags: y[i] ~ dt( mu[i], 1/ySigma[x1[i]]^2, nu)
    #mu[i] <- b0 + b1[x1[i]] + b2[x2[i]] + b3[x3[i]]
    mu <- b0 + b1[x1[i]] + b2[x2[i]] + b3[x3[i]]
    s <- sigma[x1[i]]
    pred <- mu + rt(n=1, df=nu)* s # random predicted values on a t distribution. 
    #since data were truncated at 0, simulation should be too.
    ysim[i,n] <- ifelse(pred>0, pred , 0 ) 
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

x1sd.sim <- t(apply(ysim,2,function(idsim){
  return(as.numeric(tapply(idsim,x1,sd)))
}))
colnames(x1sd.sim) <- levels(hostab$yearsSinceSoy)

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
                            plantnum = 1:64,HDIlow = HDI95.x3sim[1,],HDIhigh=HDI95.x3sim[2,])
x3forplotting <- x3forplotting[order(x3forplotting$meanest.plant),]
x3forplotting$rank.plantmean <- c(1:Nx3Lvl)
hostab <- merge(hostab, x3forplotting[,c("plantnum","rank.plantmean")],by="plantnum")

# compare simulated to actual data for all factor levels (it's just 2 simulated datasets)
plot(x1 + runif(length(x1),-0.05,0.05), hostab$sqrt.phb,xlim=c(0,5))
points(x1-0.25 + runif(length(x1),-0.05,0.05), ysim[,2],xlim=c(0,5),col="red")
points(x1+0.25 + runif(length(x1),-0.05,0.05), ysim[,1],xlim=c(0,5),col="blue")

#tiff(file="SROChostabPostPredCheck1.tiff", width=1000, height=550, res=100)
# plot plants in order of mean to help visualization
xlim=c(1,Nx3Lvl)
par(mar=c(5.1,4.1,4.1,2.1) + c(0,+1,-2,0) )
plot(hostab$rank.plantmean, hostab$sqrt.phb,xlim=xlim,col=grey(0,alpha=0.5),xlab="plantID (sorted by mean PHB)",
     ylab= expression(sqrt("PHB/cell (pg)")),cex.lab=1.4,cex.axis=1.3,pch=1,lwd=2,cex=1.5,xaxt="n")
points(meanest.plant ~rank.plantmean, data=x3forplotting,pch=17,col=grey(0,alpha=0.5),cex=2)
with(x3forplotting, segments(rank.plantmean,HDIlow.y,rank.plantmean,HDIhigh.y),
     col=grey(0,alpha=0.5))
#dev.off()

## there are a few outliers, but data are generally within the 95% HDI of the model. 
## how many outliers:
numhigher <- c()
numlower <- c()
for(i in 1:64){
  numhigher <- c(numhigher,sum(hostab[hostab$rank.plantmean==i,"sqrt.phb"] > x3forplotting$HDIhigh.y[i]))
  numlower <- c(numlower,sum(hostab[hostab$rank.plantmean==i,"sqrt.phb"] < x3forplotting$HDIlow.y[i]))
}
numhigher
sum(numhigher)
numlower
sum(numlower)
numhigher+numlower
sum(numhigher,numlower) # a total of 30 nodules are outside the range of the model (more in the high range)

#tiff("SROChostabPostPredCheck2.tiff", width=600, height=550, res=100)
yrange = c(min(x3forplotting$HDIlow)-0.01, max(x3forplotting$HDIhigh)+0.01)
plot(x3forplotting$meanfromdata.plant,x3forplotting$meanest.plant ,type="n",
     ylab = "simulated plant mean",xlab="observed plant mean",
     cex.axis=1.3, cex.lab=1.4,xlim=yrange,ylim=yrange)
with(x3forplotting, segments(meanfromdata.plant,HDIlow,meanfromdata.plant,HDIhigh, col=grey(0,alpha=0.5),lwd=2))
points(x3forplotting$meanfromdata.plant,x3forplotting$meanest.plant, col=grey(0,alpha=0.5),pch=21,lwd=2,cex=2)
abline(a=0,b=1,lty=2,lwd=2)
#dev.off()

## Simulated plant means are slightly above observed plant eans in the low range (<0.4), but
# the observed means are well within the 95% HDI. 
##Conclusion from posterior preditions: model paramters do a pretty good job of simulating the data.

#-------------------------------------------------------
####Checking residuals: 
###Make sure residuals match assumptions from the model (distribution assumptions and variance assumptions)

## Level 1 residuals (individual nodules)
for(i in 1:Ntotal){
  hostab$yhat[i] <- b0 + b1[x1[i]] + b2[x2[i]] + b3[x3[i]]
}
hostab$residual <- hostab$sqrt.phb - hostab$yhat

## Level 2 residuals:
## rotation mean:
x1hat <- b0 + b1
x1resid <- x1means - x1hat

## plot mean:
x2hat <- numeric(Nx2Lvl)
for(i in 1:Nx2Lvl){
  x2hat[i] <- b0 + b1[x1Inx2[i]] + b2[i]
}
x2resid <- x2means - x2hat

## plant mean:
x3hat <- numeric(Nx3Lvl)
for( i in 1:Nx3Lvl){
  x3hat[i] <- b0 + b1[x1Inx3[i]] + b2[x2Inx3[i]] + b3[i]
}
x3resid <- x3means - x3hat

#tiff("SROChostab_QQresid.tiff", width=800, height = 680, res=100)
layout(matrix(1:4,ncol=2,byrow=T))
for(i in 1:4){
  r <- hostab[hostab$yearsSinceSoy==x1levels[i],"residual"]
  qqplot(qt(ppoints(length(r)),df=nu),y=r/sigma[i], main=x1levels[i],
         xlab="theoretical quantiles on student's t (resid/sd)",ylab="sample quantiles",cex.lab=1.4,cex.axis=1.3)
  qqline(y=r/sigma[i], distribution = function(p) qt(p, df=nu),col="red")
}
layout(1)
#dev.off()

## Since nu is pretty high, it's close to a normal distribution.
#tiff("SROChostab_QQnorm.tiff", width=800, height = 680, res=100)
layout(matrix(1:4,ncol=2,byrow=T))
for(i in 1:4){
  r <- hostab[hostab$yearsSinceSoy==x1levels[i],"residual"]
  qqplot(qnorm(ppoints(length(r)),mean=0, sd=sigma[i]),y=r, main=x1levels[i],
         xlab="theoretical quantiles (normal)",ylab="sample quantiles",cex.lab=1.4,cex.axis=1.3)
  qqline(y=r, distribution = function(p) qnorm(p, mean=0, sd=sigma[i]),col="red")
}
layout(1)
#dev.off()

#tiff("SROChostab_histresid.tiff", width=800, height = 680, res=100)
layout(matrix(1:4,ncol=2,byrow=T))
tcomb <- qt(seq(0.001, 0.999, length.out=1000),df=nu)
for(i in 1:4){
  r <- hostab[hostab$yearsSinceSoy==x1levels[i],"residual"]
  dts <- dt(tcomb, df=nu)
  hh <- hist(r/sigma[i], plot=F,breaks=10)# scales to student's t with mean = 0, sd=1
  yspan <- c(0, max(dts,hh$density)+0.01)
  plot(hh,freq=F,main=x1levels[i],xlab="residuals/sd",ylim=yspan)
  lines(tcomb,dts, col="red")
}
layout(1)
#dev.off()

### The QQ plots and histograms look OK. Some groups have a little dispersion in the high and low
# range, but it's not terrible. 

##any trends in residuals vs. fitted values?
#tiff(file="SROChostab-residualvsfitted.tiff", width=550, height=500, res=100)
layout(matrix(1:4,ncol=2,byrow=T))
par(mar=c(4,4,3,2)+0.1)
plot(hostab$yhat, hostab$residual, xlab="model estimate",ylab="level 1 resid.",
     cex.axis=1.3, cex.lab=1.4, main="nodule level")
abline(h=0,lty=2)
plot(x3hat, x3resid, xlab="model estimate",ylab="level 2 resid.",main="plant level",
     cex.axis=1.3, cex.lab=1.4)
abline(h=0,lty=2)
par(mar=c(5,4,3,2)+0.1)
plot(x2hat, x2resid, xlab="model estimate",ylab="level 2 resid.",main="plot level",
     cex.axis=1.3, cex.lab=1.4,lwd=2)
abline(h=0,lty=2)
plot(x1hat, x1resid, xlab="model estimate",ylab="level 2 resid.",main="rotation level",
     cex.axis=1.3, cex.lab=1.4,lwd=2)
abline(h=0,lty=2)
layout(1)
#dev.off()

###There's a slight pattern in the low-range of predicted values (negative residuals are consrained). 
# That reflects the fact that measurements were constrained to zero. Aside from that, there's no 
# apparent pattern. 

##Make sure variance in residuals matches assumptions in model:
resid.sd <- tapply(hostab$residual, hostab$yearsSinceSoy, sd)
sigma
yrange <- c(min(resid.sd, sigma)-0.01, max(resid.sd,sigma)+0.01)
#tiff(file="SROChostab-residvar.tiff", width=450,height=400,res=100)
par(mar=c(5,4,2,2)+0.1)
plot(resid.sd, sigma, ylab="sd estimate from model",xlab="residual sd", 
     ylim=yrange, xlim=yrange,cex=1.5,lwd=2,cex.axis=1.3, cex.lab=1.4)
abline(a=0,b=1,lty=2,lwd=2)
#dev.off()
##residual estimate isn't exactly the same as SD estimates from model, but they are reasonably close

##Look at the HDI of the difference of means with simulated data:
x1contrasts = list(
  list(c("0"),c("1")),
  list(c("0"),c("5")),
  list(c("0"),c("30")),
  list(c("1"),c("5")),
  list(c("1"),c("30")),
  list(c("5"),c("30"))
)
x1levels <- levels(hostab$yearsSinceSoy)

##make function to calculate contrasts in a large matrix of simulate data.
getcontrasts <- function(x1means.sim,x1contrasts,mnames=colnames(x1means.sim), snames=NULL){
  contmat1 <- matrix(NA, nrow = nrow(x1means.sim), ncol=length(x1contrasts))# sqrt-transformed
  contmat2 <- matrix(NA, nrow = nrow(x1means.sim), ncol=length(x1contrasts))# sqrt transformed
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
      # undo sqrt-transformation on SD
    }else{ contmat2[,i] <- rep(NA,nrow(contmat2))}
  }
  colnames(contmat1) <- cnames
  colnames(contmat2) <- cnames
  return(list(meandiff = contmat1,sdpooled = contmat2, contrasts = contmat1/contmat2))
}
x1meansd <- cbind(x1means.sim,x1sd.sim)
colnames(x1meansd) <- c(paste0("m1[",1:4,"]"),paste0("sd[",1:4,"]"))
contmat <- getcontrasts(x1meansd, x1contrasts, mnames=colnames(x1meansd)[1:4],snames=colnames(x1meansd)[5:8] )

##  compare to summary stats from the MCMCchain:
m1names <- paste0("m1[",c(1:4),"]")
sdnames <- paste0("ySigma[",c(1:4),"]")
contmatMCMC <- getcontrasts(mcmcMat,x1contrasts,mnames=m1names,snames=sdnames)

##Figure out reasonable range of equivalancy (based on effect size from sampling from population 
# with the global mean and pooled SD from posterior model)
m <- postpars[postpars$X=="b0","Mode"]
ysigma <- postpars[postpars$X %in% c("ySigma[1]","ySigma[2]","ySigma[3]","ySigma[4]"),"Mode"]
ysigmapooled <- sqrt(sum(ysigma^2)/4)
nu <- postpars[postpars$X=="nu","Mode"]
Ntotal = nrow(hostab)
nsim=5e4
ysimnull <- matrix(NA, nrow=Ntotal, ncol=nsim)
set.seed(500)
for(n in 1:nsim){
  for(i in 1:Ntotal){
    pred <- m + rt(n=1, df=nu)* ysigmapooled # random predicted values on a t distribution. 
    #since data were truncated at 0, simulation should be too.
    ysimnull[i,n] <- ifelse(pred>0, pred , 0 ) 
  }}
## get within rotation means from null data:
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
nullcontrasts.max <- apply(nullcontrasts$contrasts,1,function(x){x[which.max(abs(x))]})
nullcontrasts.max.meandiff <- apply(nullcontrasts$meandiff,1,function(x){x[which.max(abs(x))]})

##get the HDI of null contrasts:
nullhdi <- HDIofMCMC(nullcontrasts.max,credMass=0.95)
nullhdi.meandiff <- HDIofMCMC(nullcontrasts.max.meandiff,credMass=0.95)
hist(nullcontrasts.max)
abline(v=nullhdi)
hist(nullcontrasts.max.meandiff)
rope <- round(max(abs(nullhdi)),3)* c(-1,1) # take the ROPE as the boundary of the null contrast HDI furthest from 0.
rope.meandiff <- round(max(abs(nullhdi.meandiff)),3)* c(-1,1)
##boundary for effect size is 0.334. Boundary for mean difference is 0.1

##get summary stats for contrasts:
summarizecontrasts <- function(contmat,rope= rope){
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

contsum.mcmc.mean <- summarizecontrasts(contmatMCMC$meandiff, rope = rope.meandiff) # contrasts from MCMC chain
contsum.sim.contrasts <- summarizecontrasts(contmat$contrasts, rope = rope) # contrasts from MCMC chain
contsum.mcmc.contrasts <- summarizecontrasts(contmatMCMC$contrasts,rope=rope)

# methods produce similar summary stats. 
contsum.mcmc.contrasts
contsum.sim.contrasts
contsum.sim.contrasts$propInROPE
contsum.mcmc.contrasts$propInROPE

##The simulated data and MCMC yield simular p-values for contrasts (in terms of significance at p<0.05).

srocdf <- hostab
sroccontrasts <- list(contrasts = contmatMCMC, contrastsummary = list(effectsize = contsum.mcmc.contrasts, meandiff = contsum.mcmc.mean))
sroc.sim <- list(x1means = meanofx1means.sim,x1HDI95 = HDI95.x1sim,
                 x2means = meanofx2means.sim, x2HDI95 = HDI95.x2sim,
                 x3means = meanofx3means.sim, x3HDI95 = HDI95.x3sim,
                 withinplantHDI95= rbind(ysim.HDIlow.plant,ysim.HDIhigh.plant),simdata = ysim)
## save all the objects I need for making figure 2:
#save(srocdf, sroccontrasts, sroc.sim, file="SROCobjectsForFig2.RDATA")

