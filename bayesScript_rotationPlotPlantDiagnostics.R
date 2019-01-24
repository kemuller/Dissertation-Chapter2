#source("bayesScript_rotationPlotPlant.R")
load("rotationPlotPlant-Mcmc.Rdata") # load output from MCMC sampling (commented out in script)
load("rotationPlotPlant-priors-Mcmc.Rdata") # simulated priors. 
postpars <- read.csv("rotationPlotPlant-SummaryInfo.csv",stringsAsFactors = F)
rownames(postpars) <- postpars$X
source("allaboutthebayes/DBDA2E-utilities.R")


####organize data to help run simulation
x1 <- as.numeric(hostab$yearsSinceSoy)
x2 <- hostab$plotnum
x3 <- hostab$plantnum
yobs <- hostab$sqrt.phb
x1inx2 <- tapply(x1,x2,unique)
x1inx3 <- tapply(x1,x3,unique)
x2inx3 <- tapply(x2,x3,unique)
Nx1Lvl <- length(unique(x1))
Nx2Lvl <- length(unique(x2))
Nx3Lvl <- length(unique(x3))
Ntotal <- length(yobs)



#---------------------------------------------------------------------------------------
## 1) Check convergence and autocorrelation of JAGS sampling:
## See diagnostic plots starting with rotationPlotPlant- (made in the model script). 

# Convergence:
# Clicking through diagnostic plots shows that chains are well-mixed for all parameters. The only exception
# is a1Sigma (variance in deflections for rotation treatment), where one chain sitcks out between iterations
# 2000 and 6000. But they are well mixed by the end (I re-ran the model with a longer burn-in and adaptation period
# and the chains looked more well-mixed).

# also look at variability in chain estimates (MCSE)
mcmcMat <- as.matrix(codaSamples3)
sumstats <- summary(codaSamples3)
postpars[1:99,"SD"] <- sumstats$statistics[,"SD"]
postpars$MCSE <- NA
for(i in 1:99){
  postpars$MCSE[i] <- postpars$SD[i]/sqrt(postpars$ESS[i])
}


## Shrinkage: There's some shrinkage in higher-level paramters (a1, a2, a3 sigma, )
postpars[postpars$ESS <10000, c("ESS","Mean","Median","HDIlow","HDIhigh","SD","MCSE")]

## Diagnostics for baseline mean and deflections among rotation treatments (b1) are all good. 
# There is some autocorrelation for deflection due to plot in some plot replicates ()
# a deflection (b2) for a couple plot replicates in the 30 year- corn treatment have slightly high autocorrelation (takes about 20 steps to get an independent estimate).
# Plants with more autocorrelation in deflection estimates are mixed. But ESS isn't much lower than 10,000.
x1inx2[c(14,15)] # rotation treatments in plots with more autocorrelation in deflection from rotation mean (b3)
x1inx3[c(11,25,26,35,36,50,62)] # rotation treatments for plants with more autocorrelation in deflection from plot mean (b2)

# the higher-level parameters (a1Sigma, a2Sigma, a3Sigma, etc.) have more autocorrelation, but that's OK, since
# I'm not using them to predict anything. 


#----------------------------------------------------------------------------
## Get credible parameter values from the model:
b0 <- postpars["b0","Mode"]
b1 <- postpars[grep("b1",postpars$X),"Mode"]
b2 <- postpars[grep("b2",postpars$X),"Mode"]
b3 <- postpars[grep("b3",postpars$X),"Mode"]
sigma <- postpars[grep("ySigma\\[",postpars$X),"Mode"]
nu <- postpars["nu","Mode"]


#------------------------------------------------------------------------
##Does estimated deflection among replicates (plants, plots, etc.) 
# differ among rotation treatments?
b2plusb3 <- numeric(length(b3)) 
for(i in 1:length(b3)){
  b2plusb3[i] <- b3[i] + b2[x2inx3[i]]
}

#png(file="hostab-deflectionsOfReplicates.png")
yrange <- c(min(b2,b3)-0.05, max(b2,b3)+0.05)
set.seed(529)
jig1 <- runif(length(b2),-0.1,0.1)
jig2 <- runif(length(b3),-0.1,0.1)
plot(c(1:4),c(seq(min(b3),max(b3), length.out=4)),ylim=yrange, xlim=c(0.5,4.5),type="n",xaxt="n",
     ylab="deflection from group mean",xlab="years since host crop")
axis(1,at=c(1,2,3,4), labels=c("0","1","5","30"))
abline(h=0,lty=2)
points((x1inx2-0.3)+jig1/2, b2,col=rgb(1,0,0,alpha=0.5),pch=19)
points(x1inx3 + jig2/2, b2plusb3,col=grey(0,alpha=0.5))
points((x1inx3+0.3)+jig2/2, b3,col=rgb(0,0,1,alpha=0.5))
mtext(3,at=c(1,2,3), text=c("rotation -> plot", "rotation -> plant","plot -> plant"),col=c("red","black","blue"))
#dev.off()

# the 30-year corn treatment had greater variability among plot replicates (mainly due to one replicate with very 
# low estimates). The means for two plots were higher and the means for two plots were lower.
# Does the greater variance in the 30-year corn treatment reflect something biological, or is 
# it an artifact of calibration error?
zz <- data.frame(plantnum = c(1:64),b3 = b3,b2plusb3 = b2plusb3)
hostab <- merge(hostab, zz, by="plantnum")
zz <- data.frame(plotnum = c(1:16),b2 = b2)
hostab <- merge(hostab,zz,by="plotnum")

#png(file="hostab-deflectionVsCalibrationSlope.png",width=500,height=500)
layout(matrix(1:4,ncol=2,byrow=T))
for(i in 1:4){
plot(b2plusb3 ~ phbslope, data=hostab, type="n",
     ylab="deflection from rotation mean",xlab="calibration slope",main=levels(hostab$yearsSinceSoy)[i])
  abline(h=hostab[hostab$yearsSinceSoy==levels(hostab$yearsSinceSoy)[i],"b2"],col="red")
  points(b2plusb3 ~ phbslope, data=hostab[hostab$yearsSinceSoy==levels(hostab$yearsSinceSoy)[i],])
}

# each plant was measured on the same plate (same calibration slope). it looks like the lower-than-average
# PHB plot in the 30-year corn treatment is genuine, not an artifact of calibration error 
# (plants within that plot are spread across multiple plates). However, the two much thigher-than average plants
# may have some effect of calibration error. 
#dev.off()
layout(1)

#----------------------------------------------------------
####Get mean estimates from model to look at residuals
mhat <- data.frame(x1 = x1, x2=x2, x3=x3)
for(i in 1:nrow(mhat)){
  mhat[i , "b0"] <- b0
  mhat[i , "b1"] <- b1[x1[i]]
  mhat[i , "b2"] <- b2[x2[i]]
  mhat[i , "b3"] <- b3[x3[i]]
}
mhat$yhat <- mhat$b0 + mhat$b1 + mhat$b2 + mhat$b3
hostab$yhat <- mhat$yhat

tapply(hostab$yhat,hostab$yearsSinceSoy,mean)
tapply(hostab$sqrt.phb,hostab$yearsSinceSoy,mean)

## double check that factor levels line up:
checkx1 <- cbind(as.numeric(hostab$yearsSinceSoy),mhat$x1)
sum(apply(checkx1, 1, function(x) x[1]==x[2]))/nrow(checkx1) # yep
checkx2 <- cbind(as.numeric(hostab$plot), mhat$x2)
sum(apply(checkx2, 1, function(x) x[1]==x[2]))/nrow(checkx2) #yep
checkx3 <- cbind(as.numeric(as.factor(hostab$plantID)), mhat$x3)
sum(apply(checkx3, 1, function(x) x[1]==x[2]))/nrow(checkx3) # yep

##Calculate level-1 residuals:
hostab$residual <- hostab$sqrt.phb - hostab$yhat

plot(hostab$residual ~ hostab$yhat)
abline(h=0,col="blue")
# there's slightly wider variance in residuals at larger predicted values, but there 
# doesn't appear to be a relationship. The model didn't assume group-level heteroscedacity,
# so it might makes more sense to plot these by group
library(ggplot2)

#png("hostab-L1residualsByGroup.png",width=550, height=500)
qplot(data=hostab, x=yhat, y=residual,geom="point",alpha=0.5) + 
  facet_wrap(~hostab$yearsSinceSoy) + geom_abline(aes(slope=0, intercept = 0)) + theme_bw()
#dev.off()

##Looks OK.

plot(residual ~ yearsSinceSoy, data=hostab)
# residuals show unequal variance between groups (which is consistent with the model)
# but there isn't a difference in central tendency of residuals. 

## look at effects of calibration error:
plot(residual ~ phbslope, data=hostab, xlab="slope used in calibrating PHB/cell")
qplot(data=hostab, x=phbslope, y=residual, geom="point") + facet_wrap(facets=hostab$yearsSinceSoy)
# calibration error explains some of the variation in residuals. 
## it appears to affect the groups similarly, since they were randomized across measurement plates
plot(phbslope ~yearsSinceSoy, data=hostab)
plot(yearsSinceSoy ~ phbslope, data=hostab, col=grey(seq(0,1,length=4)))

# see if residuals follow assumption of t distribution from the model:

#openGraph()
#png(file="hostab-residualsQQ.png",width=1000, height=800)
layout(matrix(c(1:4),ncol=2,byrow=T))
for(i in 1:4){
  r <- hostab[hostab$yearsSinceSoy==levels(hostab$yearsSinceSoy)[i],"residual"]
  qqplot(qt(ppoints(length(r)), df=nu) , r/sigma[i] , main = levels(hostab$yearsSinceSoy)[i],
         ylab="residuals (scaled to sd=1)")
  qqline(r/sigma[i], distribution = function(p){ qt(p, df=nu )} ,col="red")
} # the fit is much better now that I calculated my residuals correctly.
#dev.off()
layout(1)

#png(file="hostab-residuals-Thist.png",width=1000, height=800)
layout(matrix(c(1:4),ncol=2,byrow=T))
for(i in 1:4){
  r <- hostab[hostab$yearsSinceSoy==levels(hostab$yearsSinceSoy)[i],"residual"]
  rdens <- dt(rpts/sigma[i],df=nu)
  hist(r/sigma[i], freq=F,ylim=c(0 , max(rdens)+0.01),main=levels(hostab$yearsSinceSoy)[i])
  rpts <- seq(min(r),max(r),length=1000)
  #lines(rpts, dnorm(rpts,mean=0,sd=sigma[i] ),col="red")
  lines(rpts/sigma[i],rdens ,col="red")
}
#dev.off()
layout(1)

# looks OK overall.

##Homoscedacity: this wasn't one of the assumptions in the model, but I should 
# see if residuals are close to posterior group-level variance 
sdobs <- tapply(hostab$residual, hostab$yearsSinceSoy, sd)
show(rbind(sdobs,sigma))
plot(sdobs, sigma)
abline(a=0,b=1)
## within-group SD from model is very close to that of residuals

##Look at group-level residuals (data mean minus predicted group mean)
# note: this is different than just taking mean of estimated values.
meanEstx1 <- b0 + b1
x1means <- tapply(hostab$sqrt.phb, hostab$yearsSinceSoy,mean)
rot.resid <- x1means - meanEstx1
plot(meanEstx1,rot.resid,xlab="mean estimate (rotation trt)", ylab="level-2 residuals",ylim=c(-0.01,0.01))
plot(meanEstx1,rot.resid,xlab="mean estimate (rotation trt)", ylab="level-2 residuals",ylim=c(-0.1,0.1))
# there is some increase in residuals with higher mean estimates. But the difference is very small 
# relative to the scale of the data:
max(rot.resid) - min(rot.resid)

## look at group-level residuals for plant and plot:
meanEstx2 <- numeric(length(unique(hostab$plot)))
for(i in 1:length(unique(hostab$plot))){
  meanEstx2[i] <- b0 +  b1[x1inx2[i]] + b2[i] 
} 
meanEstx3 <- numeric(length(unique(hostab$plantID)))

for(i in 1:length(unique(hostab$plantID))){
  meanEstx3[i] <- b0 +  b1[x1inx3[i]] + b2[x2inx3[i]] + b3[i] 
} 

x2means <- tapply(hostab$sqrt.phb, hostab$plot, mean)
x2.resid <- x2means - meanEstx2
x3means <- tapply(hostab$sqrt.phb, hostab$plantID, mean)
x3.resid <- x3means - meanEstx3

openGraph()
layout(matrix(1:3, ncol=3))
plot(meanEstx1, rot.resid,ylab="level-2 residuals (rotation)",xlab="model estimate",ylim=c(-0.01,0.01))
abline(h=0,col="red")
plot(meanEstx2, x2.resid,ylab="level-2 residuals (plot)",xlab="model estimate")
abline(h=0,col="red")
plot(meanEstx3, x3.resid,ylab="level-3 residuals (plant)",xlab="model estimate")
abline(h=0,col="red")
layout(1)

# the model tends to slightly overestimate means in the low range (negative residuals)
# and underestimate means in the high range( positive residuals). 
# but the scale of level-2 residuals is small compared to the variation in parameter estimates
# (for example, the HDI for rotation mean (m1) spans a range of ~0.1, which is 10-fold larger than the range of the residuals.  )

plot(meanEstx3 ~ jitter(I(x1inx3-0.2),amount = 0.1),xlim=c(0.5,4.5))
points(meanEstx2 ~ jitter(I(x1inx2 + 0.2),amount=0.1),col="red")
points(meanEstx1 ~ c(1:4),col="blue",cex=2,pch=17)

