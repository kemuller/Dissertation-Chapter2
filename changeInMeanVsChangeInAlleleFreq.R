## use simulated data to help interpret model of allele frequency change in terms of phenotypic change in PHB
##I've postulated that the change in mean PHB per cell could be modeled as a change in allele frequency. 
# I just realized that I could turn it back around. 
load("phbplusplantdata_final.RDATA")
load("SROCobjectsForFig2.RDATA")
load("SROCobjectsForFig2_hostpres.RDATA")
library(truncnorm)


####------------------------------------------------------------------------
## make a graph showing simulated PHB data alongside real PHB data:
srocdata <- phbdata[phbdata$rotation%in% c("C5","S1","S5","C1","Ss","Cn"),]
mu_low <- 0.07
#mu_high <- 0.7
mu_high <- 1
#cv=0.6
cv=0.5
sd_low <- mu_low * cv
sd_high <- mu_high * cv
sd_low <- 0.1
sd_high <- 0.4
figSX1 <- function(){
sdlow_plusminus <- c(mu_low-sd_low, mu_low+sd_low)
sdhigh_plusminus <- c(mu_high-sd_high, mu_high+sd_high)
plot(density(srocdata$phbEst_zerobound,from=0),lwd=2,main="",ylab="Density",
     xlab="PHB per cell (pg)",xlim=c(0,2),
     cex.axis=1.3,cex.lab=1.4)
abline(v=c(mu_low,mu_high),lty=3,lwd=2,col=grey(0.5))
arrows(x0 = c(sdlow_plusminus[1],sdhigh_plusminus[1]),x1=c(sdlow_plusminus[2],sdhigh_plusminus[2]),
         y0=0.1,y1=0.1,code=3,angle=90,length = 0.05,lwd=2,col=grey(0.5))
}

# tiff("FigSX_densityPlotForFakeVsRealData.tiff",width=4,height=3.5,units="in",res=100)
 figSX1()
# dev.off()
 
# make a version with a separate curve for S5 and C5 (for presentation)
figSX1_a <- function(){
   sdlow_plusminus <- c(mu_low-sd_low, mu_low+sd_low)
   sdhigh_plusminus <- c(mu_high-sd_high, mu_high+sd_high)
   s5 <- srocdata[srocdata$rotation=="S5",]
   c5 <- srocdata[srocdata$rotation=="C5",]
   plot(density(s5$phbEst_zerobound,from=0),lwd=2,main="",ylab="Density",
        xlab="PHB per cell (pg)",xlim=c(0,2),
        cex.axis=1.3,cex.lab=1.4,col="darkgreen")
   points(density(c5$phbEst_zerobound,from=0),type="l",lwd=2,col="orange")
   abline(v=c(mu_low,mu_high),lty=3,lwd=2,col=grey(0.5))
   arrows(x0 = c(sdlow_plusminus[1],sdhigh_plusminus[1]),x1=c(sdlow_plusminus[2],sdhigh_plusminus[2]),
          y0=0.1,y1=0.1,code=3,angle=90,length = 0.05,lwd=2,col=grey(0.5))
   legend("topright",bty="n",legend=c("S5","C5"),text.col=c("darkgreen","orange"),cex=2)
 }
 
 
#   png("FigSX_densityPlotForFakeVsRealData.png",width=6,height=5,units="in",res=300)
 par(mar=c(5,4,2,2)+0.1)
   figSX1_a()
#  dev.off()

 
 
 
####################--------------------------------------------

####use fake data to find a proportionality constant for converting between change 
#in allele frequency and change in square-root transformed mean PHB/cell (since that's what I have 
# in my statistical models from my experimental data)

##first, make a function to generate simulated population based on allele frequency. 
makefakeVp <- function(pL=0.5, N=10^5, 
                       sdhigh=sd_high,sdlow=sd_low, mulow=mu_low, muhigh=mu_high){
  Nl = ceiling(N * pL)
  Nh = N - Nl
  ##draws from 2 zero-truncated normal distributions with different mean and same CV
  high.sim <- rtruncnorm(Nh,mean=muhigh,sd = sdhigh,a=0,b=Inf)
  low.sim <- rtruncnorm(Nl, mean=mulow, sd = sdlow,a=0,b=Inf)
  simpop <-  sample(c(low.sim,high.sim),replace=F) # just mixes it up.
  return(simpop)
}

##Now generate populations to calculate the relationship between change in phenotypic mean 
# and change in allele frequency. 

pL0 =0.5 # start at 50% (the slope is the same wherever you start)
pLchange = c(seq(-0.5,0.5,by=0.1))
#keep constant high population size:
Ntot = 10^5
zchange <- numeric(length(pLchange)) # change in mean PHB phenotype (square-root transformed)

##get initial mean:
set.seed(67)
simpop0 <- makefakeVp(pL=pL0, N=Ntot)
zmean0 <- mean(sqrt(simpop0))

for(i in 1:length(pLchange)){
  pLi <- pL0 + pLchange[i]
  simpop <- makefakeVp(pL= pLi, N=Ntot ) # assume unchanging large population size
  zchange[i] <- mean(sqrt(simpop)) - zmean0
}

plot(pLchange,zchange)
lm.pL <- lm(zchange ~ pLchange)
show(coef(lm.pL))

deltaPLtodeltaZ = round(coef(lm.pL)[2],3)
##proportionality constant is -0.669 for an increase in frequency of low-PHB allele (or +0.669 for 
# frequency of high-PHB allele). I'm sure there's an analytical solution, but this is easier.

###---------------------------------------------------------------------
##Compare phenotypic change from simulated data to changes from actual data.
## meandiffobs --> deltaPL --> meandiffsim.

##mean differences from square-root transformed data (MCMC estimates)
meandiffobs_hostab <- sroccontrasts$meandiff[,c("contrast","mean","HDIlow","HDIhigh")]
meandiffobs_hostpres <- sroccontrasts.pres$meandiff[,c("contrast","mean","HDIlow","HDIhigh")]
meandiffobs_hostab$model = "hostabsence"
meandiffobs_hostpres$model = "hostpresence"
##put the two models into one dataframe for efficiency. 
meandiffobs <- rbind(meandiffobs_hostab,meandiffobs_hostpres)
##translate difference in means to difference in allele frequency based on proportionality constant
meandiffobs$deltaPL <- meandiffobs$mean/deltaPLtodeltaZ
meandiffobs$deltaPLHigh <- meandiffobs$HDIhigh/deltaPLtodeltaZ
meandiffobs$deltaPLLow <- meandiffobs$HDIlow/deltaPLtodeltaZ
##get simulated mean change from allele frequency change 
# use change from pL=0.5 (starting value doesn't affect difference in means, so long
# as final pL isn't 0 or 1).  
set.seed(90)
referencepop <- makefakeVp(pL=0.5, N=Ntot)
for(i in 1:nrow(meandiffobs)){
  deltaPL <- meandiffobs$deltaPL[i]
  deltaPLhigh <- meandiffobs$deltaPLHigh[i]
  deltaPLlow <- meandiffobs$deltaPLLow[i]
  simpop <- makefakeVp(pL=0.5+deltaPL)
  simpopHigh <- makefakeVp(pL=0.5+deltaPLhigh)
  simpopLow <- makefakeVp(pL=0.5+deltaPLlow)
  meandiffobs$meandiffsim[i] <-
    mean(sqrt(simpop)) - mean(sqrt(referencepop))
  meandiffobs$meandiffsimHigh[i] <-
    mean(sqrt(simpopHigh)) - mean(sqrt(referencepop))
  meandiffobs$meandiffsimLow[i] <-
    mean(sqrt(simpopLow)) - mean(sqrt(referencepop))
}


plot(meandiffobs$mean~meandiffobs$meandiffsim,ylim=c(min(meandiffobs$HDIlow)-0.01,max(meandiffobs$HDIhigh)+0.01))
points(meandiffobs$HDIlow,meandiffobs$meandiffsimLow,col="red")
points(meandiffobs$HDIhigh,meandiffobs$meandiffsimHigh,col="red")
segments(x0=meandiffobs$meandiffsim,x1=meandiffobs$meandiffsim,
         y0=meandiffobs$HDIlow,y1=meandiffobs$HDIhigh)
abline(a=0,b=1)

##do bootstrapped samples from simulated population to get an estimate of variability of mean estimates
# use same sample size as in data (i.e., # of nodules measured for each group). 
meandiffobs$group1 <-  c(rep("S5",3),rep("C1",2),"C5",rep("C5",3),rep("S1",2),"S5")
meandiffobs$group2 <- c("C1","C5","Cn","C5","Cn","Cn","S1","S5","Ss","S5","Ss","Ss")
table(srocdata$rotation)
meandiffobs$N1 <- 128
meandiffobs$N2 <- 128
meandiffobs[meandiffobs$group1=="Cn","N1"] <- 125
meandiffobs[meandiffobs$group2=="Cn","N2"] <- 125
meandiffobs[meandiffobs$group1=="S1","N1"] <- 127
meandiffobs[meandiffobs$group2=="S1","N2"] <- 127

nboot = 1e4
referencepop <- makefakeVp(pL=0.5, N=Ntot) # use the same reference population
meandiffsample <- matrix(NA, nrow=nrow(meandiffobs), ncol=nboot)
rownames(meandiffsample) <- paste0(meandiffobs$group1,"v",meandiffobs$group2)
set.seed(10)
for(i in 1:nrow(meandiffobs)){
  deltaPL <- meandiffobs$deltaPL[i]
  simpop <- makefakeVp(pL=0.5+deltaPL, N=Ntot)
  N1 = meandiffobs$N1[i]
  N2 = meandiffobs$N2[i]
  for(b in 1:nboot){ # take 1000 samples of n = 128 from simulated population (125 and 127 for Cn nd S1) 
    group1sample <- sample(referencepop, size=N1,replace=T)
    group2sample <- sample(simpop, size=N2,replace=T)
    meandiffsample[i,b] <- mean(sqrt(group2sample)) - mean(sqrt(group1sample))
}}

##get 95% confidence interval for bootstrapped samples from simulated population
bsci <- t(apply(meandiffsample,1,function(x){quantile(x,c(0.025,0.975))}))
meandiffobs$BSCIlow <- bsci[,1]
meandiffobs$BSCIhigh <- bsci[,2]

##make a plot showing mean and HDI of contrasts from actual data along with mean and bootstrap CI from simulated data
figSX2 <- function(){
x1 = rev(c(seq(1,12,by=2),seq(16,26,by=2)))
x2 = x1-0.75
par(mar=c(5,5,1,2))
xtext <- expression(paste(Delta," ",bar("z")))
plot(meandiffobs$mean,x1,xlim=c(min(meandiffobs$HDIlow-0.1),max(meandiffobs$HDIhigh+0.15)),
     ylim=c(0,(max(x1)+0.5)),pch=21,lwd=2,yaxt="n",ylab="",cex=1.5,
     cex.axis=1.3,cex.lab=1.4,xlab=xtext)
axis(side=2,at=x1,labels=rep("",length(x1)))
segments(y0=x1, y1=x1, x0=meandiffobs$HDIlow,
         x1=meandiffobs$HDIhigh,lwd=2)
points(meandiffobs$meandiffsim,x2,col=grey(0.5),pch=5,lwd=2,cex=1.3)
segments(y0=x2, y1=x2, x0=meandiffobs$BSCIlow,
         x1=meandiffobs$BSCIhigh,col=grey(0.5),lwd=2)
  grp1 <- rep(c("0","0","0","1","1","5"),2)
  grp2 <- rep(c("1","5","30","5","30","30"),2)
  xl <- par()$usr[1]
  yt <- par()$usr[4]
  xr <- par()$usr[2]
  ar = 0.1
  al = 0.22
  ashrinkr = 1.5 
  ashrinkl = 0.9
  text(xpd=NA,xl-ar,x1,grp2,cex=1.5)
  text(xpd=NA,xl-al,x1, grp1,cex=1.5)
  arrows(rep(xl-al*ashrinkl,length(x1)),x1,rep(xl-ar*ashrinkr,length(x1)),x1,xpd=NA,length=0.1,lwd=3)
  abline(h=max(x1)/2)
  text(xl+0.2,yt-1,labels="years since\nlast host crop\n(maize\nmonoculture)",cex=1.1,pos=1)
  text(xr-0.2,(max(x1)/2)-1,labels="years with\nhost crop\n(soybean\nmonoculture)",cex=1.1,pos=1)
  par(mar=c(5,5,4,2)+0.1)
}

#tiff(file="FigSX2_simvsrealdata.tiff",width=4,height=5,units="in",res=100)
figSX2()
#dev.off()


##-----------------------------------------------
# the last thing I want to show is a density plot with simulated data:
figSX3 <- function(){
pLref = 0.5
referencepop <- makefakeVp(pL=pLref,N=Ntot)
#colz <- c(grey(0),grey(0.5),grey(0.5))
colz <- c("black","blue","red")
ltyz <- c(1,1,1)
simpop1 = makefakeVp(pL=pLref+0.2,N=Ntot)
simpop2 = makefakeVp(pL=pLref-0.2,N=Ntot)
plot(density(simpop1,from=0),xlim=c(0,2),lwd=2,col=colz[2],lty=ltyz[2],cex.axis=1.3,
     cex.lab=1.4,xlab="PHB per cell (pg)",main="")
points(density(simpop2,from=0),type="l",lwd=2,col=colz[3],lty=ltyz[3])
points(density(referencepop,from=0),type="l",lwd=2,col=colz[1],lty=ltyz[1])
l1 <- expression(paste("p"[" Low PHB"]," = ",0.5))
l2 <- expression(paste("p"[" Low PHB"]," = ",0.7))
l3 <- expression(paste("p"[" Low PHB"]," = ",0.3))
legend("topright",legend=c(l1,l2,l3),cex=1.4,text.col=colz,bty="n")
}

# tiff("FigSX_densityPlotwithSimData.tiff",width=4,height=3.5,units="in",res=100)
 figSX3()
# dev.off()
 
 # png("FigSX_densityPlotwithSimData.png",width=6,height=5,units="in",res=300)
 par(mar=c(5,4,2,2)+0.1)
 figSX3()
 # dev.off()

#save(meandiffobs,file="meanDiffwithPL.RData")
