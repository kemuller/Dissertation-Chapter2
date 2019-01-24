####This is the code used to generate Fig. 3 illustrating selection outside the host. 
## The idea is that changes in mean PHB phenotype (shown in Fig. 2) could be modeled with the simplifying assumption
# that the population contains two alleles (low- and high-PHB), and that the change in mean PHB phenotype is proportional
# to the change in the relative abundance of high- and low-PHB alleles. Then, use population size data reported by 
# Revellin et al. 1996 to illustrate how the composition of the rhizobia population changes over time. 

load("SROCobjectsForFig2.RDATA")

###----------------------
## assume change in mean is proportional to change in allele frequency of a high vs. low PHB
# alleles
obs.means <- sroc.sim$x1means
meandiffs <- sroccontrasts$contrastsummary$meandiff # contrasts from MCMC (using square-root transformed means)
##get change in mean between time points:
changeovertime <- meandiffs[meandiffs$contrast %in% c("0to1","1to5","5to30"),"mean"]
## divide change by # of years elapsed to get change per year:
changeperyear <- changeovertime/c(1,4,25)

##Make plot for defense slides:
#png("defenseslides_changeperyear.png",width=6,height=5,units="in",res=200)
par(mar=c(5,5,3,2)+0.1)
plot(c(1,5,30), changeperyear, pch=19,cex=2,xlab="years since last host crop",ylab= expression(Delta*"mean PHB per year"),cex.lab=1.4,cex.axis=1.3)
#dev.off()

#png("defenseslides_changeperyearLine.png",width=6,height=5,units="in",res=200)
par(mar=c(5,5,3,2)+0.1)
plot(c(1,5,30), changeperyear, pch=19,cex=2,xlab="years since last host crop",ylab= expression(Delta*"mean PHB per year"),cex.lab=1.4,cex.axis=1.3)
x <- c(0:29)
tau <- -0.27
zt <- changeperyear[1]*exp(x*tau)
#names(zt) <- years
points(zt,col="black",type="l",lwd=2,lty=2)
#dev.off()


## setting exponential decay function with tau = -0.27 does a decent job appoximating
# change in mean phenotype per year.

##Now, show change in allele frequency for a population with two strains (H and L).
#Assume that the change in mean phenotype is proportional to the change in allele frequency of high-PHB strain (H)
# change in allele frequency undergoes exponential decay. Initial change in allele frequency
# depends on initial allele frequency and selection coefficient s, which is the proportional difference in 
# fitness between the high- and low-PHB strains. (s = (Wh - Wl)/Wl)
years <- 0:30

##based on fake data (changeInMeanVsChangeInAlleleFreq.R), 
# the relationship between mean change and allele frequency cheange is deltaZ = deltaP*(-/+ 0.566) 
# (negative if you're talking about low-PHB allele, positive if it's the high-PHB allele)
pHchange0 = changeperyear[1]/0.67 # initial allele frequency change in the first year is 0.14
#s <- 0.48 # starting selection coefficient
pH <- c(0.3) # starting allele frequency. 
#pHchange0 <- s*pH*(1-pH)/(1+s*pH) #allele frequency change over 1 cycle of selection (here, a year)

pHchange <- pHchange0*exp(years*tau)
names(pHchange) <- years
for(t in 2:length(years)){
  pH <- c(pH, pH[t-1] + pHchange[t-1])
}
plot(pHchange)

#png(file="defenseSlide_pH.png",width=6,height=5,units="in",res=200)
par(mar=c(5,5,3,2)+0.1)
plot(0:30,pH,type="l",lwd=2,ylim=c(0,1),ylab="genotype frequency in population",xlab="years since last host crop",cex.lab=1.4,cex.axis=1.3)
points(0:30,1-pH,type="l",lwd=2)
#dev.off()

proportionalityConstant <- mean(pHchange[1:30]/zt) # it's only one number, but there 
# may be slight rounding error.

##calculate selection coefficient over time:
st <- c()
for(t in 1:length(years)){
  st_t <- pHchange[t]/(pH[t]*(1-pH[t] - pHchange[t]))
  st <- c(st, st_t)
}

sL <- c() # selection for low-PHB allele
pL <- 1-pH
pLchange <- -pHchange
for(t in 1:length(years)){
  st_t <- (pLchange[t])/(pL[t]*(1-pL[t] - pLchange[t]))
  sL <- c(sL, st_t)
}
plot(0:30,st)


#png(file="defenseSlide_sHnohost.png",width=6,height=5,units="in",res=200)
par(mar=c(5,5,3,2)+0.1)
plot(0:30,st,type="l",ylab="selection for the high-PHB genotype",xlab="years since last host crop",cex.lab=1.4,cex.axis=1.3,lwd=2)
abline(h=0,lty=2,lwd=2,col=grey(0.5))
points(0:30, st, type="l",lwd=2)
#dev.off()

#par(new=T)
#plot(0:30,pH)
#plot(st,pHchange) # fitness difference declines over time.

## now apply those frequencies to realistic population sizes in soil (use regression from Revellin)
## let's just assume population undergoes a linear decrease (on a log-scale) until it reaches a stable population size.
## for revellin data, show points corrsponding to the slope and intercept of the regression
revslopes <- c(-0.18, -0.49)
revintercepts <- c(5.30, 4.73)
years <- c(0:29)
revyears <- list(gr1=c(1:11), gr2=c(1:7))

log10Nt1 <- revslopes[1]*years + revintercepts[1]
log10Nt2 <- revslopes[2]*years + revintercepts[2]
log10Nt2[log10Nt2<0] <- 0
Nt1 <- 10^log10Nt1
Nh1 <- (Nt1)*pH[1:length(log10Nt1)]
Nl1 <- (Nt1)*(1-pH)[1:length(log10Nt1)]

Nh2 <- (10^log10Nt2)*pH[1:length(log10Nt2)]
Nl2 <- (10^log10Nt2)*(1-pH)[1:length(log10Nt2)]

# get time point at which population size passes below "nodulation saturation" (10^2 rhizobia/g soil)
nst <- (2 - revintercepts)/revslopes
##figure out allele frequency at that point:

nst <- round(nst,3)
round(pH[18],2)
round(pH[6],2)

## Figure 3
colz <- c(grey(0.5,alpha=0.5),grey(0,alpha=0.5),grey(0,alpha=0.5))
pchz <- c(21,24,22)
lty.l= c(3,1,1)
cex.pt <- 1.5
ylabtext <- expression("rhizobia g"^-1*"soil")

fig3_hostabNsoil <- function(){
  par(mar=c(5,5,2,4)+0.1)
  plot(revyears$gr1, log10Nt1[revyears$gr1],type="n",ylim=c(1,5.4),xlab="years since last host crop",
     ylab="",yaxt="n",cex.axis=1.3,cex.lab=1.4,xaxt="n")
  mtext(ylabtext, side=2,line=2.5,cex=1.4)
  axis(2, at=c(1:5),labels=parse(text=paste0("10","^",1:5)),las=1,cex.axis=1.3)
  axis(1,at=c(1:11),cex.axis=1.3)
  abline(h=2,lty=2,lwd=2)
  points(revyears$gr1, log10Nt1[revyears$gr1],type="l",ylim=c(1,5.4),lwd=2,col=colz[1],pch=pchz[1],cex=cex.pt,lty=lty.l[1])
  points(revyears$gr2, log10Nt2[revyears$gr2], type="l",lwd=2,col=colz[1],pch=pchz[1],cex=cex.pt,lty=lty.l[1])
  points(revyears$gr1, log10(Nh1)[revyears$gr1], type="b",lwd=2,bg=colz[2],col=colz[2],pch=pchz[2],cex=cex.pt,lty=lty.l[2])
  points(revyears$gr1, log10(Nl1)[revyears$gr1], type="b",lwd=2,bg=colz[3],col=colz[3],pch=pchz[3],cex=cex.pt,lty=lty.l[3])
  points(revyears$gr2, log10(Nh2)[revyears$gr2], type="b",lwd=2,bg=colz[2],col=colz[2],pch=pchz[2],cex=cex.pt,lty=lty.l[2])
  points(revyears$gr2, log10(Nl2)[revyears$gr2], type="b",lwd=2,bg=colz[3],col=colz[3],pch=pchz[3],cex=cex.pt,lty=lty.l[3])
  par(mar=c(5,4,4,2)+0.1)
  }

fig3_hostabselection <- function(margins=c(5,5,2,4)+0.1){
  ytext <- expression(paste("freq. of high-PHB allele (p "["H"]," )"))
  par(mar=margins)
  plot(revyears$gr1,pH[revyears$gr1],type="b", cex=1.5,lwd=2,pch=24,col=grey(0,alpha=1),bg=grey(0,alpha=0.5),
     ylab=ytext,xlab="years since last host crop",cex.lab=1.4,cex.axis=1.3,xaxt="n",
     ylim=c(0,1))
  axis(1,at=c(1:11))
  par(new=T)
  plot(revyears$gr1,st[revyears$gr1], type="b",yaxt="n",xaxt="n",cex=1.5,lwd=2,ylab="",
       xlab="",col=grey(0,alpha=0.5),pch=24,lty=2,ylim=c(0,max(st)+0.1))
  axis(4,cex.axis=1.3,at=seq(0,1,by=0.1))
  mtext(side=4, text="selection coefficient (H vs. L)",line=2,cex=1.4,col=grey(0,alpha=0.5))
  par(mar=c(5,4,4,2)+0.1)
}
## put in x indicating change in mean PHB phenotype (proportional to allele frequency.

#tiff(file="Figure3selectionoutsidehost_Nsoil.tiff", width=6,height=5, units="in",res=300)
fig3_hostabNsoil()
#dev.off()

#tiff(file="Figure3selectionoutsidehost_pL.tiff", width=6,height=5, units="in",res=300)
fig3_hostabselection()
#dev.off()


##See how it looks with exponential decay function instead.

log10.0 <- 7
log10.s <- 3
years <- c(0:29)
log10.t <- log10.s + (log10.0 - log10.s)*exp(years * -0.6)
Nt <- 10^3 + (10^7-10^3)*exp(years*-0.6)
plot(years,Nt)
Nh <- pH * Nt
Nl <- Nt-Nh
points(years,Nh, col="red")
points(years,Nl, col="blue")

plot(years,log10(Nt),ylim=c(2,7))
points(years,log10(Nh),col="red")
points(years,log10(Nl),col="blue")

