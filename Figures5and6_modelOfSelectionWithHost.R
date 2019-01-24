
##run script to make functions
source("selectionFunctions.R")

##I want to show the total selection differential between the low and high-PHB strain due to combined effects
# of selection in the symbiotic and free-living phase. The demographic variables and starting conditions are constant.
# the goal is to show the range of selection parameters that would re-create the patterns in my data. 

load("SROCobjectsForFig2_hostpres.RDATA")

p.constant <- 0.67 # proportionality constant for converting between change in mean PHB and change in 
# allele frequency. See script changeInMeanVsChangeInAlleleFreq.R
##use proportionality constant to get an "observed" change in allele frequency (roughly on par with experimental results)

meandiffobs <- sroccontrasts.pres$contrastsummary$meandiff
meandiffobs$deltaPH <- meandiffobs$mean/(0.67)
meandiffobs$deltaPH_high <- meandiffobs$HDIhigh/(0.67)
meandiffobs$deltaPH_low <- meandiffobs$HDIlow/(0.67)

meandiffobs_toplot <- data.frame(year=c(0,1,5,30),
                                 deltaPHbar = c(NA,meandiffobs[c(1,4,6),"deltaPH"]),
                                 deltaPHbar_low = c(NA,meandiffobs[c(1,4,6),"deltaPH_low"]),
                                deltaPHbar_high = c(NA,meandiffobs[c(1,4,6),"deltaPH_high"]))
                                 
meandiffobs_toplot$pH <- c(0.5,NA,NA,NA)
meandiffobs_toplot$pHlow <- c(0.5,NA,NA,NA)
meandiffobs_toplot$pHhigh <- c(0.5,NA,NA,NA)
for(i in 2:4){
  meandiffobs_toplot$pH[i] <- (meandiffobs_toplot$pH[i-1] + meandiffobs_toplot$deltaPHbar[i])
  meandiffobs_toplot$pHlow[i] <- (meandiffobs_toplot$pHlow[i-1] + meandiffobs_toplot$deltaPHbar_low[i])
  meandiffobs_toplot$pHhigh[i] <- (meandiffobs_toplot$pHhigh[i-1] + meandiffobs_toplot$deltaPHbar_high[i])
}
meandiffobs_toplot$pL <- 1- meandiffobs_toplot$pH
meandiffobs_toplot$pLlow <- 1-meandiffobs_toplot$pHlow
meandiffobs_toplot$pLhigh <- 1-meandiffobs_toplot$pHhigh


##---------------------------------
###Change in year 1: effects of population size in soil
##get change from year 0 to 1 at different starting conditions of population size:
N0test <- c(3,4,5,6,7)
#Nnodtest <- c(3,4,5,6)
Nnodtest <- 5
testcombos <- expand.grid(N0test,Nnodtest)
heatmaplist <- vector(mode="list",length=nrow(testcombos))
s.nodtest <- seq(0,2,length.out=100)
s.soiltest <- seq(0,1,length.out = 100)
for(i in 1:nrow(testcombos)){
  testout <- variablepars(N0=testcombos$Var1[i],Nnod=testcombos$Var2[i],overwinterR = 0.5)
  heatmaplist[[i]] <- testout$pLchange0to1
}

figX1heatmap0to1_N0vsNnod <- function(){
  layout(matrix(1:5,nrow=1,byrow=T),widths=c(1.5,1,1,1,1.25))
  plotimagewithneg.deltaPL(heatmaplist[[1]],plottitle = "",margins=c(5,5,5,1)+0.1,xat=c(0,"",0.5,"",1,"",1.5,"",""),xtext="",ssoillim = c(0,1))
  plotimagewithneg.deltaPL(heatmaplist[[2]],plottitle = "",margins=c(5,0,5,1)+0.1,xat=c(0,"",0.5,"",1,"",1.5,"",""),yaxt="n",ytext="",xtext="",ssoillim = c(0,1))
  plotimagewithneg.deltaPL(heatmaplist[[3]],plottitle = "",margins=c(5,0,5,1)+0.1,yaxt="n",ytext="",xat=c(0,"",0.5,"",1,"",1.5,"",""),ssoillim = c(0,1))
  plotimagewithneg.deltaPL(heatmaplist[[4]],plottitle = "",margins=c(5,0,5,1)+0.1,yaxt="n",ytext="",xat=c(0,"",0.5,"",1,"",1.5,"",""),xtext="",ssoillim = c(0,1))
  plotimagewithneg.deltaPL(heatmaplist[[5]],plottitle = "",margins=c(5,0,5,2)+0.1,yaxt="n",ytext="",xat=c(0,"",0.5,"",1,"",1.5,"",""),xtext="",ssoillim = c(0,1))
  par(mar=c(5,4,4,2)+0.1)
  layout(1)
}

####Figure 5 plot:
#tiff(file="year1change_effectofpopsize.tiff", width=8,height=3,units="in",res=100)
figX1heatmap0to1_N0vsNnod()
#dev.off()

################--------------

###Plots for Fig. 6A
s.soiltest_LM <- seq(0,0.25,length.out = 100)
s.soiltest_HM <- seq(0,1,length.out = 100)
s.nodtest <- seq(0,3.3,length.out=1000)

# this function runs the model over a range of selection paramters
pLchange_lowmortality <- variablepars(s.soilmax=max(s.soiltest_LM),dimsoil=100,dimnod=1000, s.nodmax=max(s.nodtest), N0=3, Nnod=6, overwinterR = 0.9,pL0=0.5)
pLchange_highmortality <- variablepars(s.soilmax=max(s.soiltest_HM),dimsoil=100,dimnod=1000,s.nodmax=max(s.nodtest),N0=3, Nnod=6,overwinterR = 0.2,pL0=0.5)

summary(as.vector(round(pLchange_lowmortality$pLchange0to1,3)))
summary(as.vector(pLchange_highmortality$pLchange0to1))

#####----------------plots
##any paramters that lead to fixation?
pltest <- pLchange_highmortality$pL.out #high mortality scenario
Lfix_highM <- which(round(pltest[,,],5)==1.0000,arr.ind = T)
Hfix_highM <- which(round(pltest[,,],5)==0.0000,arr.ind = T)
colnames(Lfix_highM) <- c("s.nod","s.soil","year")
colnames(Hfix_highM) <- c("s.nod","s.soil","year")
summary(Lfix_highM)
summary(Hfix_highM)

pltest <- pLchange_lowmortality$pL.out #no fixation in the low mortality scenario. 
Lfix_lowM <- which(round(pltest[,,],5)==1.0000,arr.ind = T)
Hfix_lowM <- which(round(pltest[,,],5)==0.0000,arr.ind = T)
colnames(Lfix_lowM) <- c("s.nod","s.soil","year")
colnames(Hfix_lowM) <- c("s.nod","s.soil","year")
summary(Lfix_lowM)# low-allele doesn't reach fixation
summary(Hfix_lowM)

##Get set of parameter values that produce similar allele frequency changes as the "data"
meandiffobs_toplot
#low mortality population
lowMchange0to1 <- as.matrix(which(pLchange_lowmortality$pLchange0to1 >= 0.286 &
                                    pLchange_lowmortality$pLchange0to1 <= 0.287, arr.ind = T))
lowMchange1to5 <- as.matrix(which(pLchange_lowmortality$pLchange1to5 >= -0.001 &
                                    pLchange_lowmortality$pLchange1to5 <= 0.001, arr.ind = T)) 
lowMchange5to30 <- as.matrix(which(pLchange_lowmortality$pLchange5to30 >= -0.09 &
                                    pLchange_lowmortality$pLchange5to30 <= -0.08, arr.ind = T)) 
###high mortality population
highMchange0to1 <- as.matrix(which(pLchange_highmortality$pLchange0to1 >= 0.286 &
                                    pLchange_highmortality$pLchange0to1 <= 0.287, arr.ind = T))
highMchange1to5 <- as.matrix(which(pLchange_highmortality$pLchange1to5 >= -0.001 &
                                    pLchange_highmortality$pLchange1to5 <= 0.001, arr.ind = T)) 
highMchange5to30 <- as.matrix(which(pLchange_highmortality$pLchange5to30 >= -0.09 &
                                     pLchange_highmortality$pLchange5to30 <= -0.08, arr.ind = T)) 


##exclude parameter values that lead to fixation
library(prodlim)
A <- highMchange5to30
B <- rbind(Hfix_highM[,1:2],Lfix_highM[,1:2])
mrows <- apply(A, 1, function(x){
  row.match(x,B)
})
highMchange5to30 <- A[which(is.na(mrows)),]


#tiff("Figure6_selectionWithSoylowMort_heatmap0to1.tiff",width=4,height=4, units="in",res=100)
plotimagewithneg.deltaPL(X= s.nodtest,Y=s.soiltest_LM,
                         plmat = pLchange_lowmortality$pLchange0to1,"",
                         ssoillim = c(0,max(s.soiltest_LM)),snodlim=c(0,max(s.nodtest)))
points(s.nodtest[lowMchange0to1[,1]],s.soiltest_LM[lowMchange0to1[,2]],type="l",lty=3,lwd=3,col=grey(1,alpha=1))
points(s.nodtest[lowMchange1to5[,1]],s.soiltest_LM[lowMchange1to5[,2]],type="l",lty=3,lwd=3,col=grey(0.5,alpha=1))
points(s.nodtest[lowMchange5to30[,1]],s.soiltest_LM[lowMchange5to30[,2]],type="l",lty=3,lwd=3,col=grey(0,alpha=1))
demopoints1 <- c(3.2,0.14)
points(demopoints1[1],demopoints1[2],pch=8,cex=1.5,lwd=2,col=grey(0,alpha=1))
#dev.off()

#tiff("Figure6_selectionWithSoylowMort_heatmap1to5.tiff",width=4,height=4, units="in",res=100)
plotimagewithneg.deltaPL(X= s.nodtest,Y=s.soiltest_LM,
                         plmat = pLchange_lowmortality$pLchange1to5,"",
                         ssoillim = c(0,max(s.soiltest_LM)),snodlim=c(0,max(s.nodtest)))
points(s.nodtest[lowMchange0to1[,1]],s.soiltest_LM[lowMchange0to1[,2]],type="l",lty=3,lwd=3,col=grey(1,alpha=1))
points(s.nodtest[lowMchange1to5[,1]],s.soiltest_LM[lowMchange1to5[,2]],type="l",lty=3,lwd=3,col=grey(0.5,alpha=1))
points(s.nodtest[lowMchange5to30[,1]],s.soiltest_LM[lowMchange5to30[,2]],type="l",lty=3,lwd=3,col=grey(0,alpha=1))
points(demopoints1[1],demopoints1[2],pch=8,cex=1.5,lwd=2,col=grey(0,alpha=1))
#dev.off()

#tiff("Figure6_selectionWithSoylowMort_heatmap5to30.tiff",width=4,height=4, units="in",res=100)
plotimagewithneg.deltaPL(X= s.nodtest,Y=s.soiltest_LM,
                         plmat = pLchange_lowmortality$pLchange5to30,"",
                         ssoillim = c(0,max(s.soiltest_LM)),snodlim=c(0,max(s.nodtest)))
points(s.nodtest[lowMchange0to1[,1]],s.soiltest_LM[lowMchange0to1[,2]],type="l",lty=3,lwd=3,col=grey(1,alpha=1))
points(s.nodtest[lowMchange1to5[,1]],s.soiltest_LM[lowMchange1to5[,2]],type="l",lty=3,lwd=3,col=grey(0.5,alpha=1))
points(s.nodtest[lowMchange5to30[,1]],s.soiltest_LM[lowMchange5to30[,2]],type="l",lty=3,lwd=3,col=grey(0,alpha=1))
points(demopoints1[1],demopoints1[2],pch=8,cex=1.5,lwd=2,col=grey(0,alpha=1))
#dev.off()

#tiff("Figure6_selectionWithSoyhighMort_heatmap0to1.tiff",width=4,height=4, units="in",res=100)
plotimagewithneg.deltaPL(X= s.nodtest,Y=s.soiltest_HM,
                         plmat = pLchange_highmortality$pLchange0to1,"",
                         ssoillim = c(0,max(s.soiltest_HM)),snodlim=c(0,max(s.nodtest)))
points(s.nodtest[highMchange5to30[,1]],s.soiltest_HM[highMchange5to30[,2]],type="l",lty=3,lwd=3)
points(s.nodtest[highMchange1to5[,1]],s.soiltest_HM[highMchange1to5[,2]],type="l",lty=3,lwd=3,col=grey(0.5,alpha=1))
points(s.nodtest[highMchange0to1[,1]],s.soiltest_HM[highMchange0to1[,2]],type="l",lty=3,lwd=3,col=grey(1,alpha=1))
demopoints2 <- c(1,0.69) # pick parameter values that produce nearly no change between years 1 and 30
demopoints3 <- c(3,0.08) # pick parameter values that produce 30% increase in low-PHB allele in year 1
points(demopoints2[1],demopoints2[2],pch=8,cex=1.5,lwd=2,col=grey(0,alpha=1))
points(demopoints3[1],demopoints3[2],pch=7,cex=1.5,lwd=2,col=grey(0,alpha=1))

#dev.off()

#tiff("Figure6_selectionWithSoyhighMort_heatmap1to5.tiff",width=4,height=4, units="in",res=100)
plotimagewithneg.deltaPL(X= s.nodtest,Y=s.soiltest_HM,
                         plmat = pLchange_highmortality$pLchange1to5,"",
                         ssoillim = c(0,max(s.soiltest_HM)),snodlim=c(0,max(s.nodtest)))
points(s.nodtest[highMchange5to30[,1]],s.soiltest_HM[highMchange5to30[,2]],type="l",lty=3,lwd=3)
points(s.nodtest[highMchange1to5[,1]],s.soiltest_HM[highMchange1to5[,2]],type="l",lty=3,lwd=3,col=grey(0.5,alpha=1))
points(s.nodtest[highMchange0to1[,1]],s.soiltest_HM[highMchange0to1[,2]],type="l",lty=3,lwd=3,col=grey(1,alpha=1))
points(demopoints2[1],demopoints2[2],pch=8,cex=1.5,lwd=2,col=grey(0,alpha=1))
points(demopoints3[1],demopoints3[2],pch=7,cex=1.5,lwd=2,col=grey(0,alpha=1))

#dev.off()

#tiff("Figure6_selectionWithSoyhighMort_heatmap5to30.tiff",width=4,height=4, units="in",res=100)
plotimagewithneg.deltaPL(X= s.nodtest,Y=s.soiltest_HM,
                         plmat = pLchange_highmortality$pLchange5to30,"",
                         ssoillim = c(0,max(s.soiltest_HM)),snodlim=c(0,max(s.nodtest)))
points(s.nodtest[highMchange5to30[,1]],s.soiltest_HM[highMchange5to30[,2]],type="l",lty=3,lwd=3)
points(s.nodtest[highMchange1to5[,1]],s.soiltest_HM[highMchange1to5[,2]],type="l",lty=3,lwd=2,col=grey(0,alpha=0.5))
points(s.nodtest[Lfix_highM[,1]],s.soiltest_HM[Lfix_highM[,2]],pch=19,col=grey(0.5,alpha=0.5))
points(s.nodtest[Hfix_highM[,1]],s.soiltest_HM[Hfix_highM[,2]],pch=19,col=grey(0.5,alpha=0.5))
points(s.nodtest[highMchange0to1[,1]],s.soiltest_HM[highMchange0to1[,2]],type="l",lty=3,lwd=3,col=grey(1,alpha=1))
points(demopoints2[1],demopoints2[2],pch=8,cex=1.5,lwd=2,col=grey(0,alpha=1))
points(demopoints3[1],demopoints3[2],pch=7,cex=1.5,lwd=2,col=grey(0,alpha=1))
#dev.off()

##summary of allele frequency change in first year under conditions that produce a 
# near steady-state after the first year
summary(as.vector(
  pLchange_highmortality$pLchange0to1[highMchange5to30[2,1],highMchange5to30[2,2]]))
summary(as.vector(
  pLchange_highmortality$pLchange0to1[highMchange1to5[,1],highMchange1to5[,2]]))
##77,95
pLchange_highmortality$pLchange0to1[77,95]


figX1colorkey <- function(margins=c(5,4,4,10)+0.1,axiscex= 1.3, labcex=1.4, line=1,ylim=c(-0.7,0.5),las=1){
  key <- matrix(seq(ylim[1],ylim[2],by=0.01),ncol=1)
  par(mar=margins)
  image(t(key),col=NA,xaxt="n",yaxt="n")
  yticks <- seq(min(ylim),max(ylim),by=0.1)
  axis(2,at=seq(1,length(key),length.out=length(yticks))/length(key), round(yticks,1),cex.axis=axiscex,las=las)
  image(t(key), zlim=c(-0.7,0),col=negpallet(50),add=T)
  image(t(key), zlim=c(0,0.5),col=pospallet(50),add=T)
  #labkey <- expression(paste(Delta,"p"["L"]," year 0 to 1"))
  #mtext(side=3,at=1,line=1,text=labkey,cex=1.2,xpd=NA)
  par(mar=c(5,4,4,2)+0.1)
}


#tiff("Figure6_selectionWithSoy_colors2.tiff",width=3,height=5, units="in",res=100)
  figX1colorkey(las=0)
#dev.off()

#### show allele frequency, selection, and population size for parameter values marked on heatplots.  
LM.outfig <- selectionWithSoy(Nsoil0 = 3,logNnod = 6,overwinterR = 0.9,s.nodL = demopoints1[1],
                                   s.soilH = demopoints1[2])
HM.outfig <- selectionWithSoy(Nsoil0 = 3,logNnod = 6,overwinterR = 0.2,s.nodL = demopoints2[1],
                             s.soilH = demopoints2[2])
HM.outfig2 <- selectionWithSoy(Nsoil0 = 3,logNnod = 6,overwinterR = 0.2,s.nodL = demopoints3[1],
                              s.soilH = demopoints3[2])
plot(0:30,LM.outfig$pL,ylim=c(0,1),pch=8)

plot(0:30,HM.outfig$pL,ylim=c(0,1),pch=8)
points(meandiffobs_toplot$year,meandiffobs_toplot$pL,type="b",cex=2)

##just out of curiosity, what happens if I make overwinter survival stochastic? (% survival poisson-distributed)
iter = 1e3
LM.outfig_stochasticPL <- c()
LM.outfig_stochasticNsoil <- c()
for(i in 1:iter){
  zz <- selectionWithSoy_stochastic(Nsoil0 = 3,logNnod = 6,overwinterR = 0.9,s.nodL = demopoints1[1],
                                                      s.soilH = demopoints1[2])  
  LM.outfig_stochasticPL <- cbind(LM.outfig_stochasticPL,zz$pL)
  LM.outfig_stochasticNsoil <- cbind(LM.outfig_stochasticNsoil, zz$Nsoil)
}

HM.outfig_stochasticPL <- c()
HM.outfig_stochasticNsoil <- c()
for(i in 1:iter){
  zz <- selectionWithSoy_stochastic(Nsoil0 = 3,logNnod = 6,overwinterR = 0.2,s.nodL = demopoints2[1],
                                    s.soilH = demopoints2[2])  
  HM.outfig_stochasticPL <- cbind(HM.outfig_stochasticPL,zz$pL)
  HM.outfig_stochasticNsoil <- cbind(HM.outfig_stochasticNsoil, zz$Nsoil)
}

matplot(0:30,LM.outfig_stochasticPL, pch=19,col="black",type="b",ylim=c(0,1))
matplot(0:30,LM.outfig_stochasticNsoil, pch=19,col="black",type="b")

matplot(0:30, HM.outfig_stochasticPL,pch=10, col="black",ylim=c(0,1))
###doesn't really make a big difference

#################
  
figX1selection <- function(m.outfig = HM.outfig,m.outfig2=NULL,m.pch=8,ylim.s=c(-0.11,0.7)){
############################### 
  ## plot of allele frequency and selection differential when selection in nodules goes down after year 1
  par(mar=c(5,4,2,5)+0.1)
  plot(0:30,m.outfig$pL,xlab="years of soybean monoculture", cex.axis=1.3,type="b",ylim=c(0,1),
       ylab="",pch=m.pch,lwd=2,bg=grey(0,alpha=0.5),cex=1.5,cex.lab=1.4)
  if(!is.null(m.outfig2)){
    points(0:30,m.outfig2$pL,pch=7,type="b",cex=1.5,lwd=2)
  }
  points(meandiffobs_toplot$year,meandiffobs_toplot$pL,type="b",cex=2,lwd=2)
  mtext(side=2,text=expression(paste("frequency of low-PHB allele ( ","p"["L"],")")),cex=1.2,line=2)
  par(new=T)
  plot(0:30,m.outfig$s.tot,type="b",xlab="",ylab="",yaxt="n",xaxt="n",
       cex=1.5,lty=2,pch=m.pch,col=grey(0,alpha=0.5),lwd=2,ylim=ylim.s)
  axis(4,cex.axis=1.3)
  mtext(side=4,text="selection coefficient (L vs. H)",col=grey(0,alpha=0.5),line=2,cex=1.2)
  par(mar=c(5,4,4,2)+0.1)
  }


figX1Nsoil <- function(m.outfig=HM.outfig,ylim.n=c(10^3,10^7),yinc=1e5){  
#########################-----------------------
  ###corresponding plot of population size in soil
  par(mar=c(5,4,2,5)+0.1)
  plot(0:30,m.outfig$Nsoil,type="l",ylab="",yaxt="n",xlab="years of soybean monoculture",
       cex=1.5,lty=2,pch=22,col=grey(0,alpha=0.5),lwd=2,cex.lab=1.4,cex.axis=1.3,ylim=ylim.n)
  points(0:30,(m.outfig$Nsoil)*m.outfig$pL,type="b",cex=1.5,lwd=2,pch=22,bg=grey(0,alpha=0.5),col=grey(0,alpha=1))
  points(0:30,(m.outfig$Nsoil)*(1-m.outfig$pL),type="b",cex=1.5,lwd=2,pch=24,bg=grey(0,alpha=0.5),col=grey(0,alpha=1))
  axis(2, cex.axis=1.3)
  #axis(2, at=c(10^3,seq(1e5, 1e7,by=1e5)),labels=rep("",101))
  #axis(2, at=c(10^3,seq(1e5, 1e6,by=yinc)),cex.axis=1.3)
  mtext(side=2,cex=1.2,text=expression("rhizobia g"^-1*soil),line=2)
  layout(1)
  par(mar=c(5,4,4,2)+0.1)
}


#tiff("Figure6_selectionWithSoy_selectionHM.tiff",width=5,height=4, units="in",res=100)
figX1selection(m.outfig=HM.outfig,m.outfig2=HM.outfig2,m.pch=8,ylim.s = c(-0.0034,0.7))
#dev.off()

#tiff("Figure6_selectionWithSoy_selectionLM.tiff",width=5,height=4, units="in",res=100)
figX1selection(m.outfig=LM.outfig,m.pch=8,ylim.s = c(-0.01,2.7))
#dev.off()


#tiff("Figure6_selectionWithSoy_NsoilHM.tiff",width=5,height=4, units="in",res=100,bg="transparent")
figX1Nsoil(m.outfig = HM.outfig,ylim.n=c(10^3,10^5.55))
#dev.off()

#tiff("Figure6_selectionWithSoy_NsoilLM.tiff",width=5,height=4, units="in",res=100)
figX1Nsoil(m.outfig = LM.outfig,ylim.n=c(10^3,10^7))
#dev.off()



###how does the change in allele frequency (for high-mortality population) translate to a change in 
# mean for fake-but-plausible PHB data?
## just look at change over the first year:
library(truncnorm)
Ntot <- HM.outfig$Nsoil[1:2]
pL <- HM.outfig$pL[1:2]
simpop <- vector(mode="list",length=2)
set.seed(57)
for(i in 1:2){
 Nh <-ceiling(Ntot[i]*(1-pL[i]))
 Nl <- ceiling(Ntot[i]*pL[i])
 cv=0.6
  high.sim <- rtruncnorm(Nh,mean=0.7,sd = cv*0.7,a=0,b=Inf)
  low.sim <- rtruncnorm(Nl, mean=0.07, sd = cv*0.07,a=0,b=Inf)
  simpop[[i]] <- sample(c(low.sim,high.sim),replace=F)
}
plot(density(simpop[[1]]))
points(density(simpop[[2]]),type="l")
##what's the difference in mean (from square-root transformed data)?
m1 <- mean(sqrt(simpop[[1]]))
m2 <- mean(sqrt(simpop[[2]]))
sd1 <-  sd(sqrt(simpop[[1]]))
sd2 <- sd(sqrt(simpop[[2]]))
sdpooled = sqrt((sd1^2 + sd2^2)/2)
##what's the effect size?
(m2-m1)/sdpooled

## effect size for change from year 0 to year 1 is -0.0875--but that's for 
# the whole population. For a smaller sample, it would be different. 
##does total population size make a difference (when it's high)?
##I could actually figure out the relationship between change in mean and change in allele frequency:


##make graphs for defense slides:
  ############################### 
  ## plot of allele frequency and selection differential when selection in nodules goes down after year 1
 

baseplot <- function(){ 
  par(mar=c(5,5,3,2)+0.1)
  plot(meandiffobs_toplot$year,meandiffobs_toplot$pL,type="b",cex=2,lwd=2,xlab="years of soybean monoculture", 
       cex.axis=1.3,ylim=c(0,1),ylab="",pch=21,cex.lab=1.4,las=1,lty=2)
  mtext(side=2,text=expression(paste("frequency of low-PHB genotype ( ","p"["L"],")")),cex=1.2,line=2.7)
}
baseplot()
points(0:30,HM.outfig$pL,type="l",lwd=2,col=rgb(1,0,0,alpha=1))
points(0:30,HM.outfig2$pL,type="l",lwd=2,col=rgb(1,0,0,alpha=1))
points(meandiffobs_toplot$year,meandiffobs_toplot$pL, type="b",lwd=2,cex=2,lty=2)

#png(file="defenseSlides_selectionInSoy1.png",width=5,height=4,units="in",res=200)
baseplot()
#dev.off()

#png(file="defenseSlides_selectionInSoy3.png",width=5,height=4,units="in",res=200)
baseplot()
points(0:30,HM.outfig$pL,type="l",lwd=2,col=rgb(1,0,0,alpha=1))
points(0:30,HM.outfig2$pL,type="l",lwd=2,col=rgb(1,0,0,alpha=1))
points(meandiffobs_toplot$year,meandiffobs_toplot$pL, type="b",lwd=2,cex=2,lty=2)
#dev.off()

#png(file="defenseSlides_selectionInSoy2.png",width=5,height=4,units="in",res=200)
baseplot()
points(0:30,LM.outfig$pL,type="l",lwd=2,col=rgb(0,0,1,alpha=1))
points(meandiffobs_toplot$year,meandiffobs_toplot$pL, type="b",lwd=2,cex=2,lty=2)
#dev.off()

#png(file="defenseSlides_selectionInSoy4.png",width=5,height=4,units="in",res=200)
baseplot()
points(0:30,LM.outfig$pL,type="l",lwd=2,col=rgb(0,0,1,alpha=1))
points(meandiffobs_toplot$year,meandiffobs_toplot$pL, type="b",lwd=2,cex=2,lty=2)
points(0:30,HM.outfig$pL,type="l",lwd=2,col=rgb(1,0,0,alpha=1))
points(0:30,HM.outfig2$pL,type="l",lwd=2,col=rgb(1,0,0,alpha=1))
points(meandiffobs_toplot$year,meandiffobs_toplot$pL, type="b",lwd=2,cex=2,lty=2)
#dev.off()


###now I need to look at frequency-dependent selection
pL <- seq(0,1,by=0.01)


Snodmax =3
a=40
c=0.99
y = Snodmax*c*exp(-a*(pL-0.5))/(c*exp(-a*(pL-0.5))+(1-c))

#tiff("Figure7_freqDependentNoduleSL.tiff",width=4, height=4,units="in",res=100)
plot(pL, y,ylim=c(0,3),type="l",lwd=2,cex.lab=1.2,cex.axis=1.3,
     ylab="",xlab="% of nodules with low-PHB genotype",xaxt="n")
axis(1,at=seq(0,1,by=0.1),labels=seq(0,100,by=10),cex.axis=1.3)
mtext(2,text=expression(paste("host sanctions ","(nodule S"[" L"],")")),cex=1.4,line=2.5)
#dev.off()
#abline(v=0.5)


test <- selectionWithSoy(pLsoil0=0.5,logNnod=6,overwinterR = 0.5,s.soilH=0.03,
                         freqDependentSanctions = TRUE,Smax=3,a=40,c=0.99)

#tiff("Figure7_freqDependentSelectionInSoyPL.tiff",width=5,height=4,units="in",res=100)
figX1selection(m.outfig=test, m.pch=22,ylim.s=c(0,3))
#dev.off()


#png(file="defenseSlides_selectionInSoyFreqDep.png",width=5,height=4,units="in",res=200)
baseplot()
points(0:30,test$pL,type="l",lwd=3,col=rgb(0,0.7,0.3))
points(meandiffobs_toplot$year,meandiffobs_toplot$pL, type="b",lwd=2,cex=2,lty=2)
#dev.off()



#tiff("Figure7_freqDependentSelectionInSoyN.tiff",width=5,height=4,units="in",res=100)
figX1Nsoil(m.outfig=test, ylim.n=c(10^3, 10^6))
#dev.off()