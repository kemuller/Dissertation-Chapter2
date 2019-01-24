###show PHB data with simulated posteriors from JAGS model
load("SROCobjectsForFig2.RDATA")
load("SROCobjectsForFig2_hostpres.RDATA")
load("LTARNobjectsForFig2.RDATA")


##------------------------------------------------
##Now make the main plot for Fig. 2 (squared so it's a the same scale as the data)


## graph of host absence (model estimates):
figX_hostabestimates <- function(ytext=expression(paste("mean PHB phenotype (pg cell"^-1,")")),
                                                        axis.cex=1.3, lab.cex=1.4,dbty="o"){
  par(mar=c(5,5,1,2)+0.1)
  hdi <- sroc.sim$x1HDI95
  ypt <- sroc.sim$x1means
  yrange <- c(min(hdi[1,]-0.02),max(hdi[2,]+0.01))^2
  plot(1:4, ypt^2,ylim=yrange,cex=2,pch=19,xaxt="n", ylab=ytext,
     xlab="years since last host crop",cex.lab=lab.cex, yaxt="n",bty=dbty)
  axis(1,at=1:4,labels=levels(srocdf$yearsSinceSoy),cex.axis=axis.cex)
  axis(2, at=seq(0.1,0.5,by=0.1),cex.axis=axis.cex)
  segments(1:4, hdi[1,]^2, 1:4,hdi[2,]^2,lwd=2)
  #mtext(side=1, at=c(1:4),text=c("a","b","c","c"),cex=1.5,line=-1)
  par(mar=c(5,4,4,2)+0.1)
}
figX_hostabcontrasts <- function(){
  ##effect size (host absence)
  conteff <-sroccontrasts$contrastsummary$effectsize
  rope <- c(conteff$rope1[1],conteff$rope2[1])
  grp1 <- c("0","0","0","1","1","5")
  grp2 <- c("1","5","30","5","30","30")
  par(mar=c(5,5,1,2))
  plot(conteff$mean,6:1,xlim=c(-0.2, 1),yaxt="n",
     ylab="",xlab="effect size",type="n",cex.lab=1.4, cex.axis=1.3,ylim=c(0.5,6.5))
  axis(2,at=6:1, labels=rep("",6))
  rect(ybottom=0,ytop=7,xleft = rope[1], xright=rope[2],col=grey(0.5,alpha=0.5))
  abline(v=0,lty="dashed",lwd=2)
  points(conteff$mean,6:1,cex=2,pch=19)
  with(conteff, segments(HDIlow,6:1, HDIhigh,6:1,lwd=2))
  #text(1.2,6:1, paste0(round(conteff$percentGtROPE,1),"%"),cex=1.3)
 #text(0.9,6.5,labels= "% above ROPE",cex=1.3,font=2)
  xl <- par()$usr[1] 
  text(xpd=NA,xl-0.15,c(6:1),grp2,cex=1.5)
  text(xpd=NA,xl-0.4,c(6:1), grp1,cex=1.5)
  arrows(rep(xl-0.35),6:1,rep(xl-0.25,6),6:1,xpd=NA,length=0.15,lwd=3)
  text(xl-0.35,7,labels="contrasts",xpd=NA,cex=1.5,font=2)
  par(mar=c(5,5,4,2)+0.1)
 }
figX_hostpresestimates <- function(ytext=expression(paste("mean PHB phenotype (pg cell"^-1,")"))){
  #### graph of host presence (model estimates)
  par(mar=c(5,5,1,2)+0.1)
  hdi <- sroc.sim.pres$x1HDI95
  ypt <- sroc.sim.pres$x1means
  yrange <- c(min(hdi[1,]-0.02),max(hdi[2,]+0.01))^2
  plot(1:4, ypt^2,ylim=yrange,cex=2,pch=19,xaxt="n", ylab=ytext,
     xlab="years with host crop",cex.lab=1.4, yaxt="n")
  axis(1,at=1:4,labels=levels(srocdf.pres$yearsWithSoy),cex.axis=1.3)
  axis(2, at=seq(0.1,0.5,by=0.1),cex.axis=1.3)
  #mtext(side=2,at=0.65,line=2.5,text="mean PHB phenotype (pg/cell)",xpd=NA,cex=1.4)
  segments(1:4, hdi[1,]^2, 1:4,hdi[2,]^2,lwd=2)
  par(mar=c(5,4,4,2)+0.1)
}
figX_hostprescontrasts <- function(){
  ##effect size (host presence) 
  par(mar=c(5,5,1,2))
  conteff <-sroccontrasts.pres$contrastsummary$effectsize
  grp1 <- c("0","0","0","1","1","5")
  grp2 <- c("1","5","30","5","30","30")
  plot(conteff$mean,6:1,xlim=c(-1, 0.55),yaxt="n",
     ylab="",xlab="effect size",type="n",cex.lab=1.4, cex.axis=1.3,ylim=c(0.5,6.5))
  axis(2,at=6:1, labels=rep("",6))
  rect(ybottom=0,ytop=7,xleft = conteff$rope1[1], xright=conteff$rope2[1],col=grey(0.5,alpha=0.5))
  abline(v=0,lty="dashed",lwd=2)
  points(conteff$mean,6:1,cex=2,pch=19)
  with(conteff, segments(HDIlow,6:1, HDIhigh,6:1,lwd=2))
  #text(-1.15,6:1, paste0(round(conteff$percentLtROPE,1),"%"),cex=1.3)
  #text(-0.9,6.5,labels= "% below ROPE",cex=1.3,font=2)
  xl <- par()$usr[1] 
  text(xpd=NA,xl-0.15,c(6:1),grp2,cex=1.5)
  text(xpd=NA,xl-0.45,c(6:1), grp1,cex=1.5)
  arrows(rep(xl-0.4),6:1,rep(xl-0.25,6),6:1,xpd=NA,length=0.15,lwd=3)
  text(xl-0.35,7,labels="contrasts",xpd=NA,cex=1.5,font=2)
  par(mar=c(5,4,4,2)+0.1)
  #layout(1)
}

#tiff(file="Figure2_SROChostab.tiff", width=9,height=4,units="in",res=100)
layout(matrix(1:2,ncol=2))
figX_hostabestimates()
figX_hostabcontrasts()
layout(1)
#dev.off()

#png(file="versionForFord_SROChostab.png", width=480,height=420,units="px")
ytext = "PHB (pg/cell)"
axis.cex = 2
lab.cex = 2
dbty="n"
par(mar=c(5,5,1,2)+0.1)
hdi <- sroc.sim$x1HDI95
ypt <- sroc.sim$x1means
yrange <- c(min(hdi[1,]-0.02),max(hdi[2,]+0.01))^2
plot(1:4, ypt^2,ylim=yrange,cex=2,pch=19,xaxt="n", ylab=ytext,
     xlab="years since last host crop",cex.lab=lab.cex, yaxt="n",bty=dbty)
axis(1,at=1:4,labels=levels(srocdf$yearsSinceSoy),cex.axis=axis.cex,lwd.ticks=2)
axis(2, at=seq(0.1,0.5,by=0.1),cex.axis=axis.cex,lwd.ticks=2)
segments(1:4, hdi[1,]^2, 1:4,hdi[2,]^2,lwd=2)
box("plot",lwd=2)
par(mar=c(5,4,4,2)+0.1)
#dev.off()


#tiff(file="Figure2_SROChostpres.tiff", width=9,height=4,units="in",res=100)
layout(matrix(1:2,ncol=2))
figX_hostpresestimates()
figX_hostprescontrasts()
layout(1)
#dev.off()



ltarnplot <- function(){
 layout(matrix(1:2,ncol=2))
  par(mar=c(5,5,2,2)+0.1)
  ###mean estimates:
  hdi <- ltarn.sim$x1HDI95
  ypt <- ltarn.sim$x1means
  ytext <- expression(paste("mean PHB phenotype (pg cell "^-1,")"))
  yrange <- c(min(hdi[1,]-0.02),max(hdi[2,]+0.01))^2
  plot(1:3, ypt^2,ylim=yrange,cex=2,pch=19,xaxt="n", ylab=ytext,
     xlab="years since last host crop",cex.lab=1.4, yaxt="n",xlim=c(1,3.2))
  axis(1,at=1:3,labels=levels(ltarndf$yearsSinceSoy),cex.axis=1.3)
  axis(2, at=seq(0.2,0.4,by=0.05),cex.axis=1.3)
  segments(1:3, hdi[1,]^2, 1:3,hdi[2,]^2,lwd=2)
  yloc <- seq(0.06,-0.06,length.out=7)
  text(1.2,yloc+0.33,
     labels=c("S","C","C","S","O","C","S"),font=c(rep(1,6),2),cex=1.3)
  text(2.2, yloc +0.33,
     labels=c("S","C","C","S","O","S","C"),font=c(rep(1,6),2),cex=1.3)
  text(3.2, yloc +0.33,
     labels=c("S","C","C","S","O","C","C"),font=c(rep(1,6),2),cex=1.3)
#mtext(side=1, at=c(1:4),text=c("a","b","c","c"),cex=1.5,line=-1)

## effect size
conteff <-ltarncontrasts$contrastsummary$effectsize
grp1 <- c("0","0","1")
grp2 <- c("1","3","3")
plot(conteff$mean,3:1,xlim=c(-0.42, 0.5),yaxt="n",
     ylab="",xlab="effect size",type="n",cex.lab=1.3, cex.axis=1.3,ylim=c(0.5,3.5))
axis(2,at=3:1, labels=rep("",3))
rect(ybottom=0,ytop=4,xleft = conteff$rope1[1], xright=conteff$rope2[1],col=grey(0.5,alpha=0.5))
abline(v=0,lty="dashed",lwd=2)
points(conteff$mean,3:1,cex=2,pch=19)
with(conteff, segments(HDIlow,3:1, HDIhigh,3:1,lwd=2))
#text(0.7,3:1, paste0(round(conteff$percentGtROPE,1),"%"),cex=1.3,pos=4)
#text(0.7,3.5,labels= "% above ROPE",cex=1.3,font=2)
xl <- par()$usr[1] 
text(xpd=NA,xl-0.1,c(3:1),grp2,cex=1.5)
text(xpd=NA,xl-0.33,c(3:1), grp1,cex=1.5)
arrows(rep(xl-0.27),3:1,rep(xl-0.17,6),3:1,xpd=NA,length=0.15,lwd=3)
text(xl-0.25,3.6,labels="contrasts",xpd=NA,cex=1.5,font=2)
par(mar=c(5,4,4,2)+0.1)
}



#tiff(file="Figure2_LTARNhostab.tiff", width=9, height=4, units="in",res=100 )
ltarnplot()
#dev.off()


#tiff(file="hostabcontrasts.tiff", width=650, height=550, res=100)
grp1 <- c(0,0,0,1,1,5)
grp2 <- c(1,5,30,5,30)
contmeans <- sroccontrasts$meandiff
par(mar=c(7,5,2,2))
plot(1:6,contmeans$mean,ylim=c(-0.1, 0.3),xaxt="n", xlab="",
     ylab="difference in transformed means",type="n",cex.lab=1.4, cex.axis=1.3)
axis(1,at=1:6, labels=rep("",6))
rect(xleft=0,xright=7,ybottom = contmeans$rope1[1], ytop=contmeans$rope2[2],col=grey(0.5,alpha=0.5))
abline(h=0,lty="dashed",lwd=2)
points(1:6,contmeans$mean,cex=2,pch=19)
with(contmeans, segments(1:6,HDIlow,1:6, HDIhigh,lwd=2))
text(xpd=NA,c(1:6),-0.22, grp1,cex=2)
arrows(1:6,rep(-0.19,6),1:6,rep(-0.16,6),xpd=NA,length=0.15,lwd=3)
text(xpd=NA,c(1:6),-0.14, grp2,cex=2)
par(mar=c(5,4,4,2)+0.1)
#dev.off()

####show distributions with and without square-root transformation
plot(density(srocdf$phbEst_zerobound))
points(density(sqrt(srocdf$phbEst_zerobound)))
