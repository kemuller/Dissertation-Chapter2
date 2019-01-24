load("phbplusplantdata_final.RDATA")
load("SROCobjectsForFig2.RDATA")
load("SROCobjectsForFig2_hostpres.RDATA")
load("LTARNobjectsForFig2.RDATA")
library(plyr)
##add back in nodule count and plant mass:
srocdf <- merge(srocdf, pd[,c("plantID","","nodMassG")],all.x=T)
srocdf.pres <- merge(srocdf.pres, pd[,c("plantID","nodCount","nodMassG")],all.x=T)
ltarndf <- merge(ltarndf, pd[,c("plantID","nodCount","nodMassG")],all.x=T)
ltarnpars <- read.csv("LTARN-SummaryInfo.csv",stringsAsFactors = F)

####yhat is the model estimate of PHB per cell for rhizobia from each trap plant. It's what I'm using as the 
# mean PHB phenotype within a host plant.

hostabpd <- ddply(srocdf, .(plantID),yearsSinceSoy= unique(yearsSinceSoy),yhat=mean(yhat) ,plantnum=unique(plantnum), plotnum=unique(plotnum), shootMassG= mean(shootMassG),
                  nodMassG = mean(nodMassG), nodCount = mean(nodCount),summarize)
ltarnpd <- ddply(ltarndf, .(plantID),yearsSinceSoy= unique(yearsSinceSoy), yhat=mean(yhat),plantnum=unique(plantnum), plotnum=unique(plotnum), shootMassG= mean(shootMassG),
                 nodMassG = mean(nodMassG), nodCount = mean(nodCount),summarize)
hostprespd <- ddply(srocdf.pres, .(plantID),yearsWithSoy= unique(yearsWithSoy), yhat=mean(yhat),plantnum=unique(plantnum), plotnum=unique(plotnum), shootMassG= mean(shootMassG),
                    nodMassG = mean(nodMassG), nodCount = mean(nodCount),summarize)

####code needed for results section:

#shoot dry mass per plant was highly variable within each category of host availability (spanning a two-fold range of < 1 g/plant to >2 g per plant)

#" Variation in shoot mass was highly correlated with variation in nodule mass per plant (r = 0.93, R2 = 0.87)"
## This graph is for all data combined (160 plants). This includes the LTARN and the 30-year nematode-resistant soy group, which I 
# didn't end up using for PHB measurements.

##Shoot mass vs. nodule mass:
par(mar=c(5,4,4,2)+0.1)
plot(pd$nodMassG,pd$shootMassG,ylab="shoot dry mass (g)", xlab="nodule dry mass (g)" ,lwd=2,
     cex.axis=1.4, cex.lab=1.3,cex=1.5,col=grey(0,alpha=0.5))
nodshootlm <- lm(shootMassG ~ nodMassG, data=pd) 
abline(nodshootlm,lwd=2)
summary(nodshootlm)
cor(pd$nodMassG, pd$shootMassG)

# shoot mass vs. nodule count:
# best fit is with a cubic function
nodcountlm1 <- lm(shootMassG ~ nodCount, data=pd)
nodCountsquared <- pd$nodCount^2
nodCountcubed <- pd$nodCount^3
nodcountlm2 <- lm(shootMassG ~ nodCount + nodCountsquared, data=pd)
nodcountlm3 <- lm(shootMassG ~ nodCount + nodCountsquared + nodCountcubed, data=pd)
nodcountlm3minus2 <- lm(shootMassG ~ nodCount  + nodCountcubed, data=pd)
AIC(nodcountlm3,nodcountlm3minus2,nodcountlm2,nodcountlm1)

nodcountlm.exp <- lm(shootMassG ~ nodCount*factor(experiment), data=pd)
anova(nodcountlm.exp)

xval <- seq(1,100,by=1)
pred1 <- predict(nodcountlm1,list(nodCount = xval))
pred2 <- predict(nodcountlm2, list(nodCount=xval, nodCountsquared=xval^2))
pred3 <- predict(nodcountlm3, list(nodCount=xval, nodCountsquared=xval^2, nodCountcubed=xval^3))

summary(nodcountlm3)
summary(nodcountlm2)
summary(nodcountlm1)
cor(pd$nodCount, pd$shootMassG)
#points(ltarnpd$nodMassG,ltarnpd$shootMassG,cex.axis=1.4, cex.lab=1.3,col="red")

#tiff("supplement_shootmassVnodMassCount.tiff",width=1000,height=500,res=100)
layout(matrix(1:2, nrow=1,byrow=T))
par(mar=c(5,4,2,1)+0.1)

plot(shootMassG ~ nodMassG,ylab="shoot dry mass per plant (g)", xlab="nodule dry mass per plant (g)" ,lwd=2,
     cex.axis=1.4, cex.lab=1.3, data=pd, type="n")
points(shootMassG ~ nodMassG , data=pd[pd$experiment=="SROC",],cex=1.5,col=grey(0,alpha=0.5),lwd=2)
points(shootMassG ~ nodMassG , data=pd[pd$experiment=="LTARN",],cex=1.5,col=grey(0,alpha=0.5),lwd=2,pch=17)
abline(nodshootlm,lwd=2)
legend("topleft",c(as.expression(bquote(R^2 ~ "= 0.87"))) ,bty="n",cex=1.5)
legend("bottomright",c("y = 9.29x + 0.59") ,bty="n",cex=1.5)

par(mar=c(5,2,2,2)+0.1)
plot(pd$nodCount,pd$shootMassG,ylab="", xlab="nodules per plant" ,lwd=2,
     cex.axis=1.4, cex.lab=1.3,cex=1.5,col=grey(0,alpha=0.5),type="n") # doesn't really look like a linear relationship
points(shootMassG ~ nodCount , data=pd[pd$experiment=="SROC",],cex=1.5,col=grey(0,alpha=0.5),lwd=2)
points(shootMassG ~ nodCount , data=pd[pd$experiment=="LTARN",],cex=1.5,col=grey(0,alpha=0.5),lwd=2,pch=17)
lines(pred1,lwd=2,col=grey(0,alpha=1))
#cubic regression was better, but not that much better.
#lines(pred3,col=grey(0,alpha=0.5),lty=1,lwd=2)
#legend("topleft",c(as.expression(bquote(R^2 ~ "= 0.45")),
#                   as.expression(bquote(R^2 ~ "= 0.47")) ),text.col=c("black",grey(0, alpha =0.5)), bty="n",cex=1.5)
legend("topleft",c(as.expression(bquote(R^2 ~ "= 0.45"))), bty="n",cex=1.5)
legend("bottomright",c("y = 0.75x + 0.02") ,bty="n",cex=1.5)

par(mar=c(5,4,4,2)+0.1)
layout(1)
#dev.off()

###Just out of curiosity: is the relationship between nodule count and nodule mass different between the two experiments?
plot(nodMassG ~ nodCount,ylab="nodule dry mass per plant (g)", xlab="nodule count (g)" ,lwd=2,
     cex.axis=1.4, cex.lab=1.3, data=pd, type="n")
points(nodMassG ~ nodCount , data=pd[pd$experiment=="SROC",],cex=1.5,col=grey(0,alpha=0.5),lwd=2)
points(nodMassG ~ nodCount , data=pd[pd$experiment=="LTARN",],cex=1.5,col=grey(0,alpha=0.5),lwd=2,pch=17)
abline(nodshootlm,lwd=2)


## Plant size did not appear to contribute to variation in PHB per cell measured in rhizobia
## I was going to use model estimates for PHB, but I decided it would be simpler to use raw means and show all plants together.

pd2 <- ddply(phbdata,.(plantID), meansqrtphb = mean(sqrt(phbEst_zerobound)),shootMassG = mean(shootMassG), 
             nodMassG = mean(nodMassG),nodCount = mean(nodCount), summarize)
pd2$meanPHB <- pd2$meansqrtphb^2
phbvshootmodel <- lm(meanPHB ~ shootMassG, data=pd2)
phbvnodcountmodel <- lm(meanPHB ~ nodCount, data=pd2)

summary(phbvshootmodel)
summary(phbvnodcountmodel)

#tiff("supplement_shootmassVmeanPHB.tiff",width=1000,height=500,res=100)
layout(matrix(1:2,nrow=1))
par(mar=c(5,4,2,1)+0.1)
plot(pd2$shootMassG, pd2$meanPHB, cex.lab=1.4,cex.axis=1.3,cex=1.5,lwd=2,col=grey(0, alpha=0.5),
     ylab="pg PHB/cell (mean per host plant)",xlab="shoot dry mass per plant (g)")
abline(phbvshootmodel,lwd=2)
par(mar=c(5,2,2,2)+0.1)

plot(pd2$nodCount, pd2$meanPHB, cex.lab=1.4,cex.axis=1.3,cex=1.5,lwd=2,col=grey(0, alpha=0.5),
     ylab="",xlab="nodules per plant")
abline(phbvnodcountmodel,lwd=2)
par(mar=c(5,4,4,2)+0.1)
layout(1)
#dev.off()

#Trap plants did show variation in nodule number associated with soil source {supplement},
#suggesting differences in rhizobia population size among rotation treatments. 
#Plants inoculated with soil from the 30-year maize monoculture produced  {57% to 76%} ~10 fewer nodules 
#per plant, on average, than plants inoculated with soil from other plots within the same experiment. 
#In the LTARN experiment, soil that had not had soybean for 1-3 years produced, on average,16 fewer nodules
#per plant than soil where soybean was grown the previous year.

pdsroc <- pd[pd$experiment=="SROC" & pd$rotation != "Sr",]
pdsroc$rotation <- factor(pdsroc$rotation,levels=c("Ss","S5","S1","C1","C5","Cn"))

pdsroc$plot <- factor(pdsroc$plot)

nodcountmodel1 <- glm(nodCount ~ rotation , data=pdsroc, family="quasipoisson")
nodcountmodel2 <- glm(nodCount ~ yearsSinceSoy , data=ltarnpd, family="quasipoisson")
summary(nodcountmodel1)
summary(nodcountmodel2)
anova(nodcountmodel1,test="Chisq")
anova(nodcountmodel2,test="Chisq")

nodmassmodel1 <- lm(nodMassG ~ rotation, data=pdsroc)
nodmassmodel2 <- lm(nodMassG ~ yearsSinceSoy, data=ltarnpd)
summary(nodmassmodel1)
summary(nodmassmodel2)


##Residual deviance is much bigger than degrees of freedom--so, quasipoisson.
library(multcomp)
nodcountTukey1 <- summary(glht(nodcountmodel1, mcp(rotation="Tukey")))
nodcountTukey2 <- summary(glht(nodcountmodel2, mcp(yearsSinceSoy="Tukey")))
nodcountTukey1
nodcountTukey2

preddat <- predict(nodcountmodel1, type="response",
                     newdata = data.frame(rotation=levels(pdsroc$rotation)),se.fit=T)
preddat$ci <- preddat$se.fit*qnorm(0.975)


#tiff("supplement_nodcountbygroup.tiff",width=800, height=510,res=100)
plot(factor(pdsroc$rotation), pdsroc$nodCount,ylim=c(0,80),lwd=1.5,
     border=grey(0.5),xaxt="n",cex.axis=1.3, cex.lab=1.4,ylab="nodules per plant",
     xlab="years with a host crop | years since last host crop")
text(1:6,rep(0,6), labels=c("ab","ab","b","ab","b","ac"),cex=1.3)
points(1:6, preddat$fit,cex=1.5,lwd=2,pch=17)
segments(1:6, preddat$fit-preddat$ci,1:6, preddat$fit+preddat$ci,lwd=2)
axis(1,at=1:6,labels=c("30| ","5|0","1| "," |1","0|5"," |30"),cex.axis=1.2)
#dev.off()

preddat <- predict(nodcountmodel2, type="response",
                   newdata = data.frame(yearsSinceSoy=levels(ltarnpd$yearsSinceSoy)),se.fit=T)
preddat$ci <- preddat$se.fit*qnorm(0.975)

#tiff("supplement_nodcountbygroupLTARN.tiff",width=500, height=400,res=100)
plot(factor(ltarnpd$yearsSinceSoy), ltarnpd$nodCount,ylim=c(0,100),lwd=1.5,
     border=grey(0.5),cex.axis=1.3, cex.lab=1.4,ylab="nodules per plant",xlab="years since last host crop")
points(1:3, preddat$fit,cex=1.5,lwd=2,pch=17)
segments(1:3, preddat$fit-preddat$ci,1:3, preddat$fit+preddat$ci,lwd=2)
text(1,100,labels="*",cex=2)
#dev.off()

max(preddat$fit)- preddat$fit 
(max(preddat$fit)- preddat$fit) /preddat$fit


##no difference in plant growth among rotation treatments:
shootmodel <- lm(shootMassG ~ rotation, data=pdsroc)
summary(shootmodel)
preddat <- predict(shootmodel,type="response",newdata=data.frame(rotation = levels(pdsroc$rotation)),se.fit=T)
preddat$ci <- preddat$se.fit*qnorm(0.975)
preddat$fit - max(preddat$fit)
preddat$fit - min(preddat$fit)
shootmodel2 <- lm(shootMassG ~ yearsSinceSoy, data=ltarnpd)
summary(shootmodel2)
shootTukey2 <- summary(glht(shootmodel2, mcp(yearsSinceSoy="Tukey")))
shootTukey2

#tiff("supplement_shootMass.tiff",width=800, height=510,res=100)
plot(factor(pdsroc$rotation), pdsroc$shootMassG,lwd=1.5,
     border=grey(0),xaxt="n",cex.axis=1.3, cex.lab=1.4,ylab="shoot dry mass per plant (g)",
     xlab="years with a host crop | years since last host crop")
axis(1,at=1:6,labels=c("30| ","5|0","1| "," |1","0|5"," |30"),cex.axis=1.2)
#dev.off()

preddat <- predict(shootmodel2,type="response",newdata=data.frame(yearsSinceSoy = levels(ltarnpd$yearsSinceSoy)),se.fit=T)
preddat$ci <- preddat$se.fit*qnorm(0.975)

#tiff("supplement_shootMassLTARN.tiff",width=500, height=400,res=100)
plot(factor(ltarnpd$yearsSinceSoy), ltarnpd$shootMassG,lwd=1.5,
     border=grey(0.5),cex.axis=1.3, cex.lab=1.4,ylab="shoot dry mass per plant (g)",xlab="years since last host crop",ylim=c(0,4))
points(1:3, preddat$fit,cex=1.5,lwd=2,pch=17)
segments(1:3, preddat$fit-preddat$ci,1:3, preddat$fit+preddat$ci,lwd=2)
mtext(1,at=1:3,line=-1,text=c("ab","bc","c"),cex=1.2)
#dev.off()


##power analysis: what difference in means would be needed to detect between group variation in shoot mass?

meanies <- tapply(pdsroc$shootMassG,pdsroc$rotation,mean)
meaniesltarn <- tapply(ltarnpd$shootMassG, ltarnpd$yearsSinceSoy, mean)
withingroupvar <- tapply(pdsroc$shootMassG,pdsroc$rotation,var)
withingroupvarltarn <- tapply(ltarnpd$shootMassG,ltarnpd$yearsSinceSoy,var)
withingroupsd <- tapply(pdsroc$shootMassG,pdsroc$rotation,sd)
betweengroupvar <- var(meanies)
power.anova.test(groups=6, n=16,within.var=mean(withingroupvar),power=0.8)# SCN (all 6 groups)
power.anova.test(groups=6, n=16,within.var=mean(withingroupvar),between.var = betweengroupvar)# SCN (all 6 groups)

power.anova.test(groups=3, n=16,within.var=mean(withingroupvarltarn),between.var = var(meaniesltarn)) # LTARN (only 3 groups)
##power: 0.62

###easier to communicate in t-test. 
power.t.test(n=16,sd=mean(tapply(pdsroc$shootMassG,pdsroc$rotation,sd)),power=0.8)
power.t.test(n=16,sd=mean(tapply(pdsroc$shootMassG,pdsroc$rotation,sd)),delta = 0.3)
power.t.test(n=16,sd=mean(tapply(ltarnpd$shootMassG,ltarnpd$yearsSinceSoy,sd)),power=0.9)
power.t.test(n=16,sd=mean(tapply(ltarnpd$shootMassG,ltarnpd$yearsSinceSoy,sd)),delta = 0.6)
power.t.test(n=16,sd=mean(tapply(ltarnpd$shootMassG,ltarnpd$yearsSinceSoy,sd)),delta = 0.5)

###t-test: with 16 plants and this variance, you'd need a difference of 0.57g to detect a significant difference
# in mean shoot mass. the biggest idfferences in the pdsroc was 0.29
max(meanies) - min(meanies)
max(meaniesltarn) - min(meaniesltarn)
