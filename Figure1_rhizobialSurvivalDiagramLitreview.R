##Code for figure 1: combines published reports of population size of symbiotic rhizobia.

## figuring out points to enter into spreadsheet (based on graphs)
##data from Hirsch 1996 (she also presented Nutman 1969) 
## I used image-J to estimate data points for her Fig. 1A (data from Nutman et al. 1969) and 
# Fig. 3A (her own data). I decided not to present the Nutman data because it's not clear when the last  
# host crop was (measurements came from fallow soil)

phdata <- read.csv("HirschSoilSurvivalData.csv",stringsAsFactors = F)
phdata$date <- as.Date(phdata$date, format=c("%m/%d/%Y"))
# convert image J coordinates to data:
hirsch <- phdata[phdata$source=="Hirsch1996",]
nutman <- phdata[phdata$source=="Nutman1969",]
hirsch$days <- as.numeric(hirsch$date - min(hirsch$date,na.rm=TRUE)) # peas cleared october 1987
nutman$daysWithoutHost <- as.numeric(nutman$date - min(nutman$date,na.rm=TRUE))

plot(log10(rhizobiaPerGSoil) ~ Y, data=hirsch[hirsch$pointType=="yaxis",])
plot(days ~ X, data=hirsch[hirsch$pointType=="xaxis",])

yaxis <- lm(log10(rhizobiaPerGSoil) ~ Y, data=hirsch[hirsch$pointType=="yaxis",])
xaxis <- lm(days ~ X, data=hirsch[hirsch$pointType=="xaxis",])
hirsch <- hirsch[!hirsch$pointType %in% c("xaxis","yaxis"),]
hirsch$log10rhizobiaPerGSoil <- hirsch$Y*coef(yaxis)[2] + coef(yaxis)[1]
hirsch$days <- hirsch$X * coef(xaxis)[2] + coef(xaxis)[1]
hirsch$date <- as.Date("1987-08-01",format=c("%Y-%m-%d")) + hirsch$days
hirsch$daysWithoutHost <- as.numeric(hirsch$date -as.Date("1987-10-15",format=c("%Y-%m-%d") ))
##peas were planted in 1987  
hirsch <- hirsch[,c("pointType","log10rhizobiaPerGSoil","daysWithoutHost","date")]
plot(log10rhizobiaPerGSoil ~ date, data=hirsch)

yaxis <- lm(log10(rhizobiaPerGSoil) ~ Y, data=nutman[nutman$pointType=="yaxis",])
xaxis <- lm(daysWithoutHost ~ X, data=nutman[nutman$pointType=="xaxis",])
nutman <- nutman[!nutman$pointType %in% c("yaxis","xaxis"),]
nutman$log10rhizobiaPerGSoil <- nutman$Y * coef(yaxis)[2] + coef(yaxis)[1]
nutman$daysWithoutHost <- nutman$X * coef(xaxis)[2] + coef(xaxis)[1]

##Add in linear regression from Revellin et al. 1996
revellin <- data.frame(years = c(1:12))
revellin$log10rhizobiaPerGSoil <- -0.18*revellin$years + 5.3
revellin$daysWithoutHost <- revellin$years * 360

##make a csv that combines data from different studies. 
## 
#revellin <- data.frame(source="Revellin et al. 1996", method="MPN",host="soybean",yearsSinceHost = revellin$years,
#                       log10rhizobiaPerGSoil = revellin$log10rhizobiaPerGSoil)
nutmanraw <- nutman
nutman <- data.frame(source="Nutman et al. 1969", method="MPN", host = nutmanraw$pointType,
                        yearsSinceHost = nutmanraw$daysWithoutHost/365, log10rhizobiaPerGSoil=nutman$log10rhizobiaPerGSoil)
nutman$host <- as.character(nutman$host)
nutman[nutman$host=="R.leg.trifolii","host"] <- "clover"
nutman[nutman$host=="R.leg.vicae","host"] <- "vetch"
nutman[nutman$host=="clover","pch"] <- 9
nutman[nutman$host=="vetch","pch"] <- 7
hirschraw <- hirsch
hirsch <- data.frame(source="Hirsch et al. 1996",method="MPN",host="peas",group = hirschraw$pointType,yearsSinceHost = hirschraw$daysWithoutHost/365,
                     log10rhizobiaPerGSoil = hirschraw$log10rhizobiaPerGSoil,pch=4)
## for hirsch and nutman data, just take a subset to plot (they have lots of samples over the first
# 2 years)
hirschtoplot <- hirsch[hirsch$group=="peasin1987",]
plot(hirschtoplot$yearsSinceHost,hirschtoplot$log10rhizobiaPerGSoil,type="b")
## basically, populaiton oscillates seasonally. I'll just take the first point at the high end
# of the oscillation. 
hirschtoplot <- hirschtoplot[c(2,20,21,22),]
plot(hirschtoplot$yearsSinceHost, hirschtoplot$log10rhizobiaPerGSoil, type="b")
## that works. Skips the oscillation captured by more frequent sampling while capturing long-term trends
nutman1 <- nutman[nutman$host=="clover",]
nutman2 <- nutman[nutman$host=="vetch",]
plot(nutman1$yearsSinceHost, nutman1$log10rhizobiaPerGSoil, type="b")
points(nutman2$yearsSinceHost, nutman2$log10rhizobiaPerGSoil, type="b")
points(log10rhizobiaPerGSoil ~ yearsSinceHost, data=hirsch[hirsch$group=="peasin1987",],type="b")
## I think that one's OK to present without thinning.
##but looking closer at the Hirsch paper (where I got the nutman data), there isn't any information
# on when the last host crop was. So I shouldn't include this dataset. 
obaton <- data.frame(source="Obaton et al. 2002", method="MPN",host="soybean", yearsSinceHost = 
                       c(1,4,12,2,6,16,1,20),site=c("t","t","t","m","m","m","d","d") ,log10rhizobiaPerGSoil = 
                       c(log10(c(5*10^4,10^3, 10,2*10^4, 7*10^3, 10, 3.2*10^4, 0.36*10^4))),pch=5)
##they report population size for three different fields observed over different time points. 
## I'm just treating them as replicates for now. 
narozna <- data.frame(source="Narozna et al. 1996", method="MPN", host="soybean", yearsSinceHost=15,
                      log10rhizobiaPerGSoil=3,pch=19)

## add in some data from recent crop residues (hiltbold et al. 1985): data from Fig. 1
hiltbold <- data.frame(source="Hiltbold et al. 1985", method="MPN", host="soybean", 
                       yearsSinceHost =c(0,1,2),log10rhizobiaPerGSoil=c(7,5,5) ,pch=24)


## for revellin data, show points corrsponding to the slope and intercept of the regression
revslopes <- c(-0.18, -0.49)
revintercepts <- c(5.30, 4.73)
revyears <- list(gr1=c(1:11), gr2=c(1:7))
##for the Revellin data, show their linear regresions for two soil types
ptt <- c(24,22,4,21,9)
bgc <- grey(0,alpha=0.5)
lc <- grey(0,alpha=0.5)
ptcex=2
ylabtext <- expression("rhizobia g"^-1*"soil")
 
#tiff("Figure1_rhizobiaSurvivalLitreview.tiff", width=6, height=5, units="in",res=100)
par(mar=c(5,5,2,2)+0.1)
plot(c(0:20),seq(1,7,length.out=21),type="n",
     ylab="",yaxt="n", xlab="years since last host crop or inoculation",cex.lab=1.4,cex.axis=1.3)
mtext(ylabtext, side=2,line=2.5,cex=1.4)
axis(2, at=c(1:7),labels=parse(text=paste0("10","^",1:7)),las=1,cex.axis=1.3)
abline(h=2,lty=2,lwd=2)
points(log10rhizobiaPerGSoil~yearsSinceHost, data=hirschtoplot,pch=ptt[1],type="b",bg=bgc, col=lc,cex=ptcex,lwd=2)
points(hiltbold$yearsSinceHost,hiltbold$log10rhizobiaPerGSoil,pch=ptt[2], type="b",bg=bgc, col=lc,cex=ptcex,lwd=2)
points(revyears$gr1, revslopes[1]*revyears$gr1 + revintercepts[1],pch=ptt[3],bg=bgc, col=lc,type="b",cex=ptcex,lwd=2)
points(revyears$gr2, revslopes[2]*revyears$gr2 + revintercepts[2],pch=ptt[3],bg=bgc, col=lc,type="b",cex=ptcex,lwd=2)
points(log10rhizobiaPerGSoil ~ yearsSinceHost, data=obaton[obaton$site=="t",],pch=ptt[4],type="b",bg=bgc,col=lc,cex=ptcex,lwd=2)
points(log10rhizobiaPerGSoil ~ yearsSinceHost, data=obaton[obaton$site=="m",],pch=ptt[4],type="b",bg=bgc,col=lc,cex=ptcex,lwd=2)
points(log10rhizobiaPerGSoil ~ yearsSinceHost, data=obaton[obaton$site=="d",],pch=ptt[4],type="b",bg=bgc,col=lc,cex=ptcex,lwd=2)
points(log10rhizobiaPerGSoil ~ yearsSinceHost, data=narozna,pch=ptt[5],col=lc,bg=bgc,cex=ptcex,lwd=2)
#points(log10rhizobiaPerGSoil ~ yearsSinceHost, data=nutman1, type="b")
#points(log10rhizobiaPerGSoil ~ yearsSinceHost, data=nutman2, type="b")
legend("topright",legend=c(1:5),pt.cex=2,pch=ptt,bty="n",pt.lwd=2,title = "ref.",col=lc, cex=1.3,pt.bg=bgc)
par(mar=c(5,4,4,2)+0.1)
#dev.off()

