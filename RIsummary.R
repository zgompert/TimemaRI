## summaries and graphics for RI data
library(scales)

dat<-read.csv("HostRIDat.csv",header=TRUE)
cdat<-read.csv("CountData.csv",header=TRUE)

p11<-cdat$yF1M1/cdat$nF1M1
p12<-cdat$yF1M2/cdat$nF1M2
p21<-cdat$yF2M1/cdat$nF2M1
p22<-cdat$yF2M2/cdat$nF2M2
si<-1-(p12+p21)/(p11+p22)

stand<-function(x){
        x<-(x-mean(x))/sd(x)
        return(x)
}

## get species/pop numbers
spN<-as.numeric(as.factor(c(dat$Species1,dat$Species2)))
pn1<-spN[1:42];pn2<-spN[43:84]

## host distance
hdist<-abs(dat$PctPick1-dat$PctPick2)/100

## same = 0 vs different = 1species
spd<-as.numeric(pn1!=pn2)


pdf("F3_TimemaRI.pdf",width=5,height=10)
par(mfrow=c(3,1))
par(mar=c(3,5.5,2.5,1))
boxplot(cbind(hdist,si),col="white",names=c("Host","Sexual"),ylab="Reproductive isolation",pch=NA,cex.lab=1.4,cex.axis=1.2,ylim=c(-2,1))
points(jitter(rep(1,length(hdist)),factor=2.5),hdist,pch=19,col=alpha("forestgreen",.5))
points(jitter(rep(2,length(si)),factor=2.5),si,pch=19,col=alpha("magenta4",.5))
text(.6,1,paste("N = ",length(hdist),sep=""))
text(1.6,1,paste("N = ",sum(is.na(si)==FALSE),sep=""))
title(main="(A) All pairs",cex.main=1.4)

boxplot(cbind(hdist[spd==0],si[spd==0]),col="white",names=c("Host","Sexual"),ylab="Reproductive isolation",pch=NA,cex.lab=1.4,cex.axis=1.2,ylim=c(-2,1))
points(jitter(rep(1,length(hdist[spd==0])),factor=2.5),hdist[spd==0],pch=19,col=alpha("forestgreen",.5))
points(jitter(rep(2,length(si[spd==0])),factor=2.5),si[spd==0],pch=19,col=alpha("magenta4",.5))
text(.6,1,paste("N = ",length(hdist[spd==0]),sep=""))
text(1.6,1,paste("N = ",sum(is.na(si)==FALSE & spd==0),sep=""))
title(main="(B) Within species pairs",cex.main=1.4)

boxplot(cbind(hdist[spd==1],si[spd==1]),col="white",names=c("Host","Sexual"),ylab="Reproductive isolation",pch=NA,cex.lab=1.4,cex.axis=1.2,ylim=c(-2,1))
points(jitter(rep(1,length(hdist[spd==1])),factor=2.5),hdist[spd==1],pch=19,col=alpha("forestgreen",.5))
points(jitter(rep(2,length(si[spd==1])),factor=2.5),dat$IPSI[spd==1],pch=19,col=alpha("magenta4",.5))
text(.6,1,paste("N = ",length(hdist[spd==1]),sep=""))
text(1.6,1,paste("N = ",sum(is.na(si)==FALSE & spd==1),sep=""))
title(main="(C) Between species pairs",cex.main=1.4)
dev.off()

