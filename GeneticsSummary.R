## genetics data summary
## summaries and graphics for RI data
library(scales)

dat<-read.csv("HostRIDat.csv",header=TRUE)

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


## correlation GBS vs mtDNA (not including those filled in)

obs<-which(dat$GBS_ref!="0.15897 + nuc * 19.41020")

pdf("F3_geneticData.pdf",width=10,height=5)
par(mfrow=c(1,2))
par(mar=c(4.5,5,2.5,1.5))
o<-lm(dat$GBS[obs] ~ dat$nucDNA[obs])
plot(dat$nucDNA[obs],dat$GBS[obs],pch=19,xlab="Nuclear distance",ylab="GBS distance",cex.lab=1.4,cex.axis=1.1)
abline(o$coefficients,lwd=1.5)
title(main="(A) Nuclear vs GBS",cex.main=1.3)
text(.03,.15,expression(paste(r[2]," = 0.68",sep="")),cex=1.25)

cs<-c("deepskyblue4","tomato4")

plot(sort(dat$GBS),xlab="Rank",ylab=expression(F[ST]),pch=19,col=cs[spd+1][order(dat$GBS)],cex.lab=1.4,cex.axis=1.1)
legend(1,.88,c("Within species","Between species"),pch=19,col=cs,bty='n')
title(main="(B) Genetic differentiation continuum",cex.main=1.3)
dev.off()
