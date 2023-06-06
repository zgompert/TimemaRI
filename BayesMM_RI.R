## bayesian mixed model for effects of Fst on Habitat RI and Sexual RI
## with effect of same vs. different species

dat<-read.csv("HostRIDat.csv",header=TRUE)
cdat<-read.csv("CountData.csv",header=TRUE)

stand<-function(x){
	x<-(x-mean(x))/sd(x)
	return(x)
}

## get species/pop numbers
spN<-as.numeric(as.factor(c(dat$Species1,dat$Species2)))
pn1<-spN[1:42];pn2<-spN[43:84]

## host distance
hdist<-abs(dat$PctPick1-dat$PctPick2)/100

## sexual isolation from Coyne & Orr, not using
## IPSI anymore
p11<-cdat$yF1M1/cdat$nF1M1
p12<-cdat$yF1M2/cdat$nF1M2
p21<-cdat$yF2M1/cdat$nF2M1
p22<-cdat$yF2M2/cdat$nF2M2
si<-1-(p12+p21)/(p11+p22)


## same = 0 vs different = 1species
spd<-as.numeric(pn1!=pn2)

#### preliminary analysis, shows no detectable effect of geography ###
o<-lm(hdist ~ dat$Geography)
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.20050    0.04142   4.841 1.96e-05 ***
#dat$GeographyP  0.06132    0.05722   1.072     0.29    
#Residual standard error: 0.1852 on 40 degrees of freedom
#Multiple R-squared:  0.0279,	Adjusted R-squared:  0.003602 
#F-statistic: 1.148 on 1 and 40 DF,  p-value: 0.2903

o<-lm(si~ dat$Geography)
#               Estimate Std. Error t value Pr(>|t|)
#(Intercept)      0.3152     0.1919   1.642    0.112
#dat$GeographyP  -0.4063     0.2669  -1.522    0.140
#Residual standard error: 0.7182 on 27 degrees of freedom
#  (13 observations deleted due to missingness)
#Multiple R-squared:  0.07904,	Adjusted R-squared:  0.04493 
#F-statistic: 2.317 on 1 and 27 DF,  p-value: 0.1396

cor.test(hdist,si)


library(rjags)

## R function to fit a Bayesian linear model with population effect for distance matrixes
## similar to clarke et al. 2002

## variable definitions
## Y = genetic distance matrix (dist)
## X = predictor distance matrixes (list of dist)
## n = number of populations or matrix rows (integer)
## p1 = numeric ID for species 1
## p2 = numeric ID for species 2
## nmcmc = number of mcmc iterations (integer)
## burnin = number of initial mcmc iterations to discard (integer)

## mixed model regression MCMC function
mmrMcmc<-function(Y=NA,X1=NULL,X2=NULL,X3=NULL,npop=NA,p1=NA,p2=NA,nmcmc=5000,burnin=1000,thin=5,nchain=2){

    library(rjags)
    load.module("dic")

    ## specifies the linear model for jags, single path or differences
    jagsModel="
        model {

            ## normal likelihood and linear model: mu = beta0 + beta X + alpha_i + alpha_j
            for (i in 1:n) {
                Y[i] ~ dnorm(mu[i],taue)
                mu[i]<-beta0 + inprod(X[i,],beta) + alpha[pop1[i]] + alpha[pop2[i]]
            }
            
            ## covariate priors
            for (k in 1:p) {
                beta[k] ~ dnorm(0,0.001)    
            }
            ## intercept prior
            beta0 ~ dnorm(0,0.001)
            ## population effect priors
            for (j in 1:npop){
                alpha[j] ~ dnorm(0, taua)
            }    
            
            ## residual priors
            taua ~ dgamma(1,0.01)
            taue ~ dgamma(1,0.01)
        }  
    "


    ## format data
    if (is.null(X2)){
        X<-as.matrix(X1)
        p<-1
    }
    else if (is.null(X3)){
        X<-as.matrix(cbind(X1,X2))
        p<-2
    }
    else{
        X<-as.matrix(cbind(X1,X2,X3))
        p<-3
    }        
    n<-length(Y)
    data<-list(Y=Y,X=X,n=n,p=p,npop=npop,pop1=p1,pop2=p2)
  
    ## model initiatlization
    model<-jags.model(textConnection(jagsModel),data=data,n.chains=nchain,n.adapt=1000) ## add back init
    update(model,n.iter=burnin)

    ## mcmc estimation of model parameters
    mcmcSamples<-coda.samples(model=model,variable.names=c("beta0","beta","taue","taua"), n.iter=nmcmc,thin=thin)
    dicSams<-dic.samples(model=model,n.iter=nmcmc,thin=thin,type="pD")
    return(list(param=mcmcSamples,dic=dicSams))
  
}
## mixed model regression MCMC function, null model version
mmrMcmcNull<-function(Y=NA,npop=NA,p1=NA,p2=NA,nmcmc=5000,burnin=1000,thin=5,nchain=2){

    library(rjags)
    load.module("dic")

    ## specifies the linear model for jags, single path or differences
    jagsModel="
        model {

            ## normal likelihood and linear model: mu = beta0 + beta X + alpha_i + alpha_j
            for (i in 1:n) {
                Y[i] ~ dnorm(mu[i],taue)
                mu[i]<-beta0  + alpha[pop1[i]] + alpha[pop2[i]]
            }
            
            ## intercept prior
            beta0 ~ dnorm(0,0.001)
            ## population effect priors
            for (j in 1:npop){
                alpha[j] ~ dnorm(0, taua)
            }    
            
            ## residual priors
            taua ~ dgamma(1,0.01)
            taue ~ dgamma(1,0.01)
        }  
    "


    n<-length(Y)
    data<-list(Y=Y,n=n,npop=npop,pop1=p1,pop2=p2)
  
    ## model initiatlization
    model<-jags.model(textConnection(jagsModel),data=data,n.chains=nchain,n.adapt=1000) ## add back init
    update(model,n.iter=burnin)

    ## mcmc estimation of model parameters
    mcmcSamples<-coda.samples(model=model,variable.names=c("beta0","taue","taua"), n.iter=nmcmc,thin=thin)
    dicSams<-dic.samples(model=model,n.iter=nmcmc,thin=thin,type="pD")
    return(list(param=mcmcSamples,dic=dicSams))
  
}

#################### GBS #################################

## habitat isolation

comp<-which(is.na(dat$GBS)==FALSE)
sub_dat<-dat[comp,]

outHostFull<-mmrMcmc(Y=hdist[comp],X1=stand(sub_dat$GBS),X2=spd[comp],X3=stand(sub_dat$GBS)*spd[comp],npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outHostNoInt<-mmrMcmc(Y=hdist[comp],X1=stand(sub_dat$GBS),X2=spd[comp],X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outHostGen<-mmrMcmc(Y=hdist[comp],X1=stand(sub_dat$GBS),X2=NULL,X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outHostSp<-mmrMcmc(Y=hdist[comp],X1=spd[comp],X2=NULL,X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outHostNull<-mmrMcmcNull(Y=hdist[comp],npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)

oh<-rbind(outHostFull[[1]][[1]],outHostFull[[1]][[2]],outHostFull[[1]][[3]])
bnds<-c(min(stand(sub_dat$GBS)),max(stand(sub_dat$GBS)))
sx<-seq(-1.5,1.7,.1)
Nx<-length(sx)
Lw<-matrix(NA,nrow=dim(oh)[1],ncol=Nx)
Lb<-matrix(NA,nrow=dim(oh)[1],ncol=Nx)
for(i in 1:dim(oh)[1]){
	Lw[i,]<-oh[i,4] + oh[i,2] * sx
	Lb[i,]<-oh[i,4] + oh[i,1] + (oh[i,2] +oh[i,3]) * sx
}
qLw<-apply(Lw,2,quantile,probs=c(.5,.05,.95))
qLb<-apply(Lb,2,quantile,probs=c(.5,.05,.95))
cs<-c("deepskyblue4","tomato4")

pdf("BayGBSxHI.pdf",width=5,height=5)
par(mar=c(4.5,4.5,.5,.5))
plot(stand(sub_dat$GBS),hdist[comp],type='n',xlab="Genetic distance",ylab="Habitat isolation",cex.lab=1.4)
polygon(c(sx,rev(sx)),c(qLw[2,],rev(qLw[3,])),border=NA,col=alpha(cs[1],.3))
polygon(c(sx,rev(sx)),c(qLb[2,],rev(qLb[3,])),border=NA,col=alpha(cs[2],.3))
points(stand(sub_dat$GBS),hdist[comp],col=cs[spd[comp]+1],pch=19)
lines(sx,qLw[1,],col=cs[1],lwd=1.5)
lines(sx,qLb[1,],col=cs[2],lwd=1.5)
dev.off()


## sexual isolation

comp<-which(is.na(dat$GBS+si)==FALSE)
sub_dat<-dat[comp,]

outSexFull<-mmrMcmc(Y=si[comp],X1=stand(sub_dat$GBS),X2=spd[comp],X3=stand(sub_dat$GBS)*spd[comp],npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexNoInt<-mmrMcmc(Y=si[comp],X1=stand(sub_dat$GBS),X2=spd[comp],X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexGen<-mmrMcmc(Y=si[comp],X1=stand(sub_dat$GBS),X2=NULL,X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexSp<-mmrMcmc(Y=si[comp],X1=spd[comp],X2=NULL,X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexNull<-mmrMcmcNull(Y=si[comp],npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)

oh<-rbind(outSexGen[[1]][[1]],outSexGen[[1]][[2]],outSexGen[[1]][[3]])
bnds<-c(min(stand(sub_dat$GBS)),max(stand(sub_dat$GBS)))
sx<-seq(-1.5,1.7,.1)
Nx<-length(sx)
Lw<-matrix(NA,nrow=dim(oh)[1],ncol=Nx)
Lb<-matrix(NA,nrow=dim(oh)[1],ncol=Nx)
for(i in 1:dim(oh)[1]){
        Lw[i,]<-oh[i,2] + oh[i,1] * sx
}
qLw<-apply(Lw,2,quantile,probs=c(.5,.05,.95))
cs<-c("deepskyblue4","tomato4")


pdf("BayGBSxSI.pdf",width=5,height=5)
par(mar=c(4.5,4.5,.5,.5))
plot(stand(sub_dat$GBS),si[comp],type='n',xlab="Genetic distance",ylab="Sexual isolation",cex.lab=1.4,ylim=c(-.4,1.5))
polygon(c(sx,rev(sx)),c(qLw[2,],rev(qLw[3,])),border=NA,col=alpha("darkgray",.3))
points(stand(sub_dat$GBS),si[comp],col=cs[spd[comp]+1],pch=19)
lines(sx,qLw[1,],col="black",lwd=1.5)
dev.off()


################### hist effect estimates #####

pdf("F6_GBSxRIHist.pdf",width=5,height=10)
par(mfrow=c(3,1))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.45;ca<-1.1;cm<-1.45
cs<-c("deepskyblue4","tomato4")
oh<-rbind(outHostFull[[1]][[1]],outHostFull[[1]][[2]],outHostFull[[1]][[3]])
par(mar=c(4.5,4.5,2.5,1))
hist(oh[,2],col=alpha(cs,.5)[1],border=alpha(cs)[1],xlab="Slope",ylab="Posterior probability",main="(A) Host isolation within species",cex.lab=cl,cex.axis=ca,cex.main=cm)
abline(v=0,lwd=1.5)
mean(oh[,2] > 0)
#[1] 0.3571667
hist(oh[,2] +oh[,3],col=alpha(cs,.5)[2],border=alpha(cs)[2],xlab="Slope",ylab="Posterior probability",main=" (B) Host isolation between species",cex.lab=cl,cex.axis=ca,cex.main=cm)
abline(v=0,lwd=1.5)
mean((oh[,2] +oh[,3]) > 0)
#0.9943333
oh<-rbind(outSexGen[[1]][[1]],outSexGen[[1]][[2]],outSexGen[[1]][[3]])
hist(oh[,1],col=alpha("darkgray",.5)[1],border="darkgray",xlab="Slope",ylab="Posterior probability",main="(C) Sexual isolation within and between species",cex.lab=cl,cex.axis=ca,cex.main=cm)
abline(v=0,lwd=1.5)
mean(oh[,1] > 0)
#[1] 0.9711667
dev.off()

pdf("sfig_GBSxRIHist.pdf",width=10,height=10)
layout(matrix(1:4,nrow=2,ncol=2,byrow=FALSE),widths=c(5,5),heights=c(5,5))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.45;ca<-1.1;cm<-1.3
cs<-c("deepskyblue4","tomato4")
oh<-rbind(outHostFull[[1]][[1]],outHostFull[[1]][[2]],outHostFull[[1]][[3]])
par(mar=c(4.5,4.5,2.5,1))
hist(oh[,2],col=alpha(cs,.5)[1],border=alpha(cs)[1],xlab="Slope",ylab="Posterior probability",main="(A) Host isolation within species",cex.lab=cl,cex.axis=ca,cex.main=cm)
abline(v=0,lwd=1.5)
#[1]0.3571667
hist(oh[,2] +oh[,3],col=alpha(cs,.5)[2],border=alpha(cs)[2],xlab="Slope",ylab="Posterior probability",main=" (B) Host isolation between species",cex.lab=cl,cex.axis=ca,cex.main=cm)
abline(v=0,lwd=1.5)
mean((oh[,2] +oh[,3]) > 0)
#[1]  0.9943333
oh<-rbind(outSexFull[[1]][[1]],outSexFull[[1]][[2]],outSexFull[[1]][[3]])
par(mar=c(4.5,4.5,2.5,1))
hist(oh[,2],col=alpha(cs,.5)[1],border=alpha(cs)[1],xlab="Slope",ylab="Posterior probability",main="(C) Sexual isolation within species",cex.lab=cl,cex.axis=ca,cex.main=cm)
abline(v=0,lwd=1.5)
mean((oh[,2]) > 0)
#[1]  0.3286667 
hist(oh[,2] +oh[,3],col=alpha(cs,.5)[2],border=alpha(cs)[2],xlab="Slope",ylab="Posterior probability",main=" (D) Sexual isolation between species",cex.lab=cl,cex.axis=ca,cex.main=cm)
abline(v=0,lwd=1.5)
mean((oh[,2] +oh[,3]) > 0)
#[1]  0.478
dev.off()


## trying breakpoint regression for HI and SI as a function of GBS
library(strucchange)

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
## habitat isolation

comp<-which(is.na(dat$GBS)==FALSE)
sub_dat<-dat[comp,]

o<-lm(hdist[comp] ~ stand(sub_dat$GBS))

o<-breakpoints(hdist[comp] ~ stand(sub_dat$GBS),h=.15,breaks=1)
## breakpoint at 75% percentile
summary(o) ## one breakpoint winss by BIC
#m   0        1 
#RSS   0.6940   0.4548
#BIC -29.2436 -33.7036
cos<-coef(o)
## break
bpt<-quantile(stand(sub_dat$GBS),.75)
x1<-seq(-2,bpt,.01)
x2<-seq(bpt,2,.01)
y1<-cos[1,1] + x1 * cos[1,2]
y2<-cos[2,1] + x2 * cos[2,2]
pdf("breakHabitat.pdf",width=5,height=5)
par(mar=c(5,5,.5,.5))
plot(stand(sub_dat$GBS),hdist[comp],pch=19,xlab="Genetic distance",ylab="Habitat isolation",cex.lab=1.4,cex.axis=1.1)
lines(x1,y1);lines(x2,y2)
dev.off()

## sexual isolation

comp<-which(is.na(dat$GBS+si)==FALSE)
sub_dat<-dat[comp,]

o<-lm(si[comp] ~ stand(sub_dat$GBS))
#                   Estimate Std. Error t value Pr(>|t|)  
#(Intercept)          0.0791     0.1354   0.584   0.5642  
#stand(sub_dat$GBS)   0.3084     0.1379   2.235   0.0346 *
#Residual standard error: 0.7034 on 25 degrees of freedom
#Multiple R-squared:  0.1666,	Adjusted R-squared:  0.1333 
#F-statistic: 4.997 on 1 and 25 DF,  p-value: 0.03456


o<-breakpoints(si[comp] ~ stand(sub_dat$GBS),h=.15,breaks=1)
summary(o) ## 0 breakpoints based on BIC
#m   0     1    
#RSS 12.37 11.05
#BIC 65.43 72.29


################# figures with breakpoint regression too #################

pdf("F5_GBSxRI.pdf",width=10,height=10)
layout(matrix(1:4,nrow=2,ncol=2,byrow=FALSE),widths=c(5,5),heights=c(5,5))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.45;ca<-1.1;cm<-1.3
comp<-which(is.na(dat$GBS)==FALSE)
sub_dat<-dat[comp,]
oh<-rbind(outHostFull[[1]][[1]],outHostFull[[1]][[2]],outHostFull[[1]][[3]])
bnds<-c(min(stand(sub_dat$GBS)),max(stand(sub_dat$GBS)))
sx<-seq(-1.5,2.5,.1)
Nx<-length(sx)
Lw<-matrix(NA,nrow=dim(oh)[1],ncol=Nx)
Lb<-matrix(NA,nrow=dim(oh)[1],ncol=Nx)
for(i in 1:dim(oh)[1]){
	Lw[i,]<-oh[i,4] + oh[i,2] * sx
	Lb[i,]<-oh[i,4] + oh[i,1] + (oh[i,2] +oh[i,3]) * sx
}
qLw<-apply(Lw,2,quantile,probs=c(.5,.05,.95))
qLb<-apply(Lb,2,quantile,probs=c(.5,.05,.95))
cs<-c("deepskyblue4","tomato4")

plot(stand(sub_dat$GBS),hdist[comp],type='n',xlab="Genetic distance",ylab="Habitat isolation",cex.lab=cl,cex.axis=ca)
polygon(c(sx,rev(sx)),c(qLw[2,],rev(qLw[3,])),border=NA,col=alpha(cs[1],.3))
polygon(c(sx,rev(sx)),c(qLb[2,],rev(qLb[3,])),border=NA,col=alpha(cs[2],.3))
points(stand(sub_dat$GBS),hdist[comp],col=cs[spd[comp]+1],pch=19)
lines(sx,qLw[1,],col=cs[1],lwd=1.5)
lines(sx,qLb[1,],col=cs[2],lwd=1.5)
legend(-.7,.6,c("Within species","Between species"),pch=19,col=cs,bty='n')

title(main="(A) Habitat isolation mixed model",cex.main=cm)

o<-lm(hdist[comp] ~ stand(sub_dat$GBS))

o<-breakpoints(hdist[comp] ~ stand(sub_dat$GBS),h=.15,breaks=1)
## breakpoint at 75% percentile
bpt<-quantile(stand(sub_dat$GBS),.75)
x1<-seq(-2,bpt,.01)
x2<-seq(bpt,2,.01)
y1<-cos[1,1] + x1 * cos[1,2]
y2<-cos[2,1] + x2 * cos[2,2]
sd<-stand(sub_dat$GBS)
lb<-min(sd[spd[comp]==1])
ub<-max(sd[spd[comp]==0])

plot(stand(sub_dat$GBS),hdist[comp],type='n',xlab="Genetic distance",ylab="Habitat isolation",cex.lab=cl,cex.axis=ca)
#polygon(c(lb,ub,ub,lb),c(-1,-1,1,1),col=alpha("gray",.4),border=NA)
points(stand(sub_dat$GBS),hdist[comp],col=cs[spd[comp]+1],pch=19)
lines(x1,y1);lines(x2,y2)
title(main="(B) Habitat isolation breakpoint regression",cex.main=cm)

### sexual ###

comp<-which(is.na(dat$GBS+si)==FALSE)
sub_dat<-dat[comp,]

oh<-rbind(outSexGen[[1]][[1]],outSexGen[[1]][[2]],outSexGen[[1]][[3]])
bnds<-c(min(stand(sub_dat$GBS)),max(stand(sub_dat$GBS)))
sx<-seq(-1.5,2.5,.1)
Nx<-length(sx)
Lc<-matrix(NA,nrow=dim(oh)[1],ncol=Nx)
for(i in 1:dim(oh)[1]){
	Lc[i,]<-oh[i,2] + oh[i,1]*sx 
}
qLc<-apply(Lc,2,quantile,probs=c(.5,.05,.95))
cs<-c("deepskyblue4","tomato4","darkgray")

plot(stand(sub_dat$GBS),si[comp],type='n',xlab="Genetic distance",ylab="Sexual isolation",cex.lab=cl,cex.axis=ca)
polygon(c(sx,rev(sx)),c(qLc[2,],rev(qLc[3,])),border=NA,col=alpha(cs[3],.3))
points(stand(sub_dat$GBS),si[comp],col=cs[spd[comp]+1],pch=19)
lines(sx,qLc[1,],col=cs[3],lwd=1.5)
title(main="(C) Sexual isolation mixed model",cex.main=cm)

o<-lm(si[comp] ~ stand(sub_dat$GBS))
plot(stand(sub_dat$GBS),si[comp],type='n',xlab="Genetic distance",ylab="Sexual isolation",cex.lab=cl,cex.axis=ca)
#polygon(c(lb,ub,ub,lb),c(-1,-1,1,1),col=alpha("gray",.4),border=NA)
points(stand(sub_dat$GBS),si[comp],col=cs[spd[comp]+1],pch=19)
abline(o$coefficients)
title(main="(D) Sexual isolation breakpoint regression",cex.main=cm)
dev.off()

save(list=ls(),file="bayes.rdat")
