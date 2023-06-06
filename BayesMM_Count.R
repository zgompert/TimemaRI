## bayesian mixed model for effects of Fst on  Sexual RI
## with effect of same vs. different species
## here accounting for uncertainty in SI
library(rjags)

dat<-read.csv("HostRIDat.csv",header=TRUE)
cdat<-read.csv("CountData.csv",header=TRUE)

stand<-function(x){
	x<-(x-mean(x))/sd(x)
	return(x)
}

## get species/pop numbers
spN<-as.numeric(as.factor(c(dat$Species1,dat$Species2)))
pn1<-spN[1:42];pn2<-spN[43:84]


## same = 0 vs different = 1species
spd<-as.numeric(pn1!=pn2)


## R function to fit a Bayesian linear model with population effect for distance matrixes
## similar to clarke et al. 2002

## variable definitions
## Y = response variable
## taum = precision related to measurement s.e.
## X = predictor distance matrixes (list of dist)
## n = number of populations or matrix rows (integer)
## p1 = numeric ID for species 1
## p2 = numeric ID for species 2
## nmcmc = number of mcmc iterations (integer)
## burnin = number of initial mcmc iterations to discard (integer)

## mixed model regression MCMC function
mmrMcmc<-function(Y=NA,taum=NA,X1=NULL,X2=NULL,X3=NULL,npop=NA,p1=NA,p2=NA,nmcmc=5000,burnin=1000,thin=5,nchain=2){

    library(rjags)
    load.module("dic")

    ## specifies the linear model for jags, single path or differences
    jagsModel="
        model {

            ## normal likelihood and linear model: mu = beta0 + beta X + alpha_i + alpha_j
            for (i in 1:n) {
		Y[i] ~ dnorm(V[i],taum[i])
                V[i] ~ dnorm(mu[i],taue)
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
    data<-list(Y=Y,taum=taum,X=X,n=n,p=p,npop=npop,pop1=p1,pop2=p2)
  
    ## model initiatlization
    model<-jags.model(textConnection(jagsModel),data=data,n.chains=nchain,n.adapt=1000) ## add back init
    update(model,n.iter=burnin)

    ## mcmc estimation of model parameters
    mcmcSamples<-coda.samples(model=model,variable.names=c("beta0","beta","taue","taua"), n.iter=nmcmc,thin=thin)
    dicSams<-dic.samples(model=model,n.iter=nmcmc,thin=thin,type="pD")
    return(list(param=mcmcSamples,dic=dicSams))
  
}
## mixed model regression MCMC function, null model version
mmrMcmcNull<-function(Y=NA,taum=NA,npop=NA,p1=NA,p2=NA,nmcmc=5000,burnin=1000,thin=5,nchain=2){

    library(rjags)
    load.module("dic")

    ## specifies the linear model for jags, single path or differences
    jagsModel="
        model {

            ## normal likelihood and linear model: mu = beta0 + beta X + alpha_i + alpha_j
            for (i in 1:n) {
		Y[i] ~ dnorm(V[i],taum[i])
                V[i] ~ dnorm(mu[i],taue)
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
    data<-list(Y=Y,taum=taum,n=n,npop=npop,pop1=p1,pop2=p2)
  
    ## model initiatlization
    model<-jags.model(textConnection(jagsModel),data=data,n.chains=nchain,n.adapt=1000) ## add back init
    update(model,n.iter=burnin)

    ## mcmc estimation of model parameters
    mcmcSamples<-coda.samples(model=model,variable.names=c("beta0","taue","taua"), n.iter=nmcmc,thin=thin)
    dicSams<-dic.samples(model=model,n.iter=nmcmc,thin=thin,type="pD")
    return(list(param=mcmcSamples,dic=dicSams))
  
}

#################### GBS #################################
## sexual isolation


se<-rep(NA,length(p11))
for(i in 1:length(se)){
	p11<-rbeta(1000,shape1=cdat$yF1M1[i]+1,shape2=cdat$nF1M1[i]-cdat$yF1M1[i]+1)
	p12<-rbeta(1000,shape1=cdat$yF1M2[i]+1,shape2=cdat$nF1M2[i]-cdat$yF1M2[i]+1)
	p21<-rbeta(1000,shape1=cdat$yF2M1[i]+1,shape2=cdat$nF2M1[i]-cdat$yF2M1[i]+1)
	p22<-rbeta(1000,shape1=cdat$yF2M2[i]+1,shape2=cdat$nF2M2[i]-cdat$yF2M2[i]+1)
	si<-1-(p12+p21)/(p11+p22)
	se[i]<-sd(si)
}
	
	
p11<-cdat$yF1M1/cdat$nF1M1
p12<-cdat$yF1M2/cdat$nF1M2
p21<-cdat$yF2M1/cdat$nF2M1
p22<-cdat$yF2M2/cdat$nF2M2

## from Coyne and Orr
si<-1-(p12+p21)/(p11+p22)

comp<-which(is.na(dat$GBS+si)==FALSE)
sub_dat<-dat[comp,]

outSexFull<-mmrMcmc(Y=si[comp],taum=1/(se[comp]^2),X1=stand(sub_dat$GBS),X2=spd[comp],X3=stand(sub_dat$GBS)*spd[comp],npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexNoInt<-mmrMcmc(Y=si[comp],taum=1/(se[comp]^2),X1=stand(sub_dat$GBS),X2=spd[comp],X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexGen<-mmrMcmc(Y=si[comp],taum=1/(se[comp]^2),X1=stand(sub_dat$GBS),X2=NULL,X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexSp<-mmrMcmc(Y=si[comp],taum=1/(se[comp]^2),X1=spd[comp],X2=NULL,X3=NULL,npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)
outSexNull<-mmrMcmcNull(Y=si[comp],taum=1/(se[comp]^2),npop=max(c(pn1,pn2)),p1=pn1,p2=pn2,nmcmc=10000,burnin=2000,thin=5,nchain=3)

save(list=ls(),file="bayesUncert.rdat")


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

