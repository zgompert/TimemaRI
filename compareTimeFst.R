div<-read.csv("TimemaDivTimes.csv")
sp<-as.numeric(div$Species1==div$Species1.1)
cs<-c("firebrick","black")
pdf("FstVsDivTime.pdf",width=4.5,height=4.5)
par(mar=c(5,5,1,1))
plot(div$Div.Time,div$GBS.FST,pch=19,col=cs[sp+1],xlab="Divergence time (MYA)",ylab=expression(paste(F[ST])),cex.lab=1.3)
o<-lm(div$GBS.FST ~ div$Div.Time)
abline(o$coefficients)
dev.off()

cor.test(div$GBS.FST,div$Div.Time)

#	Pearson's product-moment correlation
#
#data:  div$GBS.FST and div$Div.Time
#t = 4.548, df = 12, p-value = 0.0006685
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.4584367 0.9325015
#sample estimates:
#      cor 
#0.7955169 

o<-lm(div$GBS.FST ~ div$Div.Time)
summary(o)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   -1.980      2.229  -0.889 0.391702    
#div$GBS.FST   14.833      3.261   4.548 0.000669 ***
#Residual standard error: 2.779 on 12 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.6328,	Adjusted R-squared:  0.6023 
#F-statistic: 20.68 on 1 and 12 DF,  p-value: 0.0006685

tcen<-(div$Div.Time-mean(div$Div.Time,na.rm=TRUE))/sd(div$Div.Time,na.rm=TRUE)
o<-lm(div$GBS.FST ~ tcen)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.64429    0.03983  16.176 1.64e-09 ***
#tcen         0.18799    0.04133   4.548 0.000669 ***
#Residual standard error: 0.149 on 12 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.6328,	Adjusted R-squared:  0.6023 
#F-statistic: 20.68 on 1 and 12 DF,  p-value: 0.0006685

tcen2<-tcen^2
o<-lm(div$GBS.FST ~ tcen + tcen2)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.70294    0.03868  18.175 1.49e-09 ***
#tcen         0.23001    0.03674   6.260 6.17e-05 ***
#tcen2       -0.06316    0.02317  -2.726   0.0197 *
#Residual standard error: 0.1203 on 11 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.7809,	Adjusted R-squared:  0.741
#F-statistic:  19.6 on 2 and 11 DF,  p-value: 0.0002365

o<-glm(div$GBS.FST ~ tcen)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.64429    0.03983  16.176 1.64e-09 ***
#tcen         0.18799    0.04133   4.548 0.000669 ***
#    Null deviance: 0.72594  on 13  degrees of freedom
#Residual deviance: 0.26653  on 12  degrees of freedom
#  (5 observations deleted due to missingness)
#AIC: -9.7282

o<-glm(div$GBS.FST ~ tcen + tcen2)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.70294    0.03868  18.175 1.49e-09 ***
#tcen         0.23001    0.03674   6.260 6.17e-05 ***
#tcen2       -0.06316    0.02317  -2.726   0.0197 *
#    Null deviance: 0.72594  on 13  degrees of freedom
#Residual deviance: 0.15907  on 11  degrees of freedom
#  (5 observations deleted due to missingness)
#AIC: -14.954

##### divergence time with IM
demog<-read.table("demog_times.txt",sep=",",header=TRUE)
tauc<-(demog$T_im-mean(demog$T_im))/sd(demog$T_im)
tauc2<-tauc^2
o<-lm(demog$GBS_FST ~ tauc)
summary(o)
o<-lm(demog$GBS_FST ~ tauc+tauc2)
summary(o) ## 2nd works, but the result is silliness, drop outlier point 2
tauc<-(demog$T_im[-2]-mean(demog$T_im[-2]))/sd(demog$T_im[-2])
tauc2<-tauc^2
o<-lm(demog$GBS_FST[-2] ~ tauc)
summary(o)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.35556    0.02953  12.042 6.21e-06 ***
#tauc         0.24822    0.03132   7.926 9.68e-05 ***
#Residual standard error: 0.08858 on 7 degrees of freedom
#Multiple R-squared:  0.8997,	Adjusted R-squared:  0.8854 
#F-statistic: 62.81 on 1 and 7 DF,  p-value: 9.675e-05

o<-lm(demog$GBS_FST[-2] ~ tauc+tauc2)
summary(o) 

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.25747    0.05781   4.454  0.00431 ** 
#tauc         0.24960    0.02681   9.310  8.7e-05 ***
#tauc2        0.11034    0.05849   1.886  0.10818    
#Residual standard error: 0.0758 on 6 degrees of freedom
#Multiple R-squared:  0.9371,	Adjusted R-squared:  0.9161 
#F-statistic: 44.67 on 2 and 6 DF,  p-value: 0.0002493

o<-glm(demog$GBS_FST[-2] ~ tauc)
summary(o)
#    Null deviance: 0.547822  on 8  degrees of freedom
#Residual deviance: 0.054928  on 7  degrees of freedom
#AIC: -14.35

o<-glm(demog$GBS_FST[-2] ~ tauc+tauc2)
summary(o) 
#    Null deviance: 0.547822  on 8  degrees of freedom
#Residual deviance: 0.034478  on 6  degrees of freedom
#AIC: -16.541

##### divergence time with IM/theta
tmu<-demog$T_im/demog$Theta_im
tauc<-(tmu[-2]-mean(tmu[-2]))/sd(tmu[-2])
tauc2<-tauc^2
o<-lm(demog$GBS_FST[-2] ~ tauc)
summary(o)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.35556    0.06601   5.387  0.00102 **
#tauc         0.18485    0.07001   2.640  0.03341 * 
#Residual standard error: 0.198 on 7 degrees of freedom
#Multiple R-squared:  0.499,	Adjusted R-squared:  0.4274 
#F-statistic: 6.971 on 1 and 7 DF,  p-value: 0.03341


o<-lm(demog$GBS_FST[-2] ~ tauc+tauc2)
summary(o) 
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.65124    0.11421   5.702  0.00126 **
#tauc         0.29458    0.06278   4.693  0.00335 **
#tauc2       -0.33264    0.11731  -2.836  0.02974 * 
#Residual standard error: 0.1398 on 6 degrees of freedom
#Multiple R-squared:  0.7859,	Adjusted R-squared:  0.7145 
#F-statistic: 11.01 on 2 and 6 DF,  p-value: 0.009813


o<-glm(demog$GBS_FST[-2] ~ tauc)
summary(o)
#    Null deviance: 0.54782  on 8  degrees of freedom
#Residual deviance: 0.27447  on 7  degrees of freedom
#AIC: 0.12967

o<-glm(demog$GBS_FST[-2] ~ tauc+tauc2)
summary(o) 
#    Null deviance: 0.54782  on 8  degrees of freedom
#Residual deviance: 0.11729  on 6  degrees of freedom
#AIC: -5.5224



pdf("fig_DivTimeFst.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4.5,5,2.5,2))
wsp<-as.numeric(div$Species1!=div$Species1.1)+1
cs<-c("deepskyblue4","tomato4")
plot(div$Div.Time,div$GBS.FST,col=cs[wsp],pch=19,xlab="Divergence time (MYA)",ylab=expression(F[ST]),ylim=c(0,1),xlim=c(0,19))
x<-seq(1.3,18.3,.2)
tcen<-(x-mean(div$Div.Time,na.rm=TRUE))/sd(div$Div.Time,na.rm=TRUE)
y<-0.70294+tcen*0.23001+tcen^2*-0.06316
lines(x,y,lwd=1.5)
legend(7,.2,c("Within species","Between species"),pch=19,col=cs,bty='n')
title("(A) Phylogentic model")

plot(tmu[-2],demog$GBS_FST[-2],col=cs[1],pch=19,xlab="Divergence time (gens./mu)",ylab=expression(F[ST]),ylim=c(0,1))
x<-seq(min(tmu[-2]),max(tmu[-2]),length.out=500)
tcen<-(x-mean(tmu[-2]))/sd(tmu[-2])
y<-0.65124+tcen*0.29458+tcen^2*-0.33264
lines(x,y,lwd=1.5)
title("(B) Demographic model")
dev.off()

## version two Ne units
tauc<-(demog$T_im[-2]-mean(demog$T_im[-2]))/sd(demog$T_im[-2])
tauc2<-tauc^2

pdf("fig_DivTimeFst.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4.5,5,2.5,2))
wsp<-as.numeric(div$Species1!=div$Species1.1)+1
cs<-c("deepskyblue4","tomato4")
plot(div$Div.Time,div$GBS.FST,col=cs[wsp],pch=19,xlab="Divergence time (MYA)",ylab=expression(F[ST]),ylim=c(0,1),xlim=c(0,19))
x<-seq(1.3,18.3,.2)
tcen<-(x-mean(div$Div.Time,na.rm=TRUE))/sd(div$Div.Time,na.rm=TRUE)
y<-0.70294+tcen*0.23001+tcen^2*-0.06316
lines(x,y,lwd=1.5)
legend(7,.2,c("Within species","Between species"),pch=19,col=cs,bty='n')
title("(A) Phylogentic model")

plot(tauc,demog$GBS_FST[-2],col=cs[1],pch=19,xlab="Divergence time (2Ne gens.)",ylab=expression(F[ST]),axes=FALSE,ylim=c(0,1))
axis(2)
axis(1,at=seq(-1,1,0.5),round(seq(-1,1,0.5) * sd(demog$T_im[-2]) + mean(demog$T_im[-2]),2))
box()
x<-seq(min(tauc),max(tauc),length.out=500)
#tcen<-(x-mean(tauc))/sd(tauc)
y<-0.25747+x*0.24960 +x^2*0.11034
lines(x,y,lwd=1.5)
title("(B) Demographic model")
dev.off()

pdf("sfig_DivTimeFstExtra.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,2.5,2))

plot(demog$T_im,demog$GBS_FST,col=cs[1],pch=19,xlab="Divergence time (2Ne generations)",ylab=expression(F[ST]),ylim=c(0,1))
title("(A) IM model, with outlier")

tmu<-demog$T_im/demog$Theta_im
tauc<-(tmu[-2]-mean(tmu[-2]))/sd(tmu[-2])
tauc2<-tauc^2
plot(tmu[-2],demog$GBS_FST[-2],col=cs[1],pch=19,xlab=expression(paste("Divergence time (generation/",mu,")",sep="")),ylab=expression(F[ST]),ylim=c(0,1))
x<-seq(min(tmu[-2]),max(tmu[-2]),length.out=500)
tcen<-(x-mean(tmu[-2]))/sd(tmu[-2])
y<-0.65124+tcen*0.29458+tcen^2*-0.33264
lines(x,y,lwd=1.5)
title(expression(paste("(B) IM model, time relative to ",mu,sep="")))

plot(demog$T_best,demog$GBS_FST,col=cs[1],pch=19,xlab="Divergence time (2Ne generations)",ylab=expression(F[ST]),ylim=c(0,1))
title("(C) Best model, time relative to Ne")

plot(demog$T_best/demog$Theta_best,demog$GBS_FST,col=cs[1],pch=19,xlab=expression(paste("Divergence time (generation/",mu,")",sep="")),ylab=expression(F[ST]),ylim=c(0,1))
title(expression(paste("(D) Best model, time relative to ",mu,sep="")))


dev.off()

## gene flow figure, again exclude crazy outlier
pdf("sfig_migs.pdf",width=4.5,height=4.5)
par(mar=c(4.5,4.5,1,1))
plot(demog$T_im[-2],demog$m12[-2],pch=19,xlab="Divergence time (2Ne generations)",ylab="Migration rate (2Ne migrants)",ylim=c(0,7.8))
points(demog$T_im[-2],demog$m21[-2],pch=19)
segments(demog$T_im[-2],demog$m12[-2],demog$T_im[-2],demog$m21[-2])
dev.off()

## best model results
demog<-read.table("demog_times.txt",sep=",",header=TRUE)
tauc<-(demog$T_best-mean(demog$T_best))/sd(demog$T_best)
tauc2<-tauc^2
o<-lm(demog$GBS_FST ~ tauc)
summary(o)
o<-lm(demog$GBS_FST ~ tauc+tauc2)
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.50374    0.11838   4.255  0.00377 **
#tauc         0.17060    0.09925   1.719  0.12933   
#tauc2       -0.19304    0.10257  -1.882  0.10186   
#Residual standard error: 0.2344 on 7 degrees of freedom
#Multiple R-squared:  0.3662,	Adjusted R-squared:  0.1851 
#F-statistic: 2.022 on 2 and 7 DF,  p-value: 0.2027


tauc<-(demog$T_best/demog$Theta_best-mean(demog$T_best/demog$Theta_best))/sd(demog$T_best/demog$Theta_best)
tauc2<-tauc^2
o<-lm(demog$GBS_FST ~ tauc)
summary(o)
o<-lm(demog$GBS_FST ~ tauc+tauc2)
#(Intercept)  0.53754    0.10963   4.903  0.00175 **
#tauc         0.33764    0.14045   2.404  0.04719 *
#tauc2       -0.23060    0.09578  -2.408  0.04693 *
#Residual standard error: 0.2142 on 7 degrees of freedom
#Multiple R-squared:  0.4705,	Adjusted R-squared:  0.3192
#F-statistic:  3.11 on 2 and 7 DF,  p-value: 0.108
#
