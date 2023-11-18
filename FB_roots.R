## Data analysis with roots data ##
## Data can be found in Ridout, M.S., C.G.B. Dem´etrio, and J.P. Hinde (1998). 

Xin<-c(64, 10, 13, 15, 21, 18, 24, 21, 23, 21, 17, 12,5, 2, 3, 0,0,1) # count for each number of roots
############################ Fitting fractional binomials ########
## Fitting FB-I ##
max.n<-length(Xin)
max1<-max.n-1
p.fun<-function(X){
  p<-X[1]
  H<-X[2]
  c<-X[3]
  P<-rep(0,max1)
  P[1]<-p+c
  d<-0
  for(i in 2:max1){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
  P[i]<-p+c*i^(2*H-2)-d }; P;   }
p.fun_0<-function(X){
  p<-X[1]
  H<-X[2]
  c<-X[3]
  P<-rep(0,max1)
  P[1]<-p
  d<-0
  for(i in 2:max1){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
  P[i]<-p-d }; P;   }
max2<-max1-1
likelihoodfunction12<-function(p0,h0,c0){if(1>p0&& p0>0&&1>h0&& h0>0 &&min(.5*(-2*p0+2^(2*h0-2)+(4*p0-p0*2^(2*h0)+2^(4*h0-4))^(1/2)),1-p0)>c0&& c0>=0 ){
  theo.p<-theo.p_0<-NULL
  PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
  theo.p_00<-p.fun_0(c(p0, h0, c0))
  theo.p<-p.fun(c(p0, h0, c0))
  theo.p_0<-1-cumsum(theo.p)
  PPm[1,]<-theo.p_00
  for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
  PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
   pre.like1<-log(c(1-sum(theo.p_00),PP1))*(Xin)
  -sum((pre.like1))} else{NA} }

est<-mle2(likelihoodfunction12, method = c("Nelder-Mead"),start=list( p0=.6,h0=.3, c0=.2))
est

Coefficients:
  p0        h0        c0 
0.2973391 0.7380297 0.3062589 

Log-likelihood: -671.22 
AIC=2*671.22 +2*3=1348.44


## FB- II ##
max.n<-length(Xin)
max1<-max.n-1
p.fun2<-function(X){
  H<-X[1]
  c<-X[2]
  P<-rep(0,max1)
  P[1]<-c
  d<-0
  for(i in 2:max1){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(c*(i-j)^(2*H-2))*P[j]}
  P[i]<-c*i^(2*H-2)-d }; P;   }
p.fun20<-function(X){
  H<-X[1]
  c<-X[2]
  la<-X[3]
  p<-la*(max.n^(2*H-2))
  P<-rep(0,max1)
  P[1]<-p
  d<-0
  for(i in 2:max1){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(c*(i-j)^(2*H-2))*P[j]}
  P[i]<-p-d }; P;   }
max2<-max1-1
likelihoodfunction22<-function(h0,c0,la0){if(1>h0&& h0>0.5 && 2^(2*h0-2)>c0&& c0>=0&& c0>la0&&la0>0 ){
  theo.p<-theo.p_0<-NULL
  PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
  theo.p<-p.fun2(c(h0, c0))
  theo.p_00<-p.fun20(c(h0, c0, la0))
  theo.p_0<-1-cumsum(theo.p)
  PPm[1,]<-theo.p_00
  for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
  PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
  
  pre.like1<-log(c(1-sum(theo.p_00),PP1))*(Xin)
  -sum((pre.like1))} else{NA} }

est<-mle2(likelihoodfunction22, method = c("Nelder-Mead"),start=list( h0=.9, c0=.5, la0=.1))
est


Call:
  mle2(minuslogl = likelihoodfunction22, start = list(h0 = 0.9, 
                                                      c0 = 0.5, la0 = 0.1), method = c("Nelder-Mead"))

Coefficients:
  h0        c0       la0 
0.9203628 0.5662389 0.4696990 

Log-likelihood: -672.05 
AIC=2*672.05 +2*3= 1350.1

## FB I* #
max.n<-length(Xin)
max1<-max.n-1
p.fun<-function(X){
  p<-X[1]
  H<-X[2]
  c<-X[3]
  P<-rep(0,max1)
  P[1]<-p+c
  d<-0
  for(i in 2:max1){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
  P[i]<-p+c*i^(2*H-2)-d }; P;   }
max2<-max1-1
likelihoodfunction1<-function(p0,h0,c0){if(1>p0&& p0>0&&1>h0&& h0>0 &&min(.5*(-2*p0+2^(2*h0-2)+(4*p0-p0*2^(2*h0)+2^(4*h0-4))^(1/2)),1-p0)>c0&& c0>=0 ){
  theo.p<-theo.p_0<-NULL
  PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
  theo.p<-p.fun(c(p0, h0, c0))
  theo.p_0<-1-cumsum(theo.p)
  PPm[1,]<-theo.p
  for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
  PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
  pre.like1<-log(c(1-sum(theo.p),PP1))*(Xin)
  -sum((pre.like1))} else{NA} }

est<-mle2(likelihoodfunction1, method = c("Nelder-Mead"),start=list( p0=.7,h0=.2, c0=.1))
est

Call:
  mle2(minuslogl = likelihoodfunction1, start = list(p0 = 0.7, 
                                                     h0 = 0.2, c0 = 0.1), method = c("Nelder-Mead"))

Coefficients:
  p0           h0           c0 
9.014336e-08 7.668361e-01 7.007610e-01 

Log-likelihood: -708.12 
AIC=2*708.12 +2*3=1422.24

## fract II* ##
max.n<-length(Xin)
max1<-max.n-1
p.fun2<-function(X){
  H<-X[1]
  c<-X[2]
  P<-rep(0,max1)
  P[1]<-c
  d<-0
  for(i in 2:max1){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(c*(i-j)^(2*H-2))*P[j]}
  P[i]<-c*i^(2*H-2)-d }; P;   }

max2<-max1-1
likelihoodfunction2<-function(h0,c0){if(1>h0&& h0>0.5 && 2^(2*h0-2)>c0&& c0>=0 ){
  theo.p<-theo.p_0<-NULL
  PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
  theo.p<-p.fun2(c(h0, c0))
  theo.p_0<-1-cumsum(theo.p)
  PPm[1,]<-theo.p
  for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
  PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
  
  pre.like1<-log(c(1-sum(theo.p),PP1))*(Xin)
  -sum((pre.like1))
} else{NA} }

est<-mle2(likelihoodfunction2, method = c("Nelder-Mead"),start=list( h0=.9, c0=.2))
est

Call:
  mle2(minuslogl = likelihoodfunction2, start = list(h0 = 0.9, 
                                                     c0 = 0.6), method = c("Nelder-Mead"))

Coefficients:
  h0        c0 
0.7665931 0.7009977 

Log-likelihood: -708.12 
AIC=2*708.12 +2*2=1420.24

########################################################
############### Fitting zero inflated models #############
# zero-inflated nbinom
library("pscl")
len.<-length(Xin)-1
data<-rep(seq(0,len.,1),Xin)

m1 <- zeroinfl(data~1|1, dist ="negbin")
summary(m1)
Call:
  zeroinfl(formula = datapois ~ 1 | 1, dist = "negbin")

Pearson residuals:
  Min       1Q   Median       3Q      Max 
-1.25925 -1.01035 -0.01475  0.73195  2.97206 

Count model coefficients (negbin with log link):
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.88539    0.03538  53.293  < 2e-16 ***
  Log(theta)   2.29962    0.28214   8.151 3.62e-16 ***
  
  Zero-inflation model coefficients (binomial with logit link):
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -1.1962     0.1464  -8.169  3.1e-16 ***
  ---
  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

Theta = 9.9704 
Number of iterations in BFGS optimization: 6 
Log-likelihood: -676.3 on 3 Df

AIC=2*676.3+2*3=1358.6
 
MLE of mu= exp(1.88539)
MLE of overdispersion parameter=9.9704
MLE of zprob=   1/(1/exp(-1.1962 )+1)   

# zero-inflated Poisson ##
m1 <- zeroinfl(data~1|1, dist ="pois")
summary(m1)
Call:
  zeroinfl(formula = datapois ~ 1 | 1, dist = "pois")

Pearson residuals:
  Min       1Q   Median       3Q      Max 
-1.40498 -1.12728 -0.01646  0.81666  3.31600 

Count model coefficients (poisson with log link):
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.89043    0.02718   69.56   <2e-16 ***
  
  Zero-inflation model coefficients (binomial with logit link):
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -1.1746     0.1437  -8.173 3.02e-16 ***
  ---
  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

Number of iterations in BFGS optimization: 4 
Log-likelihood: -688.8 on 2 Df

AIC=2*688.8 +2*2=1381.6
MLEs:
zprob= 1/(1/exp(-1.1746 )+1)
lambda=exp(1.89043 )


######################################################
###################plot in Figure 8 ###
library("emdbook")
library("VGAM")
root_dist_nbinom_zero<-dzinbinom(seq(0,len.,1),6.588924 , 9.9704 , 0.2321519, log = FALSE)
root_dist_pois_zero<-dzipois(seq(0,len.,1), lambda=6.622216, pstr0 =  0.2360245, log = FALSE)

p0<-0.2973391;h0<- 0.7380297;c0<- 0.3062589
theo.p<-theo.p_0<-NULL
PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
theo.p_00<-p.fun_0(c(p0, h0, c0))
theo.p<-p.fun(c(p0, h0, c0))
theo.p_0<-1-cumsum(theo.p)
PPm[1,]<-theo.p_00
for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
root_dist<-c(1-sum(theo.p_00),PP1)

h0<-0.9203628 ;c0<-0.5662389;la0<- 0.4696990
theo.p<-theo.p_0<-NULL
PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
theo.p<-p.fun2(c(h0, c0))
theo.p_00<-p.fun20(c(h0, c0, la0))
theo.p_0<-1-cumsum(theo.p)
PPm[1,]<-theo.p_00
for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
root_dist_frac2<-c(1-sum(theo.p_00),PP1)

p0<-9.014336e-08;h0<- 7.668361e-01;c0<- 7.007610e-01 
theo.p<-theo.p_0<-NULL
PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
theo.p<-p.fun(c(p0, h0, c0))
theo.p_0<-1-cumsum(theo.p)
PPm[1,]<-theo.p
for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
root_dist_frac1_s<-c(1-sum(theo.p),PP1)

h0<-0.7665931;c0<- 0.7009977  
theo.p<-theo.p_0<-NULL
PPm<-matrix(0, ncol=max1, nrow=max1); PP1<-c()
theo.p<-p.fun2(c(h0, c0))
theo.p_0<-1-cumsum(theo.p)
PPm[1,]<-theo.p
for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  } 
PP1<-PPm[,max1]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
root_dist_frac2_s<-c(1-sum(theo.p),PP1)


par(mar=c(4,5,2,1))
data<-rep(seq(0,17,1 ), Xin)
hist(data, breaks=18, probability = TRUE)
lines(seq(0,17,1 ),(root_dist), type="b",cex.lab=2,cex.axis=2, xlab="", lwd=2, col="red", ylab="Probability")
lines(seq(0,17,1 ), (root_dist_frac1_s), type="b",cex.lab=1.5,cex.axis=1.5, xlab="", lwd=2, col="purple")
lines(seq(0,17,1 ),(root_dist_frac2), type="b",cex.lab=2,cex.axis=2, lwd=2, col="orange", xlab="", ylab="Probability")
lines(seq(0,17,1 ),(root_dist_frac2_s), type="b",cex.lab=1.5,cex.axis=1.5, xlab="", lwd=2, col="cyan")
lines(seq(0,17,1 ),root_dist_pois_zero,lwd=2, type="b", col="green")
lines(seq(0,17,1 ),root_dist_nbinom_zero,lwd=2, type="b", col="blue")


legend("topright",col=c("red","orange", "purple", "cyan", "green", "blue"),legend= c("FB-I", "FB-II",
   "FB-I*", "FB-II*",  "ZI-pois","ZI-NB"), lty=c(1,1,1),
       lwd=c(2,2,2,2,2,2), cex=2)




