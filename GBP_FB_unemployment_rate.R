# DATA from https://fred.stlouisfed.org/series/UNRATE

unrate<-read.csv("UNRATE.csv")

D<-unrate$DATE
X<-unrate$UNRATE
X.ts<-NULL
X.ts <- ts(X, start=c(1,1), freq = 12)
d.r<-decompose(X.ts)
Xd<-diff(X)
d.r.r<-c()
d.r.r<-d.r$random[7:894]
date<-D[7:894]
## Figure 5 ##
d<-NULL
d<-seq(1,888 ,by=216)
plot(X, type="l" ,xlab="Date" ,xaxt="n", cex.lab=1.5,cex.axis=1.5,ylab="Unemployment rate")
plot(d.r.r, type="l" ,xlab="Date" ,xaxt="n", ylim=c(-1.1,1.1), cex.lab=1.5,cex.axis=1.5,ylab="Detrended unemployment rate")
axis(1, at =d, labels=date[d], cex.axis=1.5)
abline(h=c(.08, .21,.28), col="red", lwd=2)
## Figure 6 ##
I<-rep(0, length(d.r.r))
I[ which(d.r.r>.08)]<-1
a<-acf(I)
plot(a, main="",  cex.lab=1.5,cex.axis=1.5) 
dif<-c()
dif<-diff(which(d.r.r>.08))
ad<-acf(dif)
plot(ad, main="",  cex.lab=1.5,cex.axis=1.5) 
I<-rep(0, length(d.r.r))
I[ which(d.r.r>.21)]<-1
a<-acf(I)
plot(a, main="",  cex.lab=1.5,cex.axis=1.5) 
dif<-c()
dif<-diff(which(d.r.r>.21))
ad<-acf(dif)
plot(ad, main="",  cex.lab=1.5,cex.axis=1.5) 
I<-rep(0, length(d.r.r))
I[ which(d.r.r>.28)]<-1
a<-acf(I)
plot(a, main="",  cex.lab=1.5,cex.axis=1.5) 
dif<-c()
dif<-diff(which(d.r.r>.28))
ad<-acf(dif)
plot(ad, main="",  cex.lab=1.5,cex.axis=1.5) 


## estimation for fracI with ##
# with cutoff .08 #
dif<-diff(which(d.r.r>.08))
Xin<-dif
est<-est.p<-est.h<-est.c<-NULL;
max.n<-max(Xin)
likelihoodfunction1<-function(p0,h0,c0){if(1>p0&& p0>0&&1>h0&& h0>0 &&1>c0&& c0>=0 ){
  p<-p0; H<-h0; c<-c0
  P<-rep(0,max.n)
  P[1]<-p+c
  d22<-c()
  f<-function(x){p+c*(x)^(2*H-2)}
  d22<-sapply( seq(max.n,1,-1), f )
  for(i in 2:max.n){ i1<-i-1; mii<-max.n-i+2 ;
  P[i]<-p+c*i^(2*H-2)-sum(d22[mii:max.n]*P[1:i1] )};
  
  pre.like1<-P[Xin];
  -sum(log(c(pre.like1)))} else{NA} }
library("bbmle")
est<-mle2(likelihoodfunction1,method = c("Nelder-Mead"),start=list(p0=.2, h0=.6, c0=.2))
est
est.p<-est@coef[1];est.h<-est@coef[2];est.c<-est@coef[3]
est.p;est.h;est.c
# with cut off .21
dif<-diff(which(d.r.r>.21))
Xin<-dif
est<-est.p<-est.h<-est.c<-NULL;
max.n<-max(Xin)
est<-mle2(likelihoodfunction1,method = c("Nelder-Mead"),start=list(p0=.2, h0=.6, c0=.2))
est
est.p<-est@coef[1];est.h<-est@coef[2];est.c<-est@coef[3]
est.p;est.h;est.c
# with cut off .28
dif<-diff(which(d.r.r>.28))
Xin<-dif
est<-est.p<-est.h<-est.c<-NULL;
max.n<-max(Xin)
est<-mle2(likelihoodfunction1,method = c("Nelder-Mead"),start=list(p0=.2, h0=.6, c0=.2))
est
est.p<-est@coef[1];est.h<-est@coef[2];est.c<-est@coef[3]
est.p;est.h;est.c

# estimation for fracII##
# with cut off .08
dif<-diff(which(d.r.r>.08))
Xin<-dif
max.n<-max(Xin)
likelihoodfunction1<-function(p0,h0,c0){if(1>h0&& h0>0 &&1>c0&& c0>=0 ){
  H<-h0; c<-c0
  P<-rep(0,max.n)
  P[1]<-c
  d22<-c()
  f<-function(x){c*(x)^(2*H-2)}
  d22<-sapply( seq(max.n,1,-1), f )
  for(i in 2:max.n){ i1<-i-1; mii<-max.n-i+2 ;
  P[i]<-c*i^(2*H-2)-sum(d22[mii:max.n]*P[1:i1] )};
  
  pre.like1<-P[Xin];
  -sum(log(c(pre.like1)))} else{NA} }
est<-mle2(likelihoodfunction1,method = c("Nelder-Mead"),start=list(h0=.9, c0=.2))
est
# with cut off .21
dif<-diff(which(d.r.r>.21))
Xin<-dif
max.n<-max(Xin)
likelihoodfunction1<-function(p0,h0,c0){if(1>h0&& h0>0 &&1>c0&& c0>=0 ){
  H<-h0; c<-c0
  P<-rep(0,max.n)
  P[1]<-c
  d22<-c()
  f<-function(x){c*(x)^(2*H-2)}
  d22<-sapply( seq(max.n,1,-1), f )
  for(i in 2:max.n){ i1<-i-1; mii<-max.n-i+2 ;
  P[i]<-c*i^(2*H-2)-sum(d22[mii:max.n]*P[1:i1] )};
  pre.like1<-P[Xin];
  -sum(log(c(pre.like1)))} else{NA} }
est<-mle2(likelihoodfunction1,method = c("Nelder-Mead"),start=list(h0=.9, c0=.2))
est
# with cut off .28
dif<-diff(which(d.r.r>.28))
Xin<-dif
max.n<-max(Xin)
likelihoodfunction1<-function(p0,h0,c0){if(1>h0&& h0>0 &&1>c0&& c0>=0 ){
  H<-h0; c<-c0
  P<-rep(0,max.n)
  P[1]<-c
  d22<-c()
  f<-function(x){c*(x)^(2*H-2)}
  d22<-sapply( seq(max.n,1,-1), f )
  for(i in 2:max.n){ i1<-i-1; mii<-max.n-i+2 ;
  P[i]<-c*i^(2*H-2)-sum(d22[mii:max.n]*P[1:i1] )};
  
  pre.like1<-P[Xin];
  -sum(log(c(pre.like1)))} else{NA} }
est<-mle2(likelihoodfunction1,method = c("Nelder-Mead"),start=list(h0=.9, c0=.2))
est

## estimation of Bernoulli  ##
 # (change the cuttoff below as .08, .21)
dif<-diff(which(d.r.r>.28))
Xin<-dif
data_X<-Xin-1
install.packages("fitdistrplus")
library("fitdistrplus")
fw <- fitdist(data_X, "geom")
summary(fw)


## estimation of Markov chain ##
# .28
table(Xin-1)/length(Xin)
0        1        2        3        4        5        8        9       10       12 
0.531250 0.093750 0.015625 0.015625 0.015625 0.015625 0.031250 0.015625 0.015625 0.046875 
15       17       18       23       25       39       40       48       68       85 
0.015625 0.015625 0.015625 0.015625 0.031250 0.015625 0.015625 0.015625 0.015625 0.015625 
86      221 
0.015625 0.015625 

table(Xin)
1   2   3   4   5   6   9  10  11  13  16  18  19  24  26  40  41  49  69  86  87 222 
34   6   1   1   1   1   2   1   1   3   1   1   1   1   2   1   1   1   1   1   1   1 
# p(1,1)= 0.531250 (MLE of p),  p(1,0)= 1-    0.531250

data_X_2<-Xin[which(Xin>1)]-2
fw <- fitdist(data_X_2, "geom")
summary(fw)
estimate  Std. Error
prob 0.03745318(MLE of q) 0.006704103
Loglikelihood:  -127.971   AIC:  257.9419   BIC:  259.3431 
# Loglikelihood:  -172.2073
-127.971+ 34*log(   0.531250)+ (length(Xin)-34)*log(1-   0.531250)

lines(seq(0,max(Xin),1), c(   0.531250 , (1-   0.531250 )*dgeom(seq(0,max(Xin)-1,1) ,  0.03745318)), type="l", col="orange", lty=1, lwd=2)

# cutoff .21 
table(Xin-1)/length(Xin)

0          1          2          3          4          7          8          9 
0.46590909 0.10227273 0.01136364 0.04545455 0.01136364 0.01136364 0.04545455 0.02272727 
10         11         12         13         14         15         16         20 
0.01136364 0.01136364 0.02272727 0.01136364 0.03409091 0.02272727 0.03409091 0.01136364 
21         25         26         29         31         36         53         57 
0.01136364 0.01136364 0.01136364 0.01136364 0.01136364 0.01136364 0.01136364 0.01136364 
58         67         92 
0.01136364 0.01136364 0.01136364 

table(Xin)

1  2  3  4  5  8  9 10 11 12 13 14 15 16 17 21 22 26 27 30 32 37 54 58 59 68 93 
41  9  1  4  1  1  4  2  1  1  2  1  3  2  3  1  1  1  1  1  1  1  1  1  1  1  1
# p(1,1)=     0.46590909,p(1,0)= 1-     0.46590909

data_X_2<-Xin[which(Xin>1)]-2
fw <- fitdist(data_X_2, "geom")
summary(fw)
estimate  Std. Error
prob 0.06048906 0.008550022
Loglikelihood:  -177.3977   AIC:  356.7954   BIC:  358.6456 
# Loglikelihood:    -238.1899
-177.3977+ 41*log(  0.46590909)+ (length(Xin)-41)*log(1-  0.46590909)

lines(seq(0,max(Xin),1), c( 0.46590909 , (1- 0.46590909 )*dgeom(seq(0,max(Xin)-1,1) ,  0.06048906)), type="l", col="orange", lty=1, lwd=2)

# cutoff .08 

table(Xin-1)/length(Xin)

0           1           2           3           4           5           6 
0.571955720 0.092250923 0.033210332 0.036900369 0.047970480 0.044280443 0.044280443 
7           8           9          10          11          12          13 
0.025830258 0.025830258 0.025830258 0.014760148 0.007380074 0.018450185 0.003690037 
14          22 
0.003690037 0.003690037 

table(Xin)

1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  23 
155  25   9  10  13  12  12   7   7   7   4   2   5   1   1   1 
# p(1,1)=   0.571955720 ,p(1,0)= 1-   0.571955720

Xin[which(Xin>1)]-2
data_X_2<-Xin[which(Xin>1)]-2
data_X_2
fw <- fitdist(data_X_2, "geom")
summary(fw)
Parameters : 
  estimate Std. Error
prob 0.1946309 0.01621701
Loglikelihood:  -293.7496   AIC:  589.4993   BIC:  592.2529 
# Loglikelihood:    -478.7764
-293.7496+ 155*log(0.571955720)+ (length(Xin)-155)*log(1-0.571955720)


