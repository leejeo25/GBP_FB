
theta_p<-0.2*pi
i=complex(real=0, imaginary=1)
W_p<- matrix( c(cos(theta_p), sin(theta_p),sin(-theta_p), cos(theta_p)), ncol=2 )
A_p<-matrix(c(3,1,2,4), nrow=2)
A<-t(W_p)%*%A_p
S1g<-array(dim=c(4,2,2))
theta<-0
i=complex(real=0, imaginary=1)
W<- matrix( c(cos(theta), sin(theta),sin(-theta), cos(theta)), ncol=2 )
C<-t(W)%*%A%*%t(A)%*%W
C[1,2]/sqrt(C[1,1]*C[2,2])
#sd11<-sd21<-sd31<-sd41<-sd51<-sd61<-sd71<-sd81<-c()
sd11<-sd21<-sd12<-sd22<-sd13<-sd23<-sd14<-sd24<-sd15<-sd25<-sd16<-sd26<-c()
h11<-h21<-h12<-h22<-h13<-h23<-h14<-h24<-h15<-h25<-h16<-h26<-c()
S<-W1<-W2<-W3<-c()
h1<-c(.2,.35,.7,.85);h2<-c(.25,.4,.75,.9);
for ( iter. in 1:4){
  h_1<-h1[iter.];h_2<-h2[iter.];
  sig_1<--(C[1,1]/2)*cos((h_1+h_1)*pi/2)*gamma(-h_1-h_1)
  sig_2<--(C[2,2]/2)*cos((h_2+h_2)*pi/2)*gamma(-h_2-h_2)
  sig_3<--(C[1,2]/2)*cos((h_1+h_2)*pi/2)*gamma(-h_1-h_2)
  SIG_0<-matrix(c(sig_1,sig_3,sig_3,sig_2),nrow=2)
  SIG<-8*W%*%SIG_0%*%t(W)
  S[iter.]<-SIG[1,2]/sqrt(SIG[1,1]*SIG[2,2])
  S1g[iter.,,]<-SIG
  n<-2^13
  f_1<-function(x){ W%*% diag( c(abs(x)^(h_1),abs(x)^(h_2) ))%*%t(W)};
  R<-function(x){(1/2)*(f_1(x+1) %*%SIG%*%f_1(x+1)+f_1(x-1) %*%SIG%*%f_1(x-1)-
                          2*f_1(x) %*%SIG%*%f_1(x))}
  R_s<-c(seq(0,n-1,1),seq(-n+2,-1,1))
  R_1<-sapply(R_s,R)
  lam1 <- Re(fft(R_1[1,]))
  lam2 <- Re(fft(R_1[2,]))
  lam3 <- Re(fft(R_1[3,]))
  lam4 <-Re( fft(R_1[4,]))
  
  W1[iter.]<-which(lam1<0.0001)
  
  W2[iter.]<-which(lam4<0.001)
  
  W3[iter.]<-which(lam3<0)
  D_12<-2*(lam1)^(-1/2)*(lam4)^(-1/2)*lam3
  D_21<-2*(lam1)^(-1/2)*(lam4)^(-1/2)*lam2
  DD<-rbind(rep(2,2*(n-1)),D_12,D_21,rep(2,2*(n-1)))
  C_new<-A_new<-array(dim=c(2*(n-1) ,2,2))
  m<-2*n-2
  for(o in 1:m){C_new[o,,]<-matrix(DD[,o],nrow=2);
  A_new[o,,]<-chol(C_new[o,,])}
  
  i0.<-30
  W_new<-array(rnorm(2*(2*n-2)*2*i0.,0,sd=sqrt(1/2)),
               dim=c(  i0.,2*(2*n-2), 2))
  W0_new<-array(,dim=c(i0.,2*(2*n-2)))
  Z<-array(,dim=c(i0.,2,2*n-2))
  for(j in 1:30){ W0_new[j,]<-W_new[j,,1]+i*W_new[j,,2] ;
  for(j_2 in 1:m){a<-2*j_2-1;b<-a+1;Z[j,,j_2]<-t(A_new[j_2,,])%*%
    rbind(W0_new[j,a],W0_new[j,b])
  }}
  final.1<-final.2<-array(,dim=c(i0.,2*n-2))
  for(j.3 in 1:30){ready<-matrix(nrow=2,ncol=2*n-2)
  ready<-rbind(lam1^(1/2),lam4^(1/2))*Z[j.3,,]
  final.1[j.3,] <-fft(ready[1,])/sqrt(2*n-2)
  final.2[j.3,]<-fft(ready[2,])/sqrt(2*n-2)  }
  IX2<-IX<-RX<-RX2<-array(dim=c(i0.,2,n ))
  final.1.B<-final.2.B<-array(dim=c(i0.,n))
  
  for(j.7 in 1:30){
    final.1.B[j.7,]<-sapply(seq(1:n), function(x)sum(final.1[j.7,1:x]))
    final.2.B[j.7,]<-sapply(seq(1:n), function(x)sum(final.2[j.7,1:x]))
  }
  IF<-RF<-IF2<-RF2<-array(dim=c(i0.,n))
  RF<-Re(final.1.B)
  IF<-Im(final.1.B)
  RF2<-Re(final.2.B)
  IF2<-Im(final.2.B)
  
  for(j.4 in 1:30){
    RX[j.4,,]<-rbind(RF[j.4,],RF2[j.4,]);
    RX2[j.4,,]<-W_p%*%f_1(1/n)%*%RX[j.4,,]
    IX[j.4,,]<-rbind(IF[j.4,],IF2[j.4,]);
    IX2[j.4,,]<-W_p%*%f_1(1/n)%*%IX[j.4,,]}
  
  ############hurst estimate ####################################
  
  fi.len<-length(RX2[1,1,])
  
  
  h.1<-h.2<-h.3<-h.4<-c()
  ei1.1<-ei1.2<-ei1.11<-ei1.12<-matrix(nrow=30,ncol=10)
  
  for(lo. in 1:30)
  {
    for(i. in 1:10){
      
      t1<-seq(1,fi.len-2*2^i.,1)
      t2<-seq(1+2^i.,fi.len-2^i.,1)
      t3<-seq(1+2*2^i.,fi.len,1)
      RXX<-matrix(nrow=2, ncol=length(t1)) 
      IXX<-matrix(nrow=2, ncol=length(t1))
      Q<-matrix(nrow=2,ncol=2)
      QI<-matrix(nrow=2,ncol=2)
      RXX[,]<-RX2[lo.,,t1]-2*RX2[lo.,,t2]+RX2[lo.,,t3]
      IXX[,]<-IX2[lo.,,t1]-2*IX2[lo.,,t2]+IX2[lo.,,t3]
      
      Q<-RXX[,]%*%t(RXX[,])/length(t1)
      QI<-IXX[,]%*%t(IXX[,])/length(t1)
      ei1.1[lo.,i.]<-max(eigen(Q)$values)
      ei1.2[lo.,i.]<-min(eigen(Q)$values) 
      ei1.11[lo.,i.]<-max(eigen(QI)$values)
      ei1.12[lo.,i.]<-min(eigen(QI)$values) 
    }
    x<-c(c(1:3)*log(2))
    h.1[lo.]<-lm(log(ei1.1[lo.,1:3])~x)$coefficients[2]/(2)
    h.2[lo.]<-lm(log(ei1.2[lo.,1:3])~x)$coefficients[2]/(2)
    h.3[lo.]<-lm(log(ei1.11[lo.,1:3])~x)$coefficients[2]/(2)
    h.4[lo.]<-lm(log(ei1.12[lo.,1:3])~x)$coefficients[2]/(2)
  }
  
  
  
  h11[iter.]<-h_1-mean(c(h.1,h.3))
  h21[iter.]<-h_2-mean(c(h.2,h.4))
  
  sd11[iter.]<-sd(c(h.1,h.3));
  sd21[iter.]<-sd(c(h.2,h.4))
  
  ######
  
  
  fi.len<-length(RX2[1,1,])
  
  
  h.1<-h.2<-h.3<-h.4<-c()
  ei1.1<-ei1.2<-ei1.11<-ei1.12<-matrix(nrow=30,ncol=10)
  
  for(lo. in 1:30)
  {
    for(i. in 1:10){
      
      t1<-seq(1,fi.len-2*2^i.,1)
      t2<-seq(1+2^i.,fi.len-2^i.,1)
      t3<-seq(1+2*2^i.,fi.len,1)
      RXX<-matrix(nrow=2, ncol=length(t1)) 
      IXX<-matrix(nrow=2, ncol=length(t1))
      Q<-matrix(nrow=2,ncol=2)
      QI<-matrix(nrow=2,ncol=2)
      RXX[,]<-RX2[lo.,,t1]-2*RX2[lo.,,t2]+RX2[lo.,,t3]
      IXX[,]<-IX2[lo.,,t1]-2*IX2[lo.,,t2]+IX2[lo.,,t3]
      
      Q<-RXX[,]%*%t(RXX[,])/length(t1)
      QI<-IXX[,]%*%t(IXX[,])/length(t1)
      ei1.1[lo.,i.]<-max(eigen(Q)$values)
      ei1.2[lo.,i.]<-min(eigen(Q)$values) 
      ei1.11[lo.,i.]<-max(eigen(QI)$values)
      ei1.12[lo.,i.]<-min(eigen(QI)$values) 
    }
    x<-c(c(3:5)*log(2))
    h.1[lo.]<-lm(log(ei1.1[lo.,3:5])~x)$coefficients[2]/(2)
    h.2[lo.]<-lm(log(ei1.2[lo.,3:5])~x)$coefficients[2]/(2)
    h.3[lo.]<-lm(log(ei1.11[lo.,3:5])~x)$coefficients[2]/(2)
    h.4[lo.]<-lm(log(ei1.12[lo.,3:5])~x)$coefficients[2]/(2)
  }
  
  
  
  h12[iter.]<-h_1-mean(c(h.1,h.3))
  h22[iter.]<-h_2-mean(c(h.2,h.4))
  
  sd12[iter.]<-sd(c(h.1,h.3));
  sd22[iter.]<-sd(c(h.2,h.4))
  ###
  
  
  
  fi.len<-length(RX2[1,1,])
  
  
  h.1<-h.2<-h.3<-h.4<-c()
  ei1.1<-ei1.2<-ei1.11<-ei1.12<-matrix(nrow=30,ncol=10)
  
  for(lo. in 1:30)
  {
    for(i. in 1:10){
      
      t1<-seq(1,fi.len-2*2^i.,1)
      t2<-seq(1+2^i.,fi.len-2^i.,1)
      t3<-seq(1+2*2^i.,fi.len,1)
      RXX<-matrix(nrow=2, ncol=length(t1)) 
      IXX<-matrix(nrow=2, ncol=length(t1))
      Q<-matrix(nrow=2,ncol=2)
      QI<-matrix(nrow=2,ncol=2)
      RXX[,]<-RX2[lo.,,t1]-2*RX2[lo.,,t2]+RX2[lo.,,t3]
      IXX[,]<-IX2[lo.,,t1]-2*IX2[lo.,,t2]+IX2[lo.,,t3]
      
      Q<-RXX[,]%*%t(RXX[,])/length(t1)
      QI<-IXX[,]%*%t(IXX[,])/length(t1)
      ei1.1[lo.,i.]<-max(eigen(Q)$values)
      ei1.2[lo.,i.]<-min(eigen(Q)$values) 
      ei1.11[lo.,i.]<-max(eigen(QI)$values)
      ei1.12[lo.,i.]<-min(eigen(QI)$values) 
    }
    x<-c(c(5:7)*log(2))
    h.1[lo.]<-lm(log(ei1.1[lo.,5:7])~x)$coefficients[2]/(2)
    h.2[lo.]<-lm(log(ei1.2[lo.,5:7])~x)$coefficients[2]/(2)
    h.3[lo.]<-lm(log(ei1.11[lo.,5:7])~x)$coefficients[2]/(2)
    h.4[lo.]<-lm(log(ei1.12[lo.,5:7])~x)$coefficients[2]/(2)
  }
  
  
  
  h13[iter.]<-h_1-mean(c(h.1,h.3))
  h23[iter.]<-h_2-mean(c(h.2,h.4))
  
  sd13[iter.]<-sd(c(h.1,h.3));
  sd23[iter.]<-sd(c(h.2,h.4))
  ####
  
  
  
  
  
  fi.len<-length(RX2[1,1,])
  
  
  h.1<-h.2<-h.3<-h.4<-c()
  ei1.1<-ei1.2<-ei1.11<-ei1.12<-matrix(nrow=30,ncol=12)
  
  for(lo. in 1:30)
  {
    for(i. in 1:10){
      
      t1<-seq(1,fi.len-2*2^i.,1)
      t2<-seq(1+2^i.,fi.len-2^i.,1)
      t3<-seq(1+2*2^i.,fi.len,1)
      RXX<-matrix(nrow=2, ncol=length(t1)) 
      IXX<-matrix(nrow=2, ncol=length(t1))
      Q<-matrix(nrow=2,ncol=2)
      QI<-matrix(nrow=2,ncol=2)
      RXX[,]<-RX2[lo.,,t1]-2*RX2[lo.,,t2]+RX2[lo.,,t3]
      IXX[,]<-IX2[lo.,,t1]-2*IX2[lo.,,t2]+IX2[lo.,,t3]
      
      Q<-RXX[,]%*%t(RXX[,])/length(t1)
      QI<-IXX[,]%*%t(IXX[,])/length(t1)
      ei1.1[lo.,i.]<-max(eigen(Q)$values)
      ei1.2[lo.,i.]<-min(eigen(Q)$values) 
      ei1.11[lo.,i.]<-max(eigen(QI)$values)
      ei1.12[lo.,i.]<-min(eigen(QI)$values) 
    }
    x<-c(c(7:9)*log(2))
    h.1[lo.]<-lm(log(ei1.1[lo.,7:9])~x)$coefficients[2]/(2)
    h.2[lo.]<-lm(log(ei1.2[lo.,7:9])~x)$coefficients[2]/(2)
    h.3[lo.]<-lm(log(ei1.11[lo.,7:9])~x)$coefficients[2]/(2)
    h.4[lo.]<-lm(log(ei1.12[lo.,7:9])~x)$coefficients[2]/(2)
  }
  
  
  
  h14[iter.]<-h_1-mean(c(h.1,h.3))
  h24[iter.]<-h_2-mean(c(h.2,h.4))
  
  sd14[iter.]<-sd(c(h.1,h.3));
  sd24[iter.]<-sd(c(h.2,h.4))
  
  ############################
  ###################################
  
  
  
  
  
}
det(S1g[1,,])
W1;W2;W3

R1<-R2<-c()
R1<-round(h11,digits=3)
R2<-round(sd11,digits=2)
Z1<-c()
for(y in 1:7){Z1<-c(Z1,R1[y],R2[y])    }
Z1

R1<-R2<-c()
R1<-round(h21,digits=3)
R2<-round(sd21,digits=2)
Z2<-c()
for(y in 1:7){Z2<-c(Z2,R1[y],R2[y])    }


R1<-round(h12,digits=3)
R2<-round(sd12,digits=2)
Z3<-c()
for(y in 1:7){Z3<-c(Z3,R1[y],R2[y])    }

R1<-R2<-c()
R1<-round(h22,digits=3)
R2<-round(sd22,digits=2)
Z4<-c()
for(y in 1:7){Z4<-c(Z4,R1[y],R2[y])    }


R1<-R2<-c()
R1<-round(h13,digits=3)
R2<-round(sd13,digits=2)
Z5<-c()
for(y in 1:7){Z5<-c(Z5,R1[y],R2[y])    }


R1<-R2<-c()
R1<-round(h23,digits=3)
R2<-round(sd23,digits=2)
Z6<-c()
for(y in 1:7){Z6<-c(Z6,R1[y],R2[y])    }


R1<-R2<-c()
R1<-round(h14,digits=3)
R2<-round(sd14,digits=2)
Z7<-c()
for(y in 1:7){Z7<-c(Z7,R1[y],R2[y])    }


R1<-R2<-c()
R1<-round(h24,digits=3)
R2<-round(sd24,digits=2)
Z8<-c()
for(y in 1:7){Z8<-c(Z8,R1[y],R2[y])    }


ZZ<-rbind(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)
DATA<-data.frame(ZZ)
table(DATA)
ZZ
DATA
