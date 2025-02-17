### Simulate the Data####
###################################################
rm(list = ls())
#ls()
########METHOD 1###########################
### Parameters #######
phi1_true<- c(0.25,0.75)
phi2_true<- c(0.1,0.9)
pi_true<- c(0.4,0.6)
N<- 1000
G<- 100


phi1<- phi1_true
phi2<- phi2_true
Pi<- pi_true
lambda1<- 2
lambda2<- 5
### 1.Generate 2 vectors of U####################
U1<- rbinom(n=G,size=1,prob = phi1[1])
U2<- rbinom(n=G,size=1,prob = phi2[1])

#U1
#U2
#U3<- rbinom(n=10,size=1,prob = 0.1)
### 2.Generate cluster assignments###
Z<- sample(2,N,replace = TRUE,prob =Pi )
#Z
### 3. Generate Y1,...,YN###

Y<- matrix(rep(0,G*N),ncol=N)
#Y
for (i in 1:N){
  if(Z[i]==1){
    for(j in 1:length(U1)){
      if(U1[j]==1){Y[j,i]<- 0}
      else {Y[j,i]<- rpois(1,lambda = lambda1)}
    }
  }
  else if(Z[i]==2){
    for(j in 1:length(U1)){
      if(U2[j]==1){Y[j,i]<- 0}
      else {Y[j,i]<- rpois(1,lambda = lambda2)}
    }
    
  }
  
}
dim(Y)

#Y[,1]
#print(Y)
Y1<- t(Y)
Y
plot(hclust(dist(Y1)))
#ls()
sum(U1)/length(U1)
sum(U2)/length(U2)
sum(Z==1)/length(Z)


rm(list=c("i","j","lambda1","lambda2","M","N","phi1","phi2","Pi","Z"))
ls()
##################################################
##################################################
### EM Algorithm####
ZIPmix2comp<- function(Y,Phi1,Phi2,Pn,L1,L2,eps=1e-8,maxiter=1000){
  phi1<- Phi1
  phi2<- Phi2
  Pi<- Pn
  lambda1<- L1
  lambda2<- L2
  M<- nrow(Y)
  N<- ncol(Y)
  dl<- 1+eps
  iter<- 0
  ll<- rep(0,maxiter+1)

  
  b2<- matrix(rep(0,N*M),nrow=M)
  
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){
      if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1))}
      else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1)}
    }
  }
  #print(b2)
  #b2
  ######################################
  b3<- matrix(rep(0,N*M),nrow=M)
  
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){
      if(Y[i,j]==0){b3[i,j]<- phi2+((1-phi2)*exp(-lambda2))}
      else {b3[i,j]<- (1-phi2)*dpois(Y[i,j],lambda2)}
    }
  }
  #print(b3)
  #b3
  ############################################
  ############################################
  b22<- Pi*apply(b2,2,prod)
  
  b33<- (1-Pi)*apply(b3,2,prod)
  
  b44<- (Pi*b2+(1-Pi)*b3)
  #print(b44)
  
  l<- sum(sum(log(b44)))
  #print(l)
  
  #b22<- Pi*b2
  
  #b33<- (1-Pi)*b3
  
  #l<- sum(log(b22+b33))
  
  #l<- sum(sum(log(b21+b31)))
  print(l)
  
  while(abs(dl)>eps & iter<maxiter){
  #while(abs(dl)<eps & iter<maxiter){
      
    iter<- iter+1
    ll[iter]<- l
    postprobs<- b22/(b22+b33)
    Pi<- mean(postprobs)
    ZZ<- cbind(postprobs,1-postprobs)
    
    
    ukg<- matrix(rep(0,M*N),nrow=M)
    for(j in 1:ncol(Y)){
      for(i in 1:nrow(Y)){
        if(Y[i,j]==0){ukg[i,j]<- (Pi*phi1+(1-Pi)*phi2)/(Pi*(phi1+(1-phi1)*exp(-lambda1))+
                                                                  (1-Pi)*(phi2+(1-phi2)*exp(-lambda2)))}
        else {ukg[i,j]<- 0}
      }
      
    }
    #ukg[,1:3]
    c<- apply(ukg,2,sum)
    #print(postprobs*c)
    phi1<- sum(postprobs*c)/(M*sum(postprobs))
    phi2<- sum((1-postprobs)*c)/(M*sum(1-postprobs))
    
    d1<- (1-ukg)*Y
    d2<- apply(d1,2,sum)
    #lambda1<- sum(postprobs*d2)/(M*sum(postprobs))
    #lambda2<- sum((1-postprobs)*d2)/(M*sum(1-postprobs))
    lambda1<- sum(postprobs*t(d1))/sum((postprobs*t(1-ukg)))
    lambda2<- sum((1-postprobs)*t(d1))/sum((1-postprobs)*t(1-ukg))
                                    
###############################################################################
    b2<- matrix(rep(0,N*M),nrow=M)
    
    for(j in 1:ncol(Y)){
      for(i in 1:nrow(Y)){
        if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1))}
        else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1)}
      }
    }
    
    #b2
    ######################################
    b3<- matrix(rep(0,N*M),nrow=M)
    
    for(j in 1:ncol(Y)){
      for(i in 1:nrow(Y)){
        if(Y[i,j]==0){b3[i,j]<- phi2+((1-phi2)*exp(-lambda2))}
        else {b3[i,j]<- (1-phi2)*dpois(Y[i,j],lambda2)}
      }
    }
    
    #b3
    ############################################
    ############################################
    b22<- Pi*apply(b2,2,prod)
    
    b33<- (1-Pi)*apply(b3,2,prod)
    
    b44<- (Pi*b2)+((1-Pi)*b3)
    #b21<- Pi*apply(b2,2,sum)
    
    #b31<- (1-Pi)*apply(b3,2,sum)
    
    oldl<- l
    l<- sum(sum(log(b44)))
    #l<- sum(log(b21+b31))
    
    dl<- l-oldl
    #print(dl)
    #print(dl>eps)
  }
  cat("number of iterations=", iter, "\n")
  iter<- iter+1
  ll[iter]<- l
  
  #postprobs<- cbind(postprobs,1-postprobs)
  #colnames(postprobs)<- c(paste("comp", ".", 1:2, sep = ""))
  myoutput<- list(Pi=c(Pi,1-Pi),Phi1=c(phi1,1-phi1),Phi2=c(phi2,1-phi2),lambda=c(lambda1,
                  lambda2),loglik=l,all.loglik=ll[1:iter],
                  restart=0,ft="ZIPmixEM")
  
  class(myoutput)<- "mixEM"
  myoutput
  
}




output1<- ZIPmix2comp(Y,Phi1=0.25,Phi2=0.15,Pn=0.4,L1=2,L2=5,eps=1e-8,maxiter=1000)

output1$Pi
pi_true

output1$Phi1
phi1_true

output1$Phi2
phi2_true


output1$lambda

lambda1
lambda2

library(mixtools)
attach(faithful)
wait1 <- normalmixEM(waiting, lambda = .5, mu = c(55, 80), sigma = 5)
wait1
###################################################################################
U1
U2

(U1==1)
a<- which(U1==1)
length(a)
30/100


b<- which(U2==1)
length(b)
12/100
