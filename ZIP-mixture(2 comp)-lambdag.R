### Simulate the Data####
###################################################
rm(list = ls())
#ls()
########METHOD 1###########################
### Parameters #######
phi1<- c(0.25,0.75)
phi2<- c(0.1,0.9)
Pi<- c(0.4,0.6)
lambda1<- rgamma(n=100,shape=0.5,rate=1/2)
lambda2<- rgamma(n=100,shape = 1,rate = 1/2)
lambda1
lambda2
### 1.Generate 2 vectors of U####################
U1<- rbinom(n=100,size=1,prob = phi1[1])
U2<- rbinom(n=100,size=1,prob = phi2[1])

#U1
#U2
#U3<- rbinom(n=10,size=1,prob = 0.1)
### 2.Generate cluster assignments###
Z<- sample(2,1000,replace = TRUE,prob =Pi )
#Z
### 3. Generate Y1,...,YN###
M<- length(U1)
N<- length(Z)
Y<- matrix(rep(0,M*N),ncol=N)
#Y
for (i in 1:N){
  if(Z[i]==1){
    for(j in 1:length(U1)){
      if(U1[j]==1){Y[j,i]<- 0}
      else {Y[j,i]<- rpois(1,lambda = lambda1[j])}
    }
  }
  else if(Z[i]==2){
    for(j in 1:length(U1)){
      if(U2[j]==1){Y[j,i]<- 0}
      else {Y[j,i]<- rpois(1,lambda = lambda2[j])}
    }
    
  }
  
}

Y[,1]
#print(Y)
#Y1<- t(Y)
#plot(hclust(dist(Y1)))
#ls()
#t(Y)
sum(U1)/length(U1)
sum(U2)/length(U2)
sum(Z==1)/length(Z)

rm(list=c("i","j","M","N","phi1","phi2","Pi","Z"))
ls()
##########################################################################

### Simulate the Data: Method 3 ####
###################################################
rm(list = ls())
#ls()
###################################################
### Function of ZIP Simulation Function###
ZIPDist<- function(n,pi,theta){
  U<- rbinom(size = 1,n=n,p=pi)
  Y<- NULL
  for(i in 1:n){
    if (U[i]==1){
      Y<- c(Y,0)
    }
    else {
      Y<- c(Y,rpois(1,lambda = theta))
    }
  }
  return(Y)
}

#sim1<- ZIPDist(n=10000,pi=0.2,theta=20)
#sim1
#hist(sim1)
####################################################################
########METHOD 1###########################
### Parameters #######
M<- 300
N<- 1000  
phi1<- c(0.25,0.75)
phi2<- c(0.1,0.9)
Pi<- c(0.4,0.6)
lambda1<- rgamma(n=M,shape=0.5,rate=1/2)
lambda2<- rgamma(n=M,shape = 1.5,rate = 1/2)
#lambda2<- rgamma(n=M,shape = 5,rate = 1/2)

#lambda1
#lambda2

#Phi1_vec <- NULL
#Phi2_vec <- NULL
#for(s in 1:20){
  #phi1<- c(0.25,0.75)
  #phi2<- c(0.1,0.9)
  #Pi<- c(0.4,0.6)
  #lambda1<- lambda11
 # lambda2<- lambda22
#  print(s)


### 1.Generate  vectors of Z For the Clusters####################
Z<- sample(2,N,replace = TRUE,prob =Pi )

Y<- matrix(rep(0,M*N),ncol=N)
#Y
for (i in 1:N){
  if(Z[i]==1){
    for(j in 1:M){
      Y[j,i]<- ZIPDist(n=1,pi=phi1[1],theta=lambda1[j])
      
    }
  }
  else if(Z[i]==2){
    for(j in 1:M){
      Y[j,i]<- ZIPDist(n=1,pi=phi2[1],theta=lambda2[j])
      
    }
    
  }
  
}

Y[,1]
#print(Y)
Y1<- t(Y)
#print(Y1)
#Y1[,1:10]
plot(hclust(dist(Y1)))
#ls()
#t(Y)
#sum(U1)/length(U1)
#sum(U2)/length(U2)
sum(Z==1)/length(Z)

rm(list=c("i","j","M","N","phi1","phi2","Pi","Z"))
#ls()
#####################################################################################
####################################################################################
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
      if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1[i]))}
      else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1[i])}
    }
  }
  
  #b2
  ######################################
  b3<- matrix(rep(0,N*M),nrow=M)
  
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){
      if(Y[i,j]==0){b3[i,j]<- phi2+((1-phi2)*exp(-lambda2[i]))}
      else {b3[i,j]<- (1-phi2)*dpois(Y[i,j],lambda2[i])}
    }
  }
  
  #b3
  ############################################
  ############################################
  b22<- Pi*apply(b2,2,prod)
  
  b33<- (1-Pi)*apply(b3,2,prod)
  
  #b21<- Pi*apply(b2,2,sum)
  
  #b31<- (1-Pi)*apply(b3,2,sum)
  
  l<- sum(log(b22+b33))
  #l<- sum(log(b21+b31))
  
  
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
        if(Y[i,j]==0){ukg[i,j]<- (Pi*phi1+(1-Pi)*phi2)/(Pi*(phi1+(1-phi1)*exp(-lambda1[i]))+
                                                          (1-Pi)*(phi2+(1-phi2)*exp(-lambda2[i])))}
        else {ukg[i,j]<- 0}
      }
      
    }
    #ukg[,1:3]
    
    comb1<- postprobs*t(ukg)
    phi1<- sum(comb1)/(M*sum(postprobs))
    #print(phi1)
    comb2<- (1-postprobs)*t(ukg)
    phi2<- sum(comb2)/(M*sum(1-postprobs))
    #print(phi2)
    
    
    
    #c<- apply(ukg,2,sum)
    #print(postprobs*c)
    #phi1<- sum(postprobs*c)/(M*sum(postprobs))
    #phi2<- sum((1-postprobs)*c)/(M*sum(1-postprobs))
    
    d1<- (1-ukg)*Y
    d12<- apply((postprobs*t(d1)),2,sum)
    d13<- apply((postprobs*t(1-ukg)),2,sum)
    
    
    d22<- apply(((1-postprobs)*t(d1)),2,sum)
    d23<- apply(((1-postprobs)*t(1-ukg)),2,sum)
    
    #for (i in 1:length(d13)){
     # if (d13[i]==0) {d13[i]<- 0.000001}
      #if (d23[i]==0) {d23[i]<- 0.000001}
    #}
    lambda1<- d12/d13
    lambda2<- d22/d23
    
      
    
    #lambda1<- sum(postprobs*d2)/(M*sum(postprobs))
    #lambda2<- sum((1-postprobs)*d2)/(M*sum(1-postprobs))
    
    #lambda1<- sum(postprobs*t(d1))/sum((postprobs*t(1-ukg)))
    #lambda2<- sum((1-postprobs)*t(d1))/sum((1-postprobs)*t(1-ukg))
    
    ###############################################################################
    b2<- matrix(rep(0,N*M),nrow=M)
    
    for(j in 1:ncol(Y)){
      for(i in 1:nrow(Y)){
        if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1[i]))}
        else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1[i])}
      }
    }
    
    #b2
    ######################################
    b3<- matrix(rep(0,N*M),nrow=M)
    
    for(j in 1:ncol(Y)){
      for(i in 1:nrow(Y)){
        if(Y[i,j]==0){b3[i,j]<- phi2+((1-phi2)*exp(-lambda2[i]))}
        else {b3[i,j]<- (1-phi2)*dpois(Y[i,j],lambda2[i])}
      }
    }
    
    #b3
    ############################################
    ############################################
    b22<- Pi*apply(b2,2,prod)
    
    b33<- (1-Pi)*apply(b3,2,prod)
    
    #b21<- Pi*apply(b2,2,sum)
    
    #b31<- (1-Pi)*apply(b3,2,sum)
    
    oldl<- l
    l<- sum(log(b22+b33))
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
  myoutput<- list(Pi=c(Pi,1-Pi),Phi1=c(phi1,1-phi1),Phi2=c(phi2,1-phi2),lambda=cbind(lambda1,
                                                                                 lambda2),loglik=l,all.loglik=ll[1:iter],
                  restart=0,ft="ZIPmixEM")
  
  class(myoutput)<- "mixEM"
  myoutput
  
}



#lambda1
out1<- ZIPmix2comp(Y,Phi1=0.25,Phi2=0.1,Pn=0.4,L1=lambda1,L2=lambda2,eps=1e-8,maxiter=1000)

#Phi1_vec <- rbind(Phi1_vec,out1$Phi1)
#Phi2_vec <- rbind(Phi2_vec,out1$Phi2)
#}

#boxplot(Phi1_vec[,1],ylim=c(0.22,0.26))
#abline(h=0.25)


#boxplot(Phi2_vec[,1],ylim=c(0.09,0.12))
#abline(h=0.1)

out1
L<- out1$lambda
L[,1]
lambda1

L[,2]
lambda2

