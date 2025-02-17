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
M<- 100
N<- 1000  
phi1<- c(0.25,0.75)
phi2<- c(0.1,0.9)
Pi<- c(0.4,0.6)
lambda1<- rgamma(n=M,shape=0.5,rate=1/2)
#lambda2<- rgamma(n=M,shape = 1.5,rate = 1/2)
lambda2<- rgamma(n=M,shape = 2,rate = 1/2)

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


### 1.Generate 2 vectors of Z For the Clusters####################
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

#Y[,1]
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

###################################################################
#### Posterior Pi & Z_nk####
posteriors<- function(Y,Pi,phi1,phi2,lambda1,lambda2){
  N<- ncol(Y)
  M<- nrow(Y)
  b2<- matrix(rep(0,N*M),nrow=M)
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){ 
      if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1[i]))}
      else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1[i])}
    }
  }
  b3<- matrix(rep(0,N*M),nrow=M)
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){
      if(Y[i,j]==0){b3[i,j]<- phi2+((1-phi2)*exp(-lambda2[i]))}
      else {b3[i,j]<- (1-phi2)*dpois(Y[i,j],lambda2[i])}
    }
  }
  b22<- Pi*apply(b2,2,prod)
  
  b33<- (1-Pi)*apply(b3,2,prod)
  postprobs<- b22/(b22+b33)
  Pi<- mean(postprobs)
  ZZ<- cbind(postprobs,1-postprobs)
  out1<- list(Pi,postprobs,ZZ)
  names(out1)<- c("Pi_hat","postprobs","ZZ")
  return(out1)
}
  
  
#posteriors(Y,Pi=0.4,phi1=0.25,phi2=0.1,lambda1=lambda1,lambda2=lambda2)
#################################################################################
#### U_kg####
IndUKG<- function(Y,Pi,phi1,phi2,lambda1,lambda2){
  N<- ncol(Y)
  M<- nrow(Y)
  ukg<- matrix(rep(0,M*N),nrow=M)
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){
      if(Y[i,j]==0){ukg[i,j]<- (Pi*phi1+(1-Pi)*phi2)/(Pi*(phi1+(1-phi1)*exp(-lambda1[i]))+
                                                        (1-Pi)*(phi2+(1-phi2)*exp(-lambda2[i])))}
      else {ukg[i,j]<- 0}
    }
    
  }
  return(ukg)
}

#IndUKG(Y,Pi=0.4,phi1=0.25,phi2=0.1,lambda1=lambda1,lambda2=lambda2)
############################################################################################
#### Phi_hats###
phi_hats<- function(Y,postprobs,ukg){
  M<- nrow(Y)
  comb1<- postprobs*t(ukg)
  phi_hat1<- sum(comb1)/(M*sum(postprobs))
  #phi_hat1
  comb2<- (1-postprobs)*t(ukg)
  phi_hat2<- sum(comb2)/(M*sum(1-postprobs))
  #phi_hat2
  return(c(phi_hat1,phi_hat2))
}  
  
#phi_hats(Y,postprobs=posteriors(Y,Pi=0.4,phi1=0.25,phi2=0.1,lambda1=lambda1,lambda2=lambda2)$postprobs

 #        ,ukg=IndUKG(Y,Pi=0.4,phi1=0.25,phi2=0.1,lambda1=lambda1,lambda2=lambda2))
  

#################################################################################################
#### Lambda_hats####
lambda_hats<- function(Y,ukg,postprobs){
  d1<- (1-ukg)*Y
  d12<- apply((postprobs*t(d1)),2,sum)
  d13<- apply((postprobs*t(1-ukg)),2,sum)
  
  d22<- apply(((1-postprobs)*t(d1)),2,sum)
  d23<- apply(((1-postprobs)*t(1-ukg)),2,sum)
  for (i in 1:length(d13)){
    if (d13[i]==0) {d13[i]<- 0.000001}
    if (d23[i]==0) {d23[i]<- 0.000001}
  }
  
  
  lambda1<- d12/d13
  lambda2<- d22/d23
  
    
  
  return(cbind(lambda1,lambda2))
}

#lambda_hats(Y,ukg=IndUKG(Y,Pi=0.4,phi1=0.25,phi2=0.1,lambda1=lambda1,lambda2=lambda2)
 #             ,postprobs=posteriors(Y,Pi=0.4,phi1=0.25,phi2=0.1,lambda1=lambda1,lambda2=lambda2)$postprobs)
  
#cbind(lambda1,lambda2)
#################################################################################
#################################################################################
#### E& M step: Iteration ####
EMZIPmixture2com<- function(Y,Pi_init,Phi1_init,Phi2_init,lambda1_init,lambda2_init,eps=1e-4,maxiter=1000){
  flag<- 0
  Pi_cur<- Pi_init
  Phi1_cur<- Phi1_init
  Phi2_cur<- Phi2_init
  lambda1_cur<- lambda1_init
  lambda2_cur<- lambda2_init
  iter<- NULL
  for (i in 1:maxiter){
    Pi_cur<- Pi_cur
    Phi1_cur<- Phi1_cur
    Phi2_cur<- Phi2_cur
    lambda1_cur<- lambda1_cur
    lambda2_cur<- lambda2_cur
    ### Pi_hat###
    posteriornew<- posteriors(Y,Pi=Pi_cur,phi1=Phi1_cur,phi2=Phi2_cur,lambda1=lambda1_cur,lambda2=lambda2_cur)
    
    Pi_new<- posteriornew$Pi_hat
    ### Phi_hat###
    Phi_hatnew<- phi_hats(Y,posteriornew$postprobs,
                          ukg= IndUKG(Y,Pi=Pi_cur,phi1=Phi1_cur,phi2=Phi2_cur,lambda1=lambda1_cur,
                                      lambda2=lambda2_cur))
    Phi1_new<- Phi_hatnew[1]
    Phi2_new<- Phi_hatnew[2]
    ### lambda_hat###
    lambda_new<- lambda_hats(Y,ukg=IndUKG(Y,Pi=Pi_cur,phi1=Phi1_cur,phi2=Phi2_cur,lambda1=lambda1_cur,
                                          lambda2=lambda2_cur)
                  ,postprobs=posteriornew$postprobs)
    lambda1_new<- lambda_new[,1]
    lambda2_new<- lambda_new[,2]
    ############################################################
    diff1<- abs(Pi_new-Pi_cur)
    diff2<- abs(Phi1_new-Phi1_cur)
    diff3<- abs(Phi2_new-Phi2_cur)
    diff4<- mean(abs(lambda1_new-lambda1_cur))
    diff5<- mean(abs(lambda2_new-lambda2_cur))
    if ((diff1<eps) & (diff2<eps) & (diff3<eps) & (diff4<eps) & ( diff5<eps)){flag<- 1; break}
    Pi_cur<- Pi_new
    Phi1_cur<- Phi1_new
    Phi2_cur<- Phi2_new
    lambda1_cur<- lambda1_new
    lambda2_cur<- lambda2_new
  }
  if (!flag) warning("Didn't converge \n")
  update_estimates<- list(Pi_cur,c(Phi1_cur,1-Phi1_cur),c(Phi2_cur,1-Phi2_cur),cbind(lambda1_cur,lambda2_cur),i)
  names(update_estimates)<- c("Pi_hat","Phi1_hats","Phi2_hats","lambda_hats","iteration")
  return(update_estimates)
  #}
  
}


myout1<- EMZIPmixture2com(Y,Pi_init=0.4,Phi1_init=0.25,Phi2_init=0.1,lambda1_init=lambda1,
                          lambda2_init=lambda2,eps=1e-4,maxiter=1000)
  


myout1
myout1$lambda_hats[,1]
lambda1
myout1$lambda_hats[,2]
lambda2
