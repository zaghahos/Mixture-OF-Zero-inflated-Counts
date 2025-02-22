rm(list = ls())
#ls()

#rgamma(10,shape=20,scale =1/2 )

#rgamma(10,shape = 10,scale=1/0.2)


## Generate data from ZIP##
ZIP_r<- function(n,pi,theta){
  Z<- rbinom(size=1,n=n,p=pi)
  Y<- NULL
  for(i in 1:n){
    if(Z[i]==1) Y<- c(Y,0)
    else Y<- c(Y,rpois(1,lambda = theta))
  }
  return(Y)
}
########################################################
#clust<- rbinom(1000,size=1,p=0.35)
N<- 1000
G<- 10
r<- 500
#pn<- c(0.35,0.65)
#phi<- c(0.15,0.2)
phi1<- 0.15
phi2<- 0.4

lambda1<- rgamma(n=G,shape=20,scale=1/2)
#lambda1

lambda2<- rgamma(n=G,shape=10,scale=1/0.2)
#lambda2
p_hat<- NULL
phi_hat1<- NULL
phi_hat2<- NULL
lambda1_hat<- matrix(0,ncol = G,nrow=r)
lambda2_hat<- matrix(0,ncol = G,nrow=r)

for(iter in 1:r){
clust<- rbinom(size=1,n=N,p=0.7)

#length(clust[clust==1])


#clust<- sample(c(1,2),size=N,replace = TRUE,prob=pn) 
#clust


cell_matrix<- matrix(rep(0,N*G),nrow=G,ncol=N)

for(i in 1:N){
  for(j in 1:G){
    if(clust[i]==1) cell_matrix[j,i]<- ZIP_r(n=1,pi=phi1,theta=lambda1[j])
    else cell_matrix[j,i]<- ZIP_r(n=1,pi=phi2,theta=lambda2[j])
  }
}

cell_matrix[,1:10]
#cell_matrix
pz<- cbind(clust,1-clust)


p_hat<- c(p_hat,mean(pz[,1]))





#dim(cell_matrix)

U<- matrix(rep(0,N*G),ncol=N,nrow=G)
for(i in 1:N){
  for(j in 1:G){
    if(cell_matrix[j,i]==0) U[j,i]<- 1
    else U[j,i]<- 0
  }
}



mult_1<- matrix(0,ncol=N,nrow=G) 
mult_2<- matrix(0,ncol=N,nrow=G) 

for(i in 1:N){
  mult_1[,i]<- (pz[,1][i])*U[,i]
  mult_2[,i]<- (pz[,2][i])*U[,i]
}



phi_hat1<- c(phi_hat1,sum(mult_1)/(G*sum(pz[,1])))
phi_hat2<- c(phi_hat2,sum(mult_2)/(G*sum(pz[,2])))
#phi_hat1
#phi_hat2
#phi1
#phi2

D<- cell_matrix*(1-U)

num_1<- matrix(0,ncol=N,nrow=G)
num_2<- matrix(0,ncol=N,nrow=G)
denom_1<- matrix(0,ncol=N,nrow=G)
denom_2<- matrix(0,ncol=N,nrow=G)

for(i in 1:N){
  num_1[,i]<- D[,i]*pz[,1][i]
  num_2[,i]<- D[,i]*pz[,2][i]
  denom_1[,i]<- (1-U)[,i]*pz[,1][i]
  denom_2[,i]<- (1-U)[,i]*pz[,2][i]
}
#num_1[,1:10]
#num_2[,1:10]
#denom_1[,1:10]
#denom_2[,1:10]


lambda1_hat[iter,]<- (apply(num_1,1,sum))/(apply(denom_1,1,sum))
#lambda1_hat
#lambda1

lambda2_hat[iter,]<- (apply(num_2,1,sum))/(apply(denom_2,1,sum))
#lambda2_hat
#lambda2

}




par(mfrow=c(2,2))

boxplot(p_hat,main="p_hat")
abline(h=0.7)

boxplot(phi_hat1,main="phi1_hat")
abline(h=phi1)

boxplot(phi_hat2,main="phi2_hat")
abline(h=phi2)

#lambda1
#colMeans(lambda1_hat)

#lambda2
#colMeans(lambda2_hat)
print(list("Lambda1_hat", colMeans(lambda1_hat), "Lambda1",lambda1,
     "Lambda2_hat", colMeans(lambda2_hat), "Lambda2",lambda2))



cell_matrix
plot(hclust(dist(t(cell_matrix))))
hclust(dist(t(cell_matrix)))
###############################################################################################################
###############################################################################################################
rm(list = ls())

ZIP_r<- function(n,pi,theta){
  Z<- rbinom(size=1,n=n,p=pi)
  Y<- NULL
  for(i in 1:n){
    if(Z[i]==1) Y<- c(Y,0)
    else Y<- c(Y,rpois(1,lambda = theta))
  }
  return(Y)
}
########################################################
#clust<- rbinom(1000,size=1,p=0.35)
N<- 1000
G<- 10
r<- 500
#pn<- c(0.35,0.65)
#phi<- c(0.15,0.2)
phi1<- 0.15
phi2<- 0.4

lambda1<- rgamma(n=G,shape=20,scale=1/2)
#lambda1

lambda2<- rgamma(n=G,shape=10,scale=1/0.2)
#lambda2

p_hat<- NULL
phi_hat1<- NULL
phi_hat2<- NULL
lambda1_hat<- matrix(0,ncol = G,nrow=r)
lambda2_hat<- matrix(0,ncol = G,nrow=r)

for(iter in 1:r){

#p_hat<- NULL
#phi_hat1<- NULL
#phi_hat2<- NULL
#lambda1_hat<- matrix(0,ncol = G,nrow=r)
#lambda2_hat<- matrix(0,ncol = G,nrow=r)

#for(iter in 1:r){
  clust<- rbinom(size=1,n=N,p=0.7)
  
  #length(clust[clust==1])
  
  
  #clust<- sample(c(1,2),size=N,replace = TRUE,prob=pn) 
  #clust
  
  
  cell_matrix<- matrix(rep(0,N*G),nrow=G,ncol=N)
  
  for(i in 1:N){
    for(j in 1:G){
      if(clust[i]==1) cell_matrix[j,i]<- ZIP_r(n=1,pi=phi1,theta=lambda1[j])
      else cell_matrix[j,i]<- ZIP_r(n=1,pi=phi2,theta=lambda2[j])
    }
  }
  
#cell_matrix[,1:10]
#cell_matrix
###########################################################################################

posteriors<- function(Y,Pi,phi1,phi2,lambda1,lambda2){
  N<- ncol(Y)
  G<- nrow(Y)
  b2<- matrix(rep(0,N*G),nrow=G)
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){ 
      if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1[i]))}
      else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1[i])}
    }
  }
  b3<- matrix(rep(0,N*G),nrow=G)
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


#output1<- posteriors(Y=cell_matrix,Pi=0.7,phi1,phi2,lambda1,lambda2)
#output1$Pi_hat  

#posteriors(Y,Pi=0.4,phi1=0.25,phi2=0.1,lambda1=lambda1,lambda2=lambda2)
#################################################################################
#### U_kg####
IndUKG<- function(Y,Pi,phi1,phi2,lambda1,lambda2){
  N<- ncol(Y)
  G<- nrow(Y)
  ukg<- matrix(rep(0,G*N),nrow=G)
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){
      if(Y[i,j]==0){ukg[i,j]<- (Pi*phi1+(1-Pi)*phi2)/(Pi*(phi1+(1-phi1)*exp(-lambda1[i]))+
                                                        (1-Pi)*(phi2+(1-phi2)*exp(-lambda2[i])))}
      else {ukg[i,j]<- 0}
    }
    
  }
  return(ukg)
}

#outpu2<- IndUKG(Y=cell_matrix,Pi=0.7,phi1,phi2,lambda1,lambda2)
#outpu2[,1:10]
#output1$postprobs[1:10]
#output1$ZZ
############################################################################################
#### Phi_hats###
phi_hats<- function(Y,postprobs,ukg){
  G<- nrow(Y)
  mult1<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
  mult2<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
  
  for(i in 1:length(postprobs)){
    mult1[,i]<- postprobs[i]*ukg[,i]
    mult2[,i]<- (1-postprobs[i])*ukg[,i]
  }
  phi_hat1<- sum(mult1)/(G*sum(postprobs))
  phi_hat2<- sum(mult2)/(G*sum(1-postprobs))
  phi_hat<- c(phi_hat1,phi_hat2)
  names(phi_hat)<- c("Phi1_hat","Phi2_hat")
  return(phi_hat)
  
}  

#phi_hats(Y=cell_matrix,postprobs=posteriors(Y=cell_matrix,Pi=0.7,phi1,phi2,lambda1=lambda1,lambda2=lambda2)$postprobs

 #       ,ukg=IndUKG(Y=cell_matrix,Pi=0.7,phi1,phi2,lambda1=lambda1,lambda2=lambda2))


#################################################################################################
lambda_hats<- function(Y,ukg,postprobs){
  mult3<- Y*(1-ukg)
  mult4<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
  mult5<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
  mult6<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
  mult7<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
  
  for(i in 1:length(postprobs)){
    mult4[,i]<- postprobs[i]*mult3[,i]
    mult5[,i]<- (1-postprobs[i])*mult3[,i]
    mult6[,i]<- postprobs[i]*(1-ukg)[,i]
    mult7[,i]<- (1-postprobs[i])*(1-ukg)[,i]
  }
  

  lambda1_hat<- apply(mult4,1,sum)/apply(mult6,1,sum)
  lambda2_hat<- apply(mult5,1,sum)/apply(mult7,1,sum)
  
  #lambda_hat<- list("lambda1_hat",lambda1_hat,"lambda1",lambda1,"lambda2_hat",lambda2_hat,"lambda2",lambda2)
  lambda_hat<- cbind(lambda1_hat,lambda2_hat)
  colnames(lambda_hat)<- c("lambda1_hat","lambda2_hat")
  return(lambda_hat)
  
}

#lambda_hats(Y=cell_matrix,ukg=IndUKG(Y=cell_matrix,Pi=0.7,phi1,phi2,lambda1=lambda1,lambda2=lambda2)
 #            ,postprobs=posteriors(Y=cell_matrix,Pi=0.7,phi1,phi2,lambda1=lambda1,lambda2=lambda2)$postprobs)

#cbind(lambda1,lambda2)
######################################################################################################################
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

Y<- cell_matrix

myout1<- EMZIPmixture2com(Y,Pi_init=0.65,Phi1_init=0.12,Phi2_init=0.37,lambda1_init=lambda1,
                          lambda2_init=lambda2,eps=1e-4,maxiter=1000)


p_hat<- c(p_hat,myout1$Pi_hat)
phi_hat1<- c(phi_hat1,myout1$Phi1_hats[1])
phi_hat2<- c(phi_hat2,myout1$Phi2_hats[1])
lambda1_hat[iter,]<- myout1$lambda_hats[,1]
lambda2_hat[iter,]<- myout1$lambda_hats[,2]

}

p_hat
phi_hat1
phi_hat2
lambda1_hat
lambda2_hat



par(mfrow=c(2,2))

boxplot(p_hat,main="p_hat")
abline(h=0.7)

boxplot(phi_hat1,main="phi1_hat")
abline(h=phi1)

boxplot(phi_hat2,main="phi2_hat")
abline(h=phi2)

#lambda1
#colMeans(lambda1_hat)

#lambda2
#colMeans(lambda2_hat)
print(list("Lambda1_hat", colMeans(lambda1_hat), "Lambda1",lambda1,
           "Lambda2_hat", colMeans(lambda2_hat), "Lambda2",lambda2))

#myout1
#myout1$lambda_hats[,1]
#lambda1
#myout1$lambda_hats[,2]
#lambda2
################################################################################################################
################################################################################################################
rm(list = ls())

ZIP_r<- function(n,pi,theta){
  Z<- rbinom(size=1,n=n,p=pi)
  Y<- NULL
  for(i in 1:n){
    if(Z[i]==1) Y<- c(Y,0)
    else Y<- c(Y,rpois(1,lambda = theta))
  }
  return(Y)
}
########################################################
#clust<- rbinom(1000,size=1,p=0.35)
N<- 3000
G<- 1000
r<- 2
#pn<- c(0.35,0.65)
#phi<- c(0.15,0.2)
phi1<- 0.15
phi2<- 0.4

lambda1<- rgamma(n=G,shape=20,scale=1/2)
#lambda1

lambda2<- rgamma(n=G,shape=10,scale=1/0.2)
#lambda2

p_hat<- NULL
phi_hat1<- NULL
phi_hat2<- NULL
lambda1_hat<- matrix(0,ncol = G,nrow=r)
lambda2_hat<- matrix(0,ncol = G,nrow=r)

for(iteration in 1:r){
  
  #p_hat<- NULL
  #phi_hat1<- NULL
  #phi_hat2<- NULL
  #lambda1_hat<- matrix(0,ncol = G,nrow=r)
  #lambda2_hat<- matrix(0,ncol = G,nrow=r)
  
  #for(iter in 1:r){
  clust<- rbinom(size=1,n=N,p=0.7)
  
  #length(clust[clust==1])
  
  
  #clust<- sample(c(1,2),size=N,replace = TRUE,prob=pn) 
  #clust
  
  
  cell_matrix<- matrix(rep(0,N*G),nrow=G,ncol=N)
  
  for(i in 1:N){
    for(j in 1:G){
      if(clust[i]==1) cell_matrix[j,i]<- ZIP_r(n=1,pi=phi1,theta=lambda1[j])
      else cell_matrix[j,i]<- ZIP_r(n=1,pi=phi2,theta=lambda2[j])
    }
  }
  
  #cell_matrix[,1:10]
  #cell_matrix






ZIP_mix2comp<- function(Y,Pi,phi1,phi2,lambda1,lambda2,eps=1e-8, maxit=1000){
  dl <- 1 + eps
  iter<-0
  ll <- rep(0, maxit+1)
  N<- ncol(Y)
  G<- nrow(Y)
  
  b2<- matrix(rep(0,N*G),nrow=G)
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){ 
      if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1[i]))}
      else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1[i])}
    }
  }
  b3<- matrix(rep(0,N*G),nrow=G)
  for(j in 1:ncol(Y)){
    for(i in 1:nrow(Y)){
      if(Y[i,j]==0){b3[i,j]<- phi2+((1-phi2)*exp(-lambda2[i]))}
      else {b3[i,j]<- (1-phi2)*dpois(Y[i,j],lambda2[i])}
    }
  }
  b22<- Pi*apply(b2,2,prod)
  b33<- (1-Pi)*apply(b3,2,prod)
  b44<- (Pi*b2)+((1-Pi)*b3)
  #print(length(b22/(b22+b33)))
  
  l<- sum(sum(log(b44)))
  print(l)
  #return(l)
  while (abs(dl)>eps && iter<maxit) { ## Here I use abs(dl) rather than dl##
    iter<-iter+1
    ll[iter] <- l
    postprobs<- b22/(b22+b33)
    Pi<- mean(postprobs)
    ukg<- matrix(rep(0,G*N),nrow=G)
      for(j in 1:ncol(Y)){
        for(i in 1:nrow(Y)){
          if(Y[i,j]==0){ukg[i,j]<- (Pi*phi1+(1-Pi)*phi2)/(Pi*(phi1+(1-phi1)*exp(-lambda1[i]))+
                                                            (1-Pi)*(phi2+(1-phi2)*exp(-lambda2[i])))}
          else {ukg[i,j]<- 0}
        }
        
      }
    mult1<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
    mult2<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
    for(i in 1:length(postprobs)){
      mult1[,i]<- postprobs[i]*ukg[,i]
      mult2[,i]<- (1-postprobs[i])*ukg[,i]}
    phi1<- sum(mult1)/(G*sum(postprobs))
    phi2<- sum(mult2)/(G*sum(1-postprobs))
    
    mult3<- Y*(1-ukg)
    mult4<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
    mult5<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
    mult6<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
    mult7<- matrix(0,ncol=ncol(Y),nrow=nrow(Y))
    
    for(i in 1:length(postprobs)){
      mult4[,i]<- postprobs[i]*mult3[,i]
      mult5[,i]<- (1-postprobs[i])*mult3[,i]
      mult6[,i]<- postprobs[i]*(1-ukg)[,i]
      mult7[,i]<- (1-postprobs[i])*(1-ukg)[,i]
    }
    
    
    lambda1<- apply(mult4,1,sum)/apply(mult6,1,sum)
    lambda2<- apply(mult5,1,sum)/apply(mult7,1,sum)
    
    ################################################################
    #b2<- matrix(rep(0,N*G),nrow=G)
    for(j in 1:ncol(Y)){
      for(i in 1:nrow(Y)){ 
        if(Y[i,j]==0){b2[i,j]<- phi1+((1-phi1)*exp(-lambda1[i]))}
        else {b2[i,j]<- (1-phi1)*dpois(Y[i,j],lambda1[i])}
      }
    }
    #b3<- matrix(rep(0,N*G),nrow=G)
    for(j in 1:ncol(Y)){
      for(i in 1:nrow(Y)){
        if(Y[i,j]==0){b3[i,j]<- phi2+((1-phi2)*exp(-lambda2[i]))}
        else {b3[i,j]<- (1-phi2)*dpois(Y[i,j],lambda2[i])}
      }
    }
    b22<- Pi*apply(b2,2,prod)
    b33<- (1-Pi)*apply(b3,2,prod)
    b44<- (Pi*b2)+((1-Pi)*b3)
    print(b44[1:10,1:10])
    
    oldl<- l
    
    
    print(log(b44[1:10,1:10]))
    l<- sum(sum(log(b44)))
    print(l)
    dl<- l-oldl
  }
    #cat("number of iterations=", iter, "\n")
    iter <- iter+1
    ll[iter] <- l
    postprobs <- cbind(postprobs, 1-postprobs)
    colnames(postprobs) <- c(paste("comp", ".", 1:2, sep = ""))
    
    
    out <- list(Pi=Pi,phi1=phi1,phi2=phi2,lambda1=lambda1,lambda2=lambda2,
                loglik=l,all.loglik=ll[1:iter],restart=0,ft="ZIPmix")
    class(out)<- "mixEM"
      out

}
    
  

output<- ZIP_mix2comp(Y=cell_matrix,Pi=0.65,phi1=0.12,phi2=0.37,lambda1,lambda2,eps=1e-8, maxit=1000)
  
p_hat<- c(p_hat,output$Pi)
phi_hat1<- c(phi_hat1,output$phi1)
phi_hat2<- c(phi_hat2,output$phi2)
lambda1_hat[iteration,]<- output$lambda1
lambda2_hat[iteration,]<- output$lambda2

}

p_hat
phi_hat1
phi_hat2
lambda1_hat
lambda2_hat

lambda1
lambda2
par(mfrow=c(2,2))

boxplot(p_hat,main="p_hat")
abline(h=0.7)

boxplot(phi_hat1,main="phi1_hat")
abline(h=phi1)

boxplot(phi_hat2,main="phi2_hat")
abline(h=phi2)

#lambda1
#colMeans(lambda1_hat)

#lambda2
#colMeans(lambda2_hat)
print(list("Lambda1_hat", colMeans(lambda1_hat), "Lambda1",lambda1,
           "Lambda2_hat", colMeans(lambda2_hat), "Lambda2",lambda2))
#}

