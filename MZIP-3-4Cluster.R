rm(list = ls())

###################################################################################################
ZIP_mix_K_comp<- function(Y,P,Phi,lambda,eps=1e-8, maxit=1000){
    dl <- 1 + eps
    iter<-0
    ll <- rep(0, maxit+1)
    G<- ncol(Y)
    N<- nrow(Y)
    K<- length(Phi)
    
    bk<- array(NA,c(N,G,K))
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0) bk[n,g,k]<- Phi[k]+((1-Phi[k])*exp(-lambda[k,g]))
          else bk[n,g,k]<- (1-Phi[k])*dpois(Y[n,g],lambda[k,g])
        }
      }
    }
    
    #bk_pord<- matrix(NA,nrow=N,ncol=K)
    #for(k in 1:K)bk_pord[,k]<- P[k]*apply(bk[,,k],1,prod)
    
    b<- array(NA,c(N,G,K))
    
    for(k in 1:K) b[,,k]<- P[k]*bk[,,k]
    bp_sum<- apply(b,c(1,2),sum)
    #print(bp_sum)
    l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
    
    Z_nk<- matrix(NA,ncol=K,nrow=N)
    thresh<- -744
    
    
    while (abs(dl)>eps && iter<maxit) { ## Here I use abs(dl) rather than dl##
      iter<-iter+1
      ll[iter] <- l
     
      
      
      Z_nk<- matrix(data=log(P),ncol=K,nrow=N,byrow=T)
      
      for(k in 1:K){
        for(g in 1:G){
          Z_nk[,k]<- Z_nk[,k]+(log(bk[,g,k]))
        }
      }
      
      if(K>1){
        v1<- which(apply(Z_nk,1,max)< thresh)
        v3<- 1:N
        len<- length(v1)
        #print(len)
        if(len>0){
          v2<- apply(array(Z_nk[v1,],dim=c(len,K)),1,order)[K,]
          ddd<- cbind(v1,v2)
          #print(ddd)
          Z_nk[v1,]<- 0
          Z_nk[ddd]<- 1
          v3<- -v1
        }
        
        Z_nk[v3,]<- exp(Z_nk[v3,])
      }
      Z_nk<- Z_nk/apply(Z_nk,1,sum)
      #print(Z_nk)
      
      epsilon<- 1e-10
      sl<- length(Z_nk[Z_nk<epsilon])
      bl<- length(Z_nk[Z_nk>1-epsilon])
      Z_nk[Z_nk<epsilon]<- rep(epsilon,sl)
      Z_nk[Z_nk>1-epsilon]<- rep(1-epsilon,bl)
      Z_nk<- Z_nk/apply(Z_nk,1,sum)
      #print(Z_nk)
      
      P<- apply(Z_nk,2,mean)
      
      #for(k in 1:K) Z_nk[,k]<- bk_pord[,k]/(rowSums(bk_pord))
      #P<- apply(Z_nk,2,mean)
      
      #print(P)
      
    
      ukg<- matrix(rep(NA,G*N),ncol=G)
      denom_U<- NULL
      for(k in 1:K){
        denom_U<- c(denom_U,P[k]*(Phi[k]+((1-Phi[k])*(exp(-lambda[k,g])))))
        
      }
      
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0){ukg[n,g]<- sum(P*Phi)/sum(denom_U)}
          else {ukg[n,g]<- 0}
        }
        
      }
      
      
      mult_Z_u<- array(NA,c(N,G,K))
      mult_Z_oneminus_U<- array(NA,c(N,G,K))
      mult_Z_oneminus_U_Y<- array(NA,c(N,G,K))
      
      
      
      for(k in 1:K){
        mult_Z_u[,,k]<- Z_nk[,k]*ukg
        mult_Z_oneminus_U[,,k]<- Z_nk[,k]*(1-ukg)
        mult_Z_oneminus_U_Y[,,k]<- Z_nk[,k]*(1-ukg)*Y
        
      }
      
      
      Phi<- apply(mult_Z_u,3,sum)/(G*apply(Z_nk,2,sum))
      #print(Phi)
      lambda<- t(apply(mult_Z_oneminus_U_Y,c(2,3),sum)/apply(mult_Z_oneminus_U,c(2,3),sum))
      
      
      for(k in 1:K){
        for(n in 1:N){
          for(g in 1:G){
            if(Y[n,g]==0) bk[n,g,k]<- Phi[k]+((1-Phi[k])*exp(-lambda[k,g]))
            else bk[n,g,k]<- (1-Phi[k])*dpois(Y[n,g],lambda[k,g])
          }
        }
      }
      
      #bk_pord<- matrix(NA,nrow=N,ncol=K)
      #for(k in 1:K)bk_pord[,k]<- P[k]*apply(bk[,,k],1,prod)
      
      #b<- array(NA,c(N,G,K))
      
      #b<- array(NA,c(N,G,K))
      for(k in 1:K) b[,,k]<- P[k]*bk[,,k]
      bp_sum<- apply(b,c(1,2),sum)
      
      oldl<- l
      
      l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
      
      ################################################################
     
      #print(l)
      dl<- l-oldl
    }
    cat("number of iterations=", iter, "\n")
    iter <- iter+1
    ll[iter] <- l
    postprobs <- Z_nk
    colnames(postprobs) <- c(paste("comp", ".", 1:4, sep = ""))
    
    
    out <- list(P=P,Phi=Phi,lambda=lambda,postprobs=postprobs,
                loglik=l,all.loglik=ll[1:iter],restart=0,ft="ZIPmix")
    class(out)<- "mixEM"
    out
    
  }
  
#########################################################################################################################
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
#N<- 1000
#G<- 100
N<- 1000
G<- 10
Phi<- c(0.1,0.1,0.2,0.2)
P<- c(0.25,0.25,0.25,0.25)
P1<- c(0.4,0.3,0.2,0.1)
#lambda1<- rep(4,G)
#lambda2<- rep(8,G)
#lambda3<- rep(15,G)
#lambda4<- rep(25,G)
lambda1<- rgamma(n=G,shape=10,scale=1/2)
#lambda1

lambda2<- rgamma(n=G,shape=10,scale=1/0.2)
#lambda2

lambda3<- rgamma(n=G,shape=20,scale=1/0.2)
#lambda3

lambda4<- rgamma(n=G,shape=20,scale=1/0.8)
#lambda4

lambda<- rbind(lambda1,lambda2,lambda3,lambda4)
#lambda


S<- 100
#S<- 3
#P_hat<- matrix(NA,ncol=length(P),nrow=S)
P_hat<- NULL
Phi_hat<- matrix(NA,ncol=length(P),nrow=S)
#Phi_hat
#phi_hat2<- NULL
lambda1_hat<- matrix(NA,ncol = G,nrow=S)
lambda2_hat<- matrix(NA,ncol = G,nrow=S)
lambda3_hat<- matrix(NA,ncol = G,nrow=S)
lambda4_hat<- matrix(NA,ncol = G,nrow=S)


#for(s in 1:S){




for(s in 1:S){


#for(iteration in 1:r){

#p_hat<- NULL
#phi_hat1<- NULL
#phi_hat2<- NULL
#lambda1_hat<- matrix(0,ncol = G,nrow=r)
#lambda2_hat<- matrix(0,ncol = G,nrow=r)

#for(iter in 1:r){
#clust<- rbinom(size=1,n=N,p=0.5)
clust<- sample(1:4,size=N,replace = TRUE,prob=P)
#clust
#length(clust[clust==1])


#clust<- sample(c(1,2),size=N,replace = TRUE,prob=pn) 
#clust


cell_matrix<- matrix(rep(0,N*G),nrow=N,ncol=G)

for(n in 1:N){
  for(g in 1:G){
    if(clust[n]==1) cell_matrix[n,g]<- ZIP_r(n=1,pi=Phi[1],theta=lambda[1,g])
    else if(clust[n]==2) cell_matrix[n,g]<- ZIP_r(n=1,pi=Phi[2],theta=lambda[2,g])
    else if(clust[n]==3) cell_matrix[n,g]<- ZIP_r(n=1,pi=Phi[3],theta=lambda[3,g])
    else cell_matrix[n,g]<- ZIP_r(n=1,pi=Phi[4],theta=lambda[4,g])
  }
}

#cell_matrix
#cell_matrix

#plot(hclust(dist(cell_matrix)))


#cell_matrix
#dim(cell_matrix)  

#########################################################################################################################
Y<- cell_matrix  
#Phi<- c(0.1,0.1,0.1,0.1)
#lambda<- rbind(lambda1,lambda2,lambda3,lambda4)
#lambda
#P<- c(0.4,0.3,0.2,0.1)

output_4comp<- ZIP_mix_K_comp(Y,P,Phi,lambda,eps=1e-8, maxit=1000)
#output_4comp  

  
#output_4comp
#output_4comp$P
#P
#output_4comp$Phi
#Phi
#output_4comp$lambda[,1:5]
#lambda[,1:5]
#####################################################################################
P_hat<- rbind(P_hat,output_4comp$P)
Phi_hat[s,]<- output_4comp$Phi
lambda1_hat[s,]<- (output_4comp$lambda)[1,]
lambda2_hat[s,]<- (output_4comp$lambda)[2,]
lambda3_hat[s,]<- (output_4comp$lambda)[3,]
lambda4_hat[s,]<- (output_4comp$lambda)[4,]

}


P_hat
Phi_hat
#lambda1_hat
#lambda2_hat
#lambda3_hat
#lambda4_hat


##################################################################################################################
pdf("Probabilities-MZIP.pdf")

par(mfrow=c(1,2))
# 
boxplot(P_hat,main="P_hat",ylim=c(0.2,0.3))
abline(h=0.25)
 
boxplot(Phi_hat,main="Phi_hat",ylim=c(0.05,0.25))
abline(h=0.1)
#abline(h=0.4)
abline(h=0.2)
#abline(0.3)

#boxplot(lambda1_hat,main="lambda1_hat")
#boxplot(lambda2_hat,main="lambda2_hat")
#boxplot(lambda3_hat,main="lambda3_hat")
#boxplot(lambda4_hat,main="lambda4_hat")
dev.off()

print(list("Lambda1_hat", colMeans(lambda1_hat)[1:10], "Lambda1",lambda1[1:10],
            "Lambda2_hat", colMeans(lambda2_hat)[1:10], "Lambda2",lambda2[1:10],
              "Lambda3_hat", colMeans(lambda3_hat)[1:10], "Lambda3",lambda3[1:10],  
              "Lambda4_hat", colMeans(lambda4_hat)[1:10], "Lambda4",lambda4[1:10] ))
# #}
# 
#colMeans(lambda1_hat)
L1<- cbind(colMeans(lambda1_hat),lambda1)
colnames(L1)<- c("Lambda1_hat","Lambda1")
L2<- cbind(colMeans(lambda2_hat),lambda2)
colnames(L2)<- c("Lambda2_hat","Lambda2")
L3<- cbind(colMeans(lambda3_hat),lambda3)
colnames(L3)<- c("Lambda3_hat","Lambda3")
L4<- cbind(colMeans(lambda4_hat),lambda4)
colnames(L4)<- c("Lambda4_hat","Lambda4")

L1
L2
L3
L4
Lambdas<- cbind(L1,L2,L3,L4)
library(xtable)
xtable(Lambdas)
