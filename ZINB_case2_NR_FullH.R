### True Z and U and Mu , Update nu###
rm(list = ls())
library(FDRSeg)
library(MASS)
library(pheatmap)

###################################################################################################
### Generate data from Zero-inflated Negative Binomial Distribution### 
ZINB_r_U<- function(n,prob,theta,disp_param){
  Z<- rbinom(size=1,n=n,p=prob)
  Y<- NULL
  for(i in 1:n){
    if(Z[i]==1) Y<- c(Y,0)
    else Y<- c(Y,rnbinom(1,mu = theta,size=disp_param))
  }
  return(output=list(Y=Y,U=Z))
}
#ZINB_r_U(50,0.5,4,10)

#################################################################################
#################################################################################################################################################
NR_ZINB<- function(Y,totals,beta_0,rho,Z_nk,U_ng,overdisp_para,criterion=1e-5,maxiter=500){
  N<- nrow(Y)
  G<- ncol(Y)
  K<- dim(Z_nk)[2]
  nu<- overdisp_para
  
  a<- 1/nu
  ### Initialization###
  mu_ng<- array(data=NA,c(N,G,K))
  diff_Y_mu<- array(data=NA,c(N,G,K))
  U_diff_Y_mu<- array(data=NA,c(N,G,K))
  U_mu<- array(data=NA,c(N,G,K))
  
  Z_diff<- array(data=NA,c(N,G,K))
  Z_mu<- array(data=NA,c(N,G,K))
  theta<- numeric((G*(K+1)))
  
  ### Assign values to theta vector###
  theta[1:G]<- beta_0
  
  for(k in 1:K){
    theta[((k*G)+1):((k*G)+G)]<- rho[k,]
  }
  
  it<- 0
  diff_theta<- 1000
  while(it<maxiter & diff_theta>criterion){
    #while(it<maxiter & all(diff_theta>criterion)){
    
    it<- it+1
    #print(it)
    ### mu###
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          mu_ng[n,g,k]<- exp(log(totals[n])+beta_0[g]+rho[k,g])
          
        }
      }
    }
    #mu_ng
    ### Y-mu###
    for(k in 1:K){
      diff_Y_mu[,,k]<- (Y-mu_ng[,,k])/(1+(a[k]*mu_ng[,,k]))
    }
    
    for(k in 1:K){
      U_diff_Y_mu[,,k]<- (1-U_ng[,,k])*diff_Y_mu[,,k]
      U_mu[,,k]<- (1-U_ng[,,k])*(mu_ng[,,k]*(1+a[k]*Y))/((1+a[k]*mu_ng[,,k])^2)
    }
    
    U_diff_Y_mu
    U_mu
    
    
    ### Z*(Y-mu)& Z*mu###
    for(k in 1:K){
      Z_diff[,,k]<- Z_nk[,k]*U_diff_Y_mu[,,k]
      Z_mu[,,k]<- Z_nk[,k]*U_mu[,,k]
    }
    
    #Z_diff
    #Z_mu
    
    ### derivatives of beta_0###
    db_0<- ddb_0<-  numeric(G)
    
    for(g in 1:G){
      db_0[g]<- sum(Z_diff[,g,])
      ddb_0[g]<- sum(Z_mu[,g,])
    }
    
    ### Derivatives of rho's###
    
    drho<- ddrho<- matrix(data=NA,nrow=K,ncol=G)
    
    for(k in 1:K){
      for(g in 1:G){
        drho[k,g]<- sum(Z_diff[,g,k])
        ddrho[k,g]<- sum(Z_mu[,g,k])
      }
    }
    ### Gradient Vecotr###
    grad<- numeric((G*(K+1)))
    
    grad[1:G]<- db_0
    for(k in 1:K){
      grad[((k*G)+1):((k*G)+G)]<- drho[k,]
    }
    ### Hessian Matrix###
    hessian<- matrix(0,ncol=(G*(K+1)),nrow=(G*(K+1)))
    hessian[1:G,1:G]<- hessian[1:G,1:G]-diag(ddb_0)
    
    for(k in 1:K){
      hessian[((k*G)+1):((k*G)+G),((k*G)+1):((k*G)+G)]<- hessian[((k*G)+1):((k*G)+G),((k*G)+1):((k*G)+G)]-diag(ddrho[k,])
      hessian[((k*G)+1):((k*G)+G),1:G]<- hessian[((k*G)+1):((k*G)+G),1:G]-diag(ddrho[k,]) 
      hessian[1:G,((k*G)+1):((k*G)+G)]<- hessian[1:G,((k*G)+1):((k*G)+G)]-diag(ddrho[k,])
    }
    
    ### Updating theta using Newton-Raphson Algorithm###
    theta_new<- theta-qr.coef(qr(hessian,tol=1e-300), grad) ### tol=1e-300
    
    for(i in 1:length(theta_new)){
      if(is.na(theta_new[i])) theta_new[i]<- theta[i]
      
    }
    
    diff_theta<- mean(abs(theta_new-theta))
    #diff_theta<- abs(theta_new-theta)
    
    #print(theta)
    #print(theta_new)
    
    theta<- theta_new
    
    #beta_0<- theta_new[1:G]
    #for(k in 1:K){
    # rho[k,]<- theta_new[((k*G)+1):((k*G)+G)]
    #}
    
    beta_0<- theta[1:G]
    for(k in 1:K){
      rho[k,]<- theta[((k*G)+1):((k*G)+G)]
    }
    
    
    
    #out<- list(beta_0=beta_0,rho=rho,thtea=theta,iteration=it)
    
    
  }
  
  out<- list(beta_0=beta_0,rho=rho,thtea=theta,iteration=it)
  
  return(out)
  
  
  #  }
  
  
}


##############################################################################################
K=2
N=100
G=10
Phi<- c(0.1,0.2)
P<- c(0.5,0.5)
P1<- c(0.7,0.3)

K<- length(P)

size<- c(5,20)

S<- 20
#S<- 50

beta_0_est<- matrix(NA,nrow=S,ncol=G)
rho1_est<- matrix(NA,nrow=S,ncol=G)
rho2_est<- matrix(NA,nrow=S,ncol=G)
iter<- NULL


for(s in 1:S){
  
  
  clust <- sample(c(1:K),size=N,replace = TRUE,prob=P)
  
  
  ### creating lambdas 
  
  T <- round(rnorm(N,mean=10,sd=0.5)) ## vector with known Totals#hist(T)
  #T
  beta_0 <- rep(0.85,G) ## vector of beta_0 for all genes (baseline)
  
  
  
  rho_1 <- c(rep(2,round(G/2)),rep(-2,round(G/2))) ## fixed effect due to cluster 1, no effect for the last 5 genes
  rho_2 <- c(rep(-2,round(G/2)),rep(2,round(G/2))) ## fixed effect due to cluster 2, no effect for the first 5 genes
  rho <- rbind(rho_1,rho_2)
  
  
  cell_matrix<- U_ng<-  matrix(rep(0,N*G),nrow=N,ncol=G)
  
  for(n in 1:N){
    for(g in 1:G){
      if(clust[n]==1) {
        out_1 = ZINB_r_U(n=1,prob=Phi[1],theta=exp(log(T[n])+beta_0[g]+rho[1,g]),disp_param=size[1])
        cell_matrix[n,g]<- out_1$Y
        U_ng[n,g]<- out_1$U}
      
      else if(clust[n]==2){
        out_2 = ZINB_r_U(n=1,prob=Phi[2],theta=exp(log(T[n])+beta_0[g]+rho[2,g]),disp_param=size[2])
        cell_matrix[n,g]<- out_2$Y
        U_ng[n,g]<- out_2$U}
      
    }
  }
  
  # cell_matrix[1:5,]
  # U_ng[1:5,]
  # Y<- cell_matrix
  #  Y <- Y[order(clust),]
  #  Y
  #  Y <- as.data.frame(Y)
  # pheatmap(Y, cluster_cols = F, cluster_rows = F)  
  #  
  ########################################################################################
  ########################################################################################
  
  
  
  
  
  
  
  
  
  
  
  ########################################################################################
  ########################################################################################
  Z_nk<- matrix(data=NA,ncol=K,nrow=N)
  
  for(i in 1:N){
    if(clust[i]==1) Z_nk[i,1]<- 1
    else Z_nk[i,1]<- 0
    if(clust[i]==2) Z_nk[i,2]<- 1
    else Z_nk[i,2]<- 0
  }
  # Z_nk
  ###################################################################################
  
  result<- NR_ZINB(Y=cell_matrix,totals=T,beta_0=beta_0,rho=rho,Z_nk=Z_nk,U_ng=U_ng,overdisp_para=size,criterion=1e-5,maxiter=500)

    #NR_Poisson_ZIP(Y=Y,T=T,beta_0=beta_0,rho=rho,Z_nk=Z_nk,U_ng=U_ng,K=K,criterion=1e-4,maxiter=500)
  #print(result)

  


  beta_0_est[s,]<- result$beta_0
  rho1_est[s,]<- result$rho[1,]
  rho2_est[s,]<- result$rho[2,]
  
  iter<- c(iter,result$iteration)
  
}
  
  
  
  
###################################################################################################################################################
##########################################################################################  
 

print(beta_0_est)
print(rho1_est)
print(rho2_est)
print(iter)


boxplot(beta_0_est,main="beta_0g")
#abline(h=0.25)
abline(h=0.85)
#abline(h=0.35)


boxplot(rho1_est,main="rho1")
abline(h=2)
#abline(h=0)
abline(h=-2)
#abline(h=1.5)
#abline(h=0.5)
#abline(h=-0.5)

boxplot(rho2_est,main="rho2")

abline(h=2)
#abline(h=0)
abline(h=-2)
#abline(h=1.5)
#abline(h=0.5)
#abline(h=-0.5)



############################################################################################  




