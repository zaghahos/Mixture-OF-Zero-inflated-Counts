rm(list=ls())
ls()
#library(pheatmap)

####################################################################################
ZIP_r<- function(n,pi,theta){
  Z<- rbinom(size=1,n=n,p=pi)
  Y<- NULL
  for(i in 1:n){
    if(Z[i]==1) Y<- c(Y,0)
    else Y<- c(Y,rpois(1,lambda = theta))
  }
  return(Y)
}
#################################################################################
K=2
N=1000
G=4
Phi<- c(0.1,0.2)
P<- c(0.5,0.5)
P1<- c(0.7,0.3)

#S<- 3
#S<- 50
#beta_0_est<- matrix(NA,nrow=S,ncol=G)
#rho1_est<- matrix(NA,nrow=S,ncol=G)
#rho2_est<- matrix(NA,nrow=S,ncol=G)

#iter<- NULL
#beta_0_est
#rho_est

#for(s in 1:S){
  
  
Z <- sample(c(1:K),size=N,replace = TRUE,prob=P)

  
  
  
  ### creating lambdas 
  
  T <- round(rnorm(N,mean=10,sd=0.5)) ## vector with known Totals#hist(T)
  #T
  beta_0 <- rep(0.25,G) ## vector of beta_0 for all genes (baseline)
  #rho_1 <- c(rgamma(1:3,shape=20,rate=10),rep(0,2)) ## fixed effect due to cluster 1, no effect for the last 5 genes
  #rho_2 <- c(rep(0,3),rgamma(1:2,shape=20,rate=10)) ## fixed effect due to cluster 2, no effect for the first 5 genes
  
  
  rho_1 <- c(rep(2,round(G/2)),rep(0,round(G/2))) ## fixed effect due to cluster 1, no effect for the last 5 genes
  rho_2 <- c(rep(0,round(G/2)),rep(2,round(G/2))) ## fixed effect due to cluster 2, no effect for the first 5 genes
  rho <- rbind(rho_1,rho_2)
  #rho
  x<- rbinom(N,1,0.5)
  beta_pg<- c(rep(1,round(G/2)),rep(0.5,round(G/2)))
  
  
  #x<- rnorm(N)
  #beta_gp<- rnorm(G,2)
  #beta_gp
  
  #beta_gp%*%t(x)
  #beta_gp[2]*x
  
  #x
  #beta<- cbind(beta_0,beta_gp)
  #beta
  #X<- cbind(1,x)
  #X
  #beta[1,]*X
  #X
  #beta[1,]*X[,1]
  #beta%*%t(X)
  
  Y <- matrix(NA,nrow=N,ncol=G)
  for(n in 1:N){
    if(Z[n] == 1){
      for(g in 1:G){
        Y[n,g] <- ZIP_r(1,pi=Phi[1],theta =exp(log(T[n]) + beta_0[g] + rho[1,g]+beta_pg[g]*x[n]))
      }
    }
    if(Z[n] == 2){
      for(g in 1:G){
        Y[n,g] <- ZIP_r(1,pi=Phi[2],theta=exp(log(T[n]) + beta_0[g] + rho[2,g]+beta_pg[g]*x[n]))
      }
    }
  }
  Y
  
  #Y <- Y[order(Z),]
  #Y
  #Y <- as.data.frame(Y)
  #pheatmap(Y, cluster_cols = F, cluster_rows = F)
  
  
  
  ########################################################################################
  ########################################################################################
  Z_nk<- matrix(data=NA,ncol=K,nrow=N)
  #Z_nk
  for(i in 1:N){
    if(Z[i]==1) Z_nk[i,1]<- 1
    else Z_nk[i,1]<- 0
    if(Z[i]==2) Z_nk[i,2]<- 1
    else Z_nk[i,2]<- 0
  }
  Z_nk
  ###################################################################################
  U_ng<- matrix(data=NA,ncol=G,nrow=N)
  for(n in 1:N){
    for(g in 1:G){
      if(Y[n,g]==0) U_ng[n,g]<- 1
      else U_ng[n,g]<- 0
    }
  }
  U_ng
  #1-U_ng
  
  #########################################################################################
  NR_Poisson_ZIP_1X<- function(Y,T,beta_0,rho,beta_pg,x,Z_nk,U_ng,criterion=1e-4,maxiter=500){
    N<- nrow(Y)
    G<- ncol(Y)
    K<- dim(Z_nk)[2]
    V<- 1
    
    
    ### Initialization###
    mu<- array(data=NA,c(N,G,K))
    diff_ym<- diff_ymx<- Z_diffym<- Z_diffymx <- array(data=NA,c(N,G,K))
    U_diff_Y_mu<- U_diff_Y_mu_X<- array(data=NA,c(N,G,K))
    M_U<- MX_U<- MX2_U<- array(NA,c(N,G,K))
    
    
    ZM<- ZMX<- ZMX2<- array(NA,c(N,G,K))
    #grad<- numeric(((G*K)+G+(G*V)))
    #Hessian<- matrix(0,nrow=(G*(K+1+V)),ncol=(G*(K+1+V)))
    theta<- numeric(((G*K)+G+(G*V)))
  
    
    ### Assign values to theta vector###
    theta[1:G]<- beta_0
    for(k in 1:K){
      theta[((k*G)+1):((k*G)+G)]<- rho[k,]
    }
    theta[(G+(K*G)+1):(G+(G*K)+(G*V))]<- beta_pg
    #print(theta)
    #BX<- x%*%t(beta_pg)
    
    it<- 0
    diff_theta<- 1000
    while(it<maxiter & all(diff_theta>criterion)){
    
    #while(it<maxiter & diff_theta>criterion){
      it<- it+1
      #print(it)
      ### mu###
      mu<- array(data=NA,c(N,G,K))
      #print(beta_0)
      #print(rho)
      #print(beta_pg)
      for(n in 1:N){
        for(g in 1:G){
          for(k in 1:K){
            mu[n,g,k]<- exp(log(T[n])+beta_0[g]+rho[k,g]+beta_pg[g]*x[n])
          }
        }
      }
      #print(mu[1:10,,])
      
      ##(1-U)*(Y-mu), Z*(1-U)*(Y-mu),Z*(1-U)*(Y-mu)*X##
      for(k in 1:K){
        diff_ym[,,k]<- (Y-mu[,,k])
        diff_ymx[,,k]<- x*(Y-mu[,,k])
        
      }
      
      for(k in 1:K){
        U_diff_Y_mu[,,k]<- (1-U_ng[,,k])*diff_ym[,,k]
        U_diff_Y_mu_X[,,k]<- (1-U_ng[,,k])*diff_ymx[,,k]
        #U_mu[,,k]<- (1-U_ng)*mu_ng[,,k]
      }
      
      for(k in 1:K){
        Z_diffym[,,k]<- Z_nk[,k]*U_diff_Y_mu[,,k]
        Z_diffymx[,,k]<- Z_nk[,k]*U_diff_Y_mu_X[,,k]
        
      }
      
    
      
      ##Z*(1-U)*mu,Z*(1-U)*mu*x,Z*(1-U)*mu*X^2 ##
      for(k in 1:K){
        M_U[,,k]<- (1-U_ng[,,k])*mu[,,k]
        MX_U[,,k]<- (1-U_ng[,,k])*mu[,,k]*x
        MX2_U[,,k]<- (1-U_ng[,,k])*mu[,,k]*(x^2)
        
        
      }
     

      for(k in 1:K){
        ZM[,,k]<- Z_nk[,k]*M_U[,,k]
        ZMX[,,k]<- Z_nk[,k]*MX_U[,,k]
        ZMX2[,,k]<- Z_nk[,k]*MX2_U[,,k]
        
        
      }
    
      
      #print(ZM[1:10,,])
      #print(ZMX[1:10,,])
      #print(ZMX2[1:10,,])
      #return(ZM[1:10,,])
      #####################################################################################
      ### derivatives of beta_0###
      db0<- apply(Z_diffym,2,sum)
      ddb0<- apply(ZM,2,sum)
      #print(db0)
      #print(ddb0)
      ### derivatives of beta_p###
      db1<- apply(Z_diffymx,2,sum)
      ddb1<- apply(ZMX2,2,sum)
      #print(db1)
      #print(ddb1)
      ### derivatives of rho###
      drho_k<- apply(Z_diffym,c(3,2),sum)
      ddrho_k<- apply(ZM,c(3,2),sum)
      #print(drho_k)
      #print(ddrho_k)
      ### derivatives of beta_0,beta_p###
      ddb10<- apply(ZMX,2,sum)
      #print(ddb10)
      ### derivatives of beta_1,rho###
      ddb1rho<- apply(ZMX,c(3,2),sum)
      #print(ddb1rho)
      ### derivatives of beta_0 , rho###
      ddb0rho<- apply(ZM,c(3,2),sum)
      
      
      ### Gradient Vecotr###
      grad<- numeric(((G*K)+G+(G*V)))
      
      grad[1:G]<- db0
      
      for(k in 1:K){
        #  for(g in 1:G){
        grad[((G*k)+1):((G*k)+G)]<- drho_k[k,]
        
        
        # }
        
      }
      grad[((G*K)+G+1):((G*K)+G+(G*V))]<- db1
      #return(grad)  
      ### Hessian Matrix###
      Hessian<- matrix(0,nrow=(G*(K+1+V)),ncol=(G*(K+1+V)))
      #print(Hessian)
      #print(ddb0)
      #print(ddb1)
      #print(ddrho_k)
      #print(ddb10)
      
      Hessian[(1:G),(1:G)]<- Hessian[(1:G),(1:G)]-diag(ddb0)
      Hessian[(G+(G*K)+1):(G+(G*K)+(G*V)),(G+(G*K)+1):(G+(G*K)+(G*V))]<- 
        Hessian[(G+(G*K)+1):(G+(G*K)+(G*V)),(G+(G*K)+1):(G+(G*K)+(G*V))]-diag(ddb1)
      #return(Hessian)
      
      
      for(k in 1:K){
        Hessian[((k*G)+1):((k*G)+G),((k*G)+1):((k*G)+G)]<- Hessian[((k*G)+1):((k*G)+G),((k*G)+1):((k*G)+G)]-diag(ddrho_k[k,])
        Hessian[1:G,((k*G)+1):((k*G)+G)]<- Hessian[1:G,((k*G)+1):((k*G)+G)]-diag(ddrho_k[k,])
        Hessian[((k*G)+1):((k*G)+G),1:G]<- Hessian[((k*G)+1):((k*G)+G),1:G]-diag(ddrho_k[k,])
        Hessian[((k*G)+1):((k*G)+G),((K*G)+G+1):((K*G)+G+(G*V))]<- Hessian[((k*G)+1):((k*G)+G),((K*G)+G+1):((K*G)+G+(G*V))]-diag(ddb1rho[k,])
        Hessian[((K*G)+G+1):((K*G)+G+(G*V)),((k*G)+1):((k*G)+G)]<- Hessian[((K*G)+G+1):((K*G)+G+(G*V)),((k*G)+1):((k*G)+G)]-diag(ddb1rho[k,])
        
        
        
        
      }
      
      #return(Hessian)
      Hessian[1:G,((K*G)+G+1):((K*G)+G+(G*V))]<- Hessian[1:G,((K*G)+G+1):((K*G)+G+(G*V))]-diag(ddb10)
      Hessian[((K*G)+G+1):((K*G)+G+(G*V)),(1:G)]<- Hessian[((K*G)+G+1):((K*G)+G+(G*V)),(1:G)]-diag(ddb10)
      
      #Hessian[1:G,((K*G)+G+1):((K*G)+G+(G*V))]<- Hessian[1:G,((K*G)+G+1):((K*G)+G+(G*V))]-diag(ddb10)
      #Hessian[((K*G)+G+1):((K*G)+G+(G*V)),(1:G)]<- Hessian[((K*G)+G+1):((K*G)+G+(G*V)),(1:G)]-diag(ddb10)
      
      
      #print(Hessian)
      #return(Hessian)
      
      ##############################################################################################
      
      
      ### Updating theta using Newton-Raphson Algorithm###
      theta_new<- theta-qr.coef(qr(Hessian,tol=1e-300), grad) ### tol=1e-300
      
      for(i in 1:length(theta_new)){
        if(is.na(theta_new[i])==TRUE) theta_new[i]<- theta[i]
        #if(is.nan(theta_new[i])==TRUE) theta_new[i]<- theta[i]
        
      }
      #return(theta_new)
      #diff_theta<- mean(abs(theta_new-theta))
      diff_theta<- abs(theta_new-theta)
      
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
      
      beta_pg<- theta[(G+(K*G)+1):(G+(G*K)+(G*V))]
      
      
      
      #out<- list(beta_0=beta_0,rho=rho,beta_pg=beta_pg,theta=theta,iteration=it)
      
      
    }
    
      out<- list(beta_0=beta_0,rho=rho,beta_P=beta_pg,theta=theta,iteration=it)
    
    return(out)
    
    
    #  }
    
    
  }
  
  
  ##############################################################################################
  result<- NR_Poisson_ZIP_1X(Y,T,beta_0,rho,beta_pg,x,Z_nk,U_ng,K,criterion=1e-4,maxiter=500)
#result  
print(result)
  #st<- system.time(NR_Poisson_ZIP(Y=Y,T=T,beta_0=beta_0,rho=rho,Z_nk=Z_nk,U_ng=U_ng,K=K,criterion=1e-5,maxiter=500))  
  
########################################################################################################################################  
###########################################################################################################################  
###########################################################################################################################  
### EM Function###
#(Y,T,beta_0,rho,beta_pg,x,Z_nk,K,criterion=1e-5,maxiter=500)

ZIP_Mix_EM_1X <- function(Y,totals,beta_0g,rho_gk,Beta_pg,X,P,Phi,eps=1e-8,maxiter_EM=1000){  
  N<- nrow(Y)
  G<- ncol(Y)
  K<- length(P)
  dl<- 1+eps
  iter<- 0
  ll<- rep(0,maxiter_EM+1)
  
  beta_0g_cur<- beta_0g
  rho_cur<- rho_gk
  beta_pg_cur<- Beta_pg
  
  Z_nk<- matrix(data=NA,nrow=N,ncol=K)
  U_ngk<- array(data=NA,c(N,G,K))
  
  ZU<- array(data=NA,c(N,G,K))
  
  thresh<- -744
  
  
  lambda<- array(data=NA,c(N,G,K))
  #print(beta_0)
  #print(rho)
  #print(beta_pg)
  for(n in 1:N){
    for(g in 1:G){
      for(k in 1:K){
        lambda[n,g,k]<- exp(log(totals[n])+beta_0g[g]+rho_gk[k,g]+Beta_pg[g]*X[n])
      }
    }
  }

  
  
  b<- array(NA,c(N,G,K))

  for(k in 1:K){
    for(n in 1:N){
      for(g in 1:G){
        if(Y[n,g]==0) {b[n,g,k]<- P[k]*(Phi[k]+(1-Phi[k])*exp(-lambda[n,g,k]))}
        else {b[n,g,k]<- P[k]*((1-Phi[k])*dpois(Y[n,g],lambda = lambda[n,g,k]))}
        #b[n,g,k]<- P[k]*dpois(Y[n,g],lambda=lambda[n,g,k])
      }
    }
  }
  
  #print(b)
  
  b1<- array(NA,c(N,G,K))
  
  for(k in 1:K){
    for(n in 1:N){
      for(g in 1:G){
        if(Y[n,g]==0) {b1[n,g,k]<- (Phi[k]+(1-Phi[k])*exp(-lambda[n,g,k]))}
        else {b1[n,g,k]<- ((1-Phi[k])*dpois(Y[n,g],lambda = lambda[n,g,k]))}
        #b1[n,g,k]<- dpois(Y[n,g],lambda=lambda[n,g,k])
      }
    }
  }
  
  ##############################################################################
  d<- array(data=NA,c(N,G,K))
  
  for(n in 1:N){
    for(g in 1:G){
      for(k in 1:K){
        if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*exp(-lambda[n,g,k])))
        else d[n,g,k]<- 0
        
      }
    }
  }
  ##############################################################################
  
  #print(b1)
  
  b_sum<- apply(b,c(1,2),sum)
  
  l<- sum(rowSums(log(b_sum))) ## observed data log-likelihood
  
  
  ratio<- 1000
  oldl<- -1000000
  #print(b1)
  #return(l)
  #  Theta_cur<- Theta
  
  
  while(abs(dl)>eps & iter<maxiter_EM){
    #while(iter<maxiter_EM){
    
    iter<- iter+1
    ll[iter]<- l
    #pi_Update <- colSums(Zhat)/length(Y)
    
    ##E-step##
    ###Estimate Z_nk###
    
    Z_nk<- matrix(data=log(P),ncol=K,nrow=N,byrow=T)
    
    for(k in 1:K){
      for(g in 1:G){
        Z_nk[,k]<- Z_nk[,k]+(log(b1[,g,k]))
      }
    }
    #print(Z_nk)
    #Z_nk<- exp(Z_nk)
    #Z_nk<- Z_nk/apply(Z_nk,1,sum)
    #print(Z_nk)
    
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
    
    #P_hat<- apply(Z_nk,2,mean)
    #print(P_hat)
    
    
    epsilon<- 1e-10
    sl<- length(Z_nk[Z_nk<epsilon])
    bl<- length(Z_nk[Z_nk>1-epsilon])
    Z_nk[Z_nk<epsilon]<- rep(epsilon,sl)
    Z_nk[Z_nk>1-epsilon]<- rep(1-epsilon,bl)
    Z_nk<- Z_nk/apply(Z_nk,1,sum)
    #print(Z_nk)
    ###############################################
    ###Estimate U_ng###
    
    # denom<- apply(d,c(1,2),sum)
    # 
    # for(n in 1:N){
    #   for(g in 1:G){
    #     if(Y[n,g]==0) U_ng[n,g]<- sum(P*Phi)/denom[n,g]
    #     else U_ng[n,g]<- 0
    #     
    #   }
    # }
    
    
    
    U_ngk<- array(data=NA,c(N,G,K))
    
    #denom<- apply(d,c(1,2),sum)
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          #print(P[k]*Phi[k])
          
          if(Y[n,g]==0) U_ngk[n,g,k]<- (P[k]*Phi[k])/d[n,g,k]
          else U_ngk[n,g,k]<- 0
          
        }
      }
    }
    
  
    
    
    #####################################################
    ### M-step###
  
    ##Update P###
    P<- apply(Z_nk,2,mean)
    
    
    #print(P)
    #return(P)  
    
    ### Update Phi###  
    # for(k in 1:K){
    #   ZU[,,k]<- Z_nk[,k]*U_ng
    # }
    
    for(k in 1:K){
      ZU[,,k]<- Z_nk[,k]*U_ngk[,,k]
    }
    
    Phi<- (apply(ZU,3,sum))/(G*apply(Z_nk,2,sum))
    #print(Phi_hat)  
    #return(P)  
    
    ### Updaet Parameters of lambda: rho_gk, Beta_0g### 
    
    ### NR POisson##
    #Newton_Update<- NR_Poisson(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,K=K,criterion=1e-4,maxiter=500)
    
    #Newton_Update<- NR_Poisson_G(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,K=K,criterion=1e-5,maxiter=500)
    
    ### NR ZIP##
    #Newton_Update<- NR_Poisson_G_ZIP(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,U_ng=U_ng,K,criterion=1e-5,maxiter=500)
    #Newton_Update<- NR_Poisson_ZIP(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,Z_nk=Z_nk,U_ng=U_ng,K,criterion=1e-5,maxiter=500)
    Newton_Update<- NR_Poisson_ZIP_1X(Y=Y,T=T,beta_0=beta_0g_cur,rho=rho_cur,beta_pg=beta_pg_cur,x=X,Z_nk=Z_nk,U_ng=U_ngk,criterion=1e-4,maxiter=500)
    
    
    beta_0g<- Newton_Update$beta_0
    rho_gk<- Newton_Update$rho
    Beta_pg<- Newton_Update$beta_P
    ##############################################################################################################
    #lambda<- array(data=NA,c(N,G,K))
    
    for(n in 1:N){
      for(g in 1:G){
        for(k in 1:K){
          lambda[n,g,k]<- exp(log(totals[n])+beta_0g[g]+rho_gk[k,g]+Beta_pg[g]*X[n])
        }
      }
    }
    
    #b<- array(NA,c(N,G,K))
    #print(b)
    
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0) {b[n,g,k]<- P[k]*(Phi[k]+(1-Phi[k])*exp(-lambda[n,g,k]))}
          else {b[n,g,k]<- P[k]*((1-Phi[k])*dpois(Y[n,g],lambda = lambda[n,g,k]))}
          #b[n,g,k]<- P[k]*dpois(Y[n,g],lambda=lambda[n,g,k])
        }
      }
    }
    
    #print(b)
    
    #b1<- array(NA,c(N,G,K))
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0) {b1[n,g,k]<- (Phi[k]+(1-Phi[k])*exp(-lambda[n,g,k]))}
          else {b1[n,g,k]<- ((1-Phi[k])*dpois(Y[n,g],lambda = lambda[n,g,k]))}
          #b1[n,g,k]<- dpois(Y[n,g],lambda=lambda[n,g,k])
        }
      }
    }
    
    
    for(n in 1:N){
      for(g in 1:G){
        for(k in 1:K){
          if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*exp(-lambda[n,g,k])))
          else d[n,g,k]<- 0
          
        }
      }
    }
    
    
    # beta_0g_cur<- Newton_Update$beta_0
    # rho_cur<- Newton_Update$rho
    # beta_pg_cur<- Newton_Update$beta_P
    
    
    beta_0g_cur<- beta_0g
    rho_cur<- rho_gk
    beta_pg_cur<- Beta_pg
    
    
    
    oldl<- l
    
    b_sum<- apply(b,c(1,2),sum)
    
    l<- sum(rowSums(log(b_sum))) 
    
    
    
    #Pi_Update <- colSums(Z_hat)/N
    
    ## updating Zhat or p_ik's
    
    
    ################################################################
    
    #print(l)
    dl<- l-oldl
  }
  cat("number of iterations=", iter, "\n")
  iter <- iter+1
  ll[iter] <- l
  postprobs <- Z_nk
  colnames(postprobs) <- c(paste("comp", ".", 1:K, sep = ""))
  
  
  out <- list(P=P,phi=Phi,beta_0g=beta_0g_cur,rho_gk=rho_cur,beta_pg_update=beta_pg_cur,no_iter=iter,
              loglik=l,all.loglik=ll[1:iter],restart=0,ft="ZIPmix")
  class(out)<- "mixEM"
  out
  
}

##############################################################################################################################
my_out<- ZIP_Mix_EM_1X(Y=Y,totals=T,beta_0g=beta_0,rho_gk=rho,Beta_pg=beta_pg,X=x,P=P,Phi=Phi,K=K,eps=1e-8,maxiter_EM=1000)
print(my_out)






########################################################################################################################################

  #beta_0_init <- rep(0.3,G) 
  #rho_1_init <- c(rep(2.5,round(G/2)),rep(0.5,round(G/2))) ## fixed effect due to cluster 1, no effect for the last 5 genes
  #rho_2_init <- c(rep(0.5,round(G/2)),rep(2.50,round(G/2))) ## fixed effect due to cluster 2, no effect for the first 5 genes
  #rho_init <- rbind(rho_1_init,rho_2_init) 
  
  #result<- NR_Poisson_ZIP(Y=Y,T=T,beta_0=beta_0_init,rho=rho_init,Z_nk=Z_nk,U_ng=U_ng,K=K,criterion=1e-4,maxiter=500)
  
  #st<- system.time(NR_Poisson_ZIP(Y=Y,T=T,beta_0=beta_0_init,rho=rho_init,Z_nk=Z_nk,U_ng=U_ng,K=K,criterion=1e-5,maxiter=500))
  print(st)
  
  beta_0_est[s,]<- result$beta_0
  rho1_est[s,]<- result$rho[1,]
  rho2_est[s,]<- result$rho[2,]
  
  iter<- c(iter,result$iteration)
  
#}

print(beta_0_est)
print(rho1_est)
print(rho2_est)
print(iter)
print(st)


boxplot(beta_0_est,main="beta_0g")
abline(h=0.25)
#abline(h=0.15)
#abline(h=0.35)


boxplot(rho1_est,main="rho1")
abline(h=2)
abline(h=0)
#abline(h=2.5)
#abline(h=1.5)
#abline(h=0.5)
#abline(h=-0.5)

boxplot(rho2_est,main="rho2")
abline(h=2)
abline(h=0)
#abline(h=2.5)
#abline(h=1.5)
#abline(h=0.5)
#abline(h=-0.5)


