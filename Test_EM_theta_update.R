rm(list = ls())

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
#################################################################################################################
#################################################################################################################

##################################################################################################################
###################################################################################################################
disp_para<- function (y, mu, n = sum(w), w, limit = 500, eps = .Machine$double.eps^0.25, trace = FALSE) {
  
  score <- function(n, th, mu, y, w) sum(w * (digamma(th +y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +th)/(mu + th)))
  info <- function(n, th, mu, y, w) sum(w * (-trigamma(th +y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +th)^2))
  if (inherits(y, "lm")) {
    mu <- y$fitted.values
    y <- if (is.null(y$y)) 
      mu + residuals(y)
    else y$y
  }
  #if (missing(weights)) 
  if (missing(w)) 
    
    #  weights <- rep(1, length(y))
    w <- rep(1, length(y))
  
  #t0 <- n/sum(weights * (y/mu - 1)^2)
  t0 <- n/sum(w * (y/mu - 1)^2)
  
  it <- 0
  del <- 1
  if (trace) 
    message(sprintf("theta.ml: iter %d 'theta = %f'", it, 
                    signif(t0)), domain = NA)
  while ((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    #del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, mu, y, weights))
    
    del <- score(n, t0, mu, y, w)/(i <- info(n, t0, mu, y, w))
    t0 <- t0 + del
    if (trace) 
      message("theta.ml: iter", it, " theta =", signif(t0))
  }
  if (t0 < 0) {
    t0 <- 0
    warning("estimate truncated at zero")
    attr(t0, "warn") <- gettext("estimate truncated at zero")
  }
  if (it == limit) {
    warning("iteration limit reached")
    attr(t0, "warn") <- gettext("iteration limit reached")
  }
  attr(t0, "SE") <- sqrt(1/i)
  t0
}
#######################################################################################################
#######################################################################################################
To_vec_transfer<- function(Y,Z,U,MU){
  N<- nrow(Y)
  G<- ncol(Y)
  K<- ncol(Z)
  
  vec_y<- as.vector(t(Y))
  ZU<- array(NA,c(N,G,K))
  for(k in 1:K){
    ZU[,,k]<- Z[,k]*(1-U)
  }
  
  
  mu_mat<- matrix(NA,nrow=N*G,ncol=K)
  for(k in 1:K){
    mu_mat[,k]<- rep(MU[k,],N)
  }
  
  weight_mat<- matrix(NA,nrow=N*G,ncol=K)
  for(k in 1:K){
    weight_mat[,k]<- as.vector(t(ZU[,,k]))
  }
  
  out<- list(Y_vec=vec_y,WW=weight_mat,MU_vec=mu_mat)
  return(out)
  
} 


########################################################
##### THE EM Function for Mixture of MIZNB####
###################################################################################################
ZINB_mix_K_comp<- function(Y,P,Phi,clust,U_ng,mu,dispresiopn_parameter,eps=1e-5, maxit=1000){
  dl <- 1 + eps
  iter<-0
  ll <- rep(0, maxit+1)
  G<- ncol(Y)
  N<- nrow(Y)
  K<- length(Phi)
  
  nu<- dispresiopn_parameter
  
  bk<- array(NA,c(N,G,K))
  for(k in 1:K){
    for(n in 1:N){
      for(g in 1:G){
        if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k]))) 
        
        
        else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])
      }
    }
  }
  
  
  b<- array(NA,c(N,G,K))
  
  for(k in 1:K) 
  {b[,,k]<- P[k]*bk[,,k]}
  
  bp_sum<- apply(b,c(1,2),sum)
  l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
  Z_nk<- matrix(NA,ncol=K,nrow=N)
  thresh<- -744
  
  
  #U_ng<- matrix(data=NA,nrow=N,ncol=G)
  ZU<- array(data=NA,c(N,G,K))
  
  
  
  d<- array(data=NA,c(N,G,K))
  
  for(n in 1:N){
    for(g in 1:G){
      for(k in 1:K){
        if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))
        
        else d[n,g,k]<- 0
        
      }
    }
  }
  
  
  
  
  
  while (abs(dl)>eps & iter<maxit) { ## Here I use abs(dl) rather than dl#
    iter<-iter+1
    ll[iter] <- l
    
    
    # Z_nk<- matrix(data=log(P),ncol=K,nrow=N,byrow=T)
    # 
    # for(k in 1:K){
    #   for(g in 1:G){
    #     Z_nk[,k]<- Z_nk[,k]+(log(bk[,g,k]))
    #   }
    # }
    # 
    # if(K>1){
    #   v1<- which(apply(Z_nk,1,max)< thresh)
    #   v3<- 1:N
    #   len<- length(v1)
    #   #print(len)
    #   if(len>0){
    #     v2<- apply(array(Z_nk[v1,],dim=c(len,K)),1,order)[K,]
    #     ddd<- cbind(v1,v2)
    #     #print(ddd)
    #     Z_nk[v1,]<- 0
    #     Z_nk[ddd]<- 1
    #     v3<- -v1
    #   }
    #   
    #   Z_nk[v3,]<- exp(Z_nk[v3,])
    # }
    # Z_nk<- Z_nk/apply(Z_nk,1,sum)
    # #print(Z_nk)
    # 
    # epsilon<- 1e-10
    # sl<- length(Z_nk[Z_nk<epsilon])
    # bl<- length(Z_nk[Z_nk>1-epsilon])
    # Z_nk[Z_nk<epsilon]<- rep(epsilon,sl)
    # Z_nk[Z_nk>1-epsilon]<- rep(1-epsilon,bl)
    # Z_nk<- Z_nk/apply(Z_nk,1,sum)
    # #print(Z_nk)
    # 
    # 
    # 
    
    #Z_nk#
    Z_nk<- matrix(0,nrow=N,ncol=length(P))
    for(i in 1:N){
      if(clust[i]==1)Z_nk[i,1]<- 1
      else Z_nk[i,2]<- 1
    }
    #Z
    ##################################################
    ##################################################
    ## U_ng##
    # U_ng<- matrix(NA,nrow=N,ncol=G)
    # for(n in 1:N){
    #   for(g in 1:G){
    #     if(Y[n,g]==0) U_ng[n,g]<- 1
    #     else U_ng[n,g]<-0
    #   }
    # }
    
    
    P<- apply(Z_nk,2,mean)
    #print(P)
    
    U_ng<- U_ng
    #print(U_ng)
    
    # U_ng<- matrix(data=NA,nrow=N,ncol=G)
    # 
    # denom<- apply(d,c(1,2),sum)
    # 
    # for(n in 1:N){
    #   for(g in 1:G){
    #     if(Y[n,g]==0) U_ng[n,g]<- sum(P*Phi)/denom[n,g]
    #     else U_ng[n,g]<- 0
    #     
    #   }
    # }
    # 
    
    for(k in 1:K){
      ZU[,,k]<- Z_nk[,k]*U_ng
    }
    
    Phi<- (apply(ZU,3,sum))/(G*apply(Z_nk,2,sum))
    ###############################################################################################
    ################################################################################################
    
    # 
    # mult_Z_u<- array(NA,c(N,G,K))
    # mult_Z_oneminus_U<- array(NA,c(N,G,K))
    # mult_Z_oneminus_U_Y<- array(NA,c(N,G,K))
    # 
    # 
    # 
    # for(k in 1:K){
    #   mult_Z_u[,,k]<- Z_nk[,k]*U_ng
    #   mult_Z_oneminus_U[,,k]<- Z_nk[,k]*(1-U_ng)
    #   mult_Z_oneminus_U_Y[,,k]<- Z_nk[,k]*(1-U_ng)*Y
    #   
    # }
    # 
    # 
    # mu<- t(apply(mult_Z_oneminus_U_Y,c(2,3),sum)/apply(mult_Z_oneminus_U,c(2,3),sum))
    
    mu<- mu
    ##########################################################################################
    ##########################################################################################
    ### To estimate size_MLE###
    result<- To_vec_transfer(Y=cell_matrix,Z=Z_nk,U=U_ng,MU=mu)
    
    size_MLE<- rep(NA,K)
    for(k in 1:K){
      size_MLE[k]<- disp_para(y=result$Y_vec, mu=result$MU_vec[,k], n = sum(result$WW[,k]), w=result$WW[,k], limit = 500, eps = .Machine$double.eps^0.25, trace = FALSE) 
      
    }
    
    nu<- size_MLE
    print(nu)
    ################################################################################################################
    ################################################################################################################
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k]))) 
          else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])
        }
      }
    }
    
    ############################################################
    for(n in 1:N){
      for(g in 1:G){
        for(k in 1:K){
          if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))
          else d[n,g,k]<- 0
          
        }
      }
    }
    
    #########################################################
    
    
    
    for(k in 1:K){
      b[,,k]<- P[k]*bk[,,k]
    }
    bp_sum<- apply(b,c(1,2),sum)
    #print(bp_sum)
    oldl<- l
    
    l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
    #print(l)
    
    ################################################################
    
    #print(l)
    dl<- l-oldl
  }
  cat("number of iterations=", iter, "\n")
  iter <- iter+1
  ll[iter] <- l
  postprobs <- Z_nk
  colnames(postprobs) <- c(paste("comp", ".", 1:K, sep = ""))
  
  
  out <- list(P=P,Phi=Phi,MU=mu,size_est=nu,a_hat=1/nu,postprobs=postprobs,
              loglik=l,all.loglik=ll[1:iter],restart=0,ft="ZINBmix")
  class(out)<- "mixEM"
  out
  
}
###################################################################################################
###################################################################################################
### Data sim_1###
### Generate Y Matrix_Matrix of raw counts###
### case_1###

N<- 60
G<- 120
P<- c(0.5,0.5)
Phi<- c(0.1,0.1)
K<- length(P)
#size<- c(10,5)
#size<- c(5,10)
size<- c(5,10)

#P<- c(0.25,0.25,0.25,0.25)
#N<- 1000
#G<- 10
#Phi<- c(0.1,0.1,0.2,0.2)
#P<- c(0.25,0.25,0.25,0.25)
#P1<- c(0.4,0.3,0.2,0.1)


#mu_1<- rgamma(n=G,shape=10,scale=1/2)
#mu_2<- rgamma(n=G,shape=10,scale=1/1.5)


#mu_1<- c(rep(5,(G/2)),rep(10,G/2))
#mu_2<- c(rep(10,(G/2)),rep(5,G/2))

mu_1<- c(rep(5,G))
mu_2<- c(rep(10,G))

MU<- rbind(mu_1,mu_2)
#MU

#S<- 256
#S<- 128
#S<- 50
S<- 1
P_hat<- matrix(NA,ncol=length(P),nrow=S)
#P_hat<- NULL
Phi_hat<- matrix(NA,ncol=length(P),nrow=S)
overdisp_hat<- matrix(NA,nrow=S,ncol=K)
size_hat<- matrix(NA,nrow=S,ncol=K)

#Phi_hat
#phi_hat2<- NULL
mu1_hat<- matrix(NA,ncol = G,nrow=S)
mu2_hat<- matrix(NA,ncol = G,nrow=S)
#lambda3_hat<- matrix(NA,ncol = G,nrow=S)
#lambda4_hat<- matrix(NA,ncol = G,nrow=S)

#MU_1<- MU_2<- matrix(NA,nrow=S,ncol=G)
#V_M<- rep(NA,S)
elapsed<- rep(NA,S)

for(s in 1:S){
  
  # mu_1<- rgamma(n=G,shape=10,scale=1/2)
  
  
  #mu_2<- rgamma(n=G,shape=10,scale=1/0.2)
  
  #MU<- rbind(mu_1,mu_2)
  
  
  
  #for(s in 1:S){
  
  
  #for(iteration in 1:r){
  
  #p_hat<- NULL
  #phi_hat1<- NULL
  #phi_hat2<- NULL
  #lambda1_hat<- matrix(0,ncol = G,nrow=r)
  #lambda2_hat<- matrix(0,ncol = G,nrow=r)
  
  #for(iter in 1:r){
  #clust<- rbinom(size=1,n=N,p=0.5)
  clust<- sample(1:2,size=N,replace = TRUE,prob=P)
  #clust
  #length(clust[clust==1])
  
  
  #clust<- sample(c(1,2),size=N,replace = TRUE,prob=pn) 
  #clust
  
  
  cell_matrix<- U_ng<-  matrix(rep(0,N*G),nrow=N,ncol=G)
  
  for(n in 1:N){
    for(g in 1:G){
      if(clust[n]==1) {
        out_1 = ZINB_r_U(n=1,prob=Phi[1],theta=MU[1,g],disp_param=size[1])
        cell_matrix[n,g]<- out_1$Y
        U_ng[n,g]<- out_1$U}
      
      else if(clust[n]==2){
        out_2 = ZINB_r_U(n=1,prob=Phi[2],theta=MU[2,g],disp_param=size[2])
        cell_matrix[n,g]<- out_2$Y
        U_ng[n,g]<- out_2$U}
        
    }
  }
  
  #cell_matrix[1:10,1:5]
  
  #U_ng[1:10,1:5]
  #print(cell_matrix)
  #print(U_ng)
  #result_test<- ZINB_mix_K_comp(Y=cell_matrix,P=P,Phi=Phi,mu=MU,dispresiopn_parameter=size,eps=1e-5, maxit=1000)
  elapsed[s]<- system.time({result_test<- ZINB_mix_K_comp(Y=cell_matrix,P=P,Phi=Phi,U_ng=U_ng,mu=MU,clust=clust,dispresiopn_parameter=size,eps=1e-5, maxit=1000)})["elapsed"]
  P_hat[s,]<- result_test$P
  Phi_hat[s,]<- result_test$Phi
  mu1_hat[s,]<- result_test$MU[1,]
  mu2_hat[s,]<- result_test$MU[2,]
  overdisp_hat[s,]<- result_test$a_hat
  size_hat[s,]<- result_test$size_est
  # Z0<- clust
  # zhat<- result_test$postprobs
  # zhat<- apply(zhat,1,which.max)
  # V_M[s]<- v_measure(Z0,zhat)
  
  
  
}



###############################################################
L1<- cbind(colMeans(mu1_hat),mu_1)
colnames(L1)<- c("mu1_hat","mu1")
#print(L1)
L2<- cbind(colMeans(mu2_hat),mu_2)
colnames(L2)<- c("mu2_hat","mu2")
#print(L2)

dif_mu1<- dif_mu2<- matrix(NA,nrow=S,ncol=G)
mad_mu1<- mad_mu2<- matrix(NA,nrow=S,ncol=G)

for(s in 1:S){
  dif_mu1[s,]<- (mu1_hat[s,]-mu_1)^2
  dif_mu2[s,]<- (mu2_hat[s,]-mu_2)^2
  mad_mu1[s,]<- median(abs(mu1_hat-mu_1))
  mad_mu2[s,]<- median(abs(mu2_hat-mu_2))
  
  
  
}

#dif_mu1
#dif_mu2

MSE_mu1<- apply(dif_mu1,2,mean)
#MSE_mu1
#mean(MSE_mu1)

MSE_mu2<- apply(dif_mu2,2,mean)

#MSE_mu2
MSE_MU<- rbind(mean(MSE_mu1),mean(MSE_mu2))
#MSE_MU
#xtable(MSE_MU)
#mean(MSE_mu1)

#########################################################################################


P_sim1_n60<- write.csv(P_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/P_sim1_n60.csv")
Phi_sim1_n60<- write.csv( Phi_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/Phi_sim1_n60.csv")
overdisp_sim1_n60<- write.csv(overdisp_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/overdisp_sim1_n60.csv")
size_sim1_n60<- write.csv(size_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/size_sim1_n60.csv")
MU1_sim1_n60<- write.csv(MSE_mu1,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/MU1_sim1_n60.csv")
MU2_sim1_n60<- write.csv(MSE_mu2,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/MU2_sim1_n60.csv")
V_M_sim1_n60<- write.csv(V_M,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/v_M_sim1_n60.csv")
elapsed_sim1_n60<- write.csv(elapsed,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/elapsed_sim1_n60.csv")


```




``` {r Scenario1-sim2,eval=TRUE,echo=FALSE}
rm(list = ls())

###################################################################################################
### Generate data from Zero-inflated Negative Binomial Distribution### 
ZINB_r<- function(n,prob,theta,disp_param){
  Z<- rbinom(size=1,n=n,p=prob)
  Y<- NULL
  for(i in 1:n){
    if(Z[i]==1) Y<- c(Y,0)
    else Y<- c(Y,rnbinom(1,mu = theta,size=disp_param))
  }
  return(Y)
}
#ZINB_r(50,0.5,4,10)
#################################################################################################################
#################################################################################################################

##################################################################################################################
###################################################################################################################
disp_para<- function (y, mu, n = sum(w), w, limit = 500, eps = .Machine$double.eps^0.25, trace = FALSE) {
  
  score <- function(n, th, mu, y, w) sum(w * (digamma(th +y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +th)/(mu + th)))
  info <- function(n, th, mu, y, w) sum(w * (-trigamma(th +y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +th)^2))
  if (inherits(y, "lm")) {
    mu <- y$fitted.values
    y <- if (is.null(y$y)) 
      mu + residuals(y)
    else y$y
  }
  #if (missing(weights)) 
  if (missing(w)) 
    
    #  weights <- rep(1, length(y))
    w <- rep(1, length(y))
  
  #t0 <- n/sum(weights * (y/mu - 1)^2)
  t0 <- n/sum(w * (y/mu - 1)^2)
  
  it <- 0
  del <- 1
  if (trace) 
    message(sprintf("theta.ml: iter %d 'theta = %f'", it, 
                    signif(t0)), domain = NA)
  while ((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    #del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, mu, y, weights))
    
    del <- score(n, t0, mu, y, w)/(i <- info(n, t0, mu, y, w))
    t0 <- t0 + del
    if (trace) 
      message("theta.ml: iter", it, " theta =", signif(t0))
  }
  if (t0 < 0) {
    t0 <- 0
    warning("estimate truncated at zero")
    attr(t0, "warn") <- gettext("estimate truncated at zero")
  }
  if (it == limit) {
    warning("iteration limit reached")
    attr(t0, "warn") <- gettext("iteration limit reached")
  }
  attr(t0, "SE") <- sqrt(1/i)
  t0
}
#######################################################################################################
#######################################################################################################
To_vec_transfer<- function(Y,Z,U,MU){
  N<- nrow(Y)
  G<- ncol(Y)
  K<- ncol(Z)
  
  vec_y<- as.vector(t(Y))
  ZU<- array(NA,c(N,G,K))
  for(k in 1:K){
    ZU[,,k]<- Z[,k]*(1-U)
  }
  
  
  mu_mat<- matrix(NA,nrow=N*G,ncol=K)
  for(k in 1:K){
    mu_mat[,k]<- rep(MU[k,],N)
  }
  
  weight_mat<- matrix(NA,nrow=N*G,ncol=K)
  for(k in 1:K){
    weight_mat[,k]<- as.vector(t(ZU[,,k]))
  }
  
  out<- list(Y_vec=vec_y,WW=weight_mat,MU_vec=mu_mat)
  return(out)
  
} 


########################################################
##### THE EM Function for Mixture of MIZNB####
###################################################################################################
ZINB_mix_K_comp<- function(Y,P,Phi,mu,dispresiopn_parameter,eps=1e-5, maxit=1000){
  dl <- 1 + eps
  iter<-0
  ll <- rep(0, maxit+1)
  G<- ncol(Y)
  N<- nrow(Y)
  K<- length(Phi)
  
  nu<- dispresiopn_parameter
  
  bk<- array(NA,c(N,G,K))
  for(k in 1:K){
    for(n in 1:N){
      for(g in 1:G){
        if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k]))) 
        
        
        else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])
      }
    }
  }
  
  
  b<- array(NA,c(N,G,K))
  
  for(k in 1:K) 
  {b[,,k]<- P[k]*bk[,,k]}
  
  bp_sum<- apply(b,c(1,2),sum)
  l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
  Z_nk<- matrix(NA,ncol=K,nrow=N)
  thresh<- -744
  
  
  U_ng<- matrix(data=NA,nrow=N,ncol=G)
  ZU<- array(data=NA,c(N,G,K))
  
  
  
  d<- array(data=NA,c(N,G,K))
  
  for(n in 1:N){
    for(g in 1:G){
      for(k in 1:K){
        if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))
        
        else d[n,g,k]<- 0
        
      }
    }
  }
  
  
  
  
  
  while (abs(dl)>eps & iter<maxit) { ## Here I use abs(dl) rather than dl#
    iter<-iter+1
    ll[iter] <- l
    
    
    # Z_nk<- matrix(data=log(P),ncol=K,nrow=N,byrow=T)
    # 
    # for(k in 1:K){
    #   for(g in 1:G){
    #     Z_nk[,k]<- Z_nk[,k]+(log(bk[,g,k]))
    #   }
    # }
    # 
    # if(K>1){
    #   v1<- which(apply(Z_nk,1,max)< thresh)
    #   v3<- 1:N
    #   len<- length(v1)
    #   #print(len)
    #   if(len>0){
    #     v2<- apply(array(Z_nk[v1,],dim=c(len,K)),1,order)[K,]
    #     ddd<- cbind(v1,v2)
    #     #print(ddd)
    #     Z_nk[v1,]<- 0
    #     Z_nk[ddd]<- 1
    #     v3<- -v1
    #   }
    #   
    #   Z_nk[v3,]<- exp(Z_nk[v3,])
    # }
    # Z_nk<- Z_nk/apply(Z_nk,1,sum)
    # #print(Z_nk)
    # 
    # epsilon<- 1e-10
    # sl<- length(Z_nk[Z_nk<epsilon])
    # bl<- length(Z_nk[Z_nk>1-epsilon])
    # Z_nk[Z_nk<epsilon]<- rep(epsilon,sl)
    # Z_nk[Z_nk>1-epsilon]<- rep(1-epsilon,bl)
    # Z_nk<- Z_nk/apply(Z_nk,1,sum)
    #print(Z_nk)
    
    
    #Z_nk#
    Z_nk<- matrix(0,nrow=N,ncol=length(P))
    for(i in 1:N){
      if(clust[i]==1)Z_nk[i,1]<- 1
      else Z_nk[i,2]<- 1
    }
    #Z
    ##################################################
    ##################################################
    ## U_ng##
    U_ng<- matrix(NA,nrow=N,ncol=G)
    for(n in 1:N){
      for(g in 1:G){
        if(cell_matrix[n,g]==0) U_ng[n,g]<- 1
        else U_ng[n,g]<-0
      }
    }
    
    
    P<- apply(Z_nk,2,mean)
    #print(P)
    
    
    # U_ng<- matrix(data=NA,nrow=N,ncol=G)
    # 
    # denom<- apply(d,c(1,2),sum)
    # 
    # for(n in 1:N){
    #   for(g in 1:G){
    #     if(Y[n,g]==0) U_ng[n,g]<- sum(P*Phi)/denom[n,g]
    #     else U_ng[n,g]<- 0
    #     
    #   }
    # }
    
    
    for(k in 1:K){
      ZU[,,k]<- Z_nk[,k]*U_ng
    }
    
    Phi<- (apply(ZU,3,sum))/(G*apply(Z_nk,2,sum))
    ###############################################################################################
    ################################################################################################
    
    
    # mult_Z_u<- array(NA,c(N,G,K))
    # mult_Z_oneminus_U<- array(NA,c(N,G,K))
    # mult_Z_oneminus_U_Y<- array(NA,c(N,G,K))
    # 
    # 
    # 
    # for(k in 1:K){
    #   mult_Z_u[,,k]<- Z_nk[,k]*U_ng
    #   mult_Z_oneminus_U[,,k]<- Z_nk[,k]*(1-U_ng)
    #   mult_Z_oneminus_U_Y[,,k]<- Z_nk[,k]*(1-U_ng)*Y
    #   
    # }
    # 
    # 
    # mu<- t(apply(mult_Z_oneminus_U_Y,c(2,3),sum)/apply(mult_Z_oneminus_U,c(2,3),sum))
    mu<- mu
    ##########################################################################################
    ##########################################################################################
    ### To estimate size_MLE###
    result<- To_vec_transfer(Y=cell_matrix,Z=Z_nk,U=U_ng,MU=mu)
    
    size_MLE<- rep(NA,K)
    for(k in 1:K){
      size_MLE[k]<- disp_para(y=result$Y_vec, mu=result$MU_vec[,k], n = sum(result$WW[,k]), w=result$WW[,k], limit = 500, eps = .Machine$double.eps^0.25, trace = FALSE) 
      
    }
    
    nu<- size_MLE
    
    ################################################################################################################
    ################################################################################################################
    for(k in 1:K){
      for(n in 1:N){
        for(g in 1:G){
          if(Y[n,g]==0) bk[n,g,k]<- (Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k]))) 
          else bk[n,g,k]<- (1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])
        }
      }
    }
    
    ############################################################
    for(n in 1:N){
      for(g in 1:G){
        for(k in 1:K){
          if(Y[n,g]==0) d[n,g,k]<- P[k]*(Phi[k]+((1-Phi[k])*dnbinom(x=Y[n,g],mu=mu[k,g],size=nu[k])))
          else d[n,g,k]<- 0
          
        }
      }
    }
    
    #########################################################
    
    
    
    for(k in 1:K){
      b[,,k]<- P[k]*bk[,,k]
    }
    bp_sum<- apply(b,c(1,2),sum)
    #print(bp_sum)
    oldl<- l
    
    l<- sum(rowSums(log(bp_sum))) #Observed data log-likelihood
    #print(l)
    
    ################################################################
    
    #print(l)
    dl<- l-oldl
  }
  cat("number of iterations=", iter, "\n")
  iter <- iter+1
  ll[iter] <- l
  postprobs <- Z_nk
  colnames(postprobs) <- c(paste("comp", ".", 1:K, sep = ""))
  
  
  out <- list(P=P,Phi=Phi,MU=mu,size_est=nu,a_hat=1/nu,postprobs=postprobs,
              loglik=l,all.loglik=ll[1:iter],restart=0,ft="ZINBmix")
  class(out)<- "mixEM"
  out
  
}
###################################################################################################
###################################################################################################
### Data sim_2###
### Generate Y Matrix_Matrix of raw counts###
### case_1###

N<- 120
G<- 120
P<- c(0.5,0.5)
Phi<- c(0.1,0.1)
K<- length(P)
#size<- c(10,5)
#size<- c(5,10)
size<- c(5,10)

#P<- c(0.25,0.25,0.25,0.25)
#N<- 1000
#G<- 10
#Phi<- c(0.1,0.1,0.2,0.2)
#P<- c(0.25,0.25,0.25,0.25)
#P1<- c(0.4,0.3,0.2,0.1)


#mu_1<- rgamma(n=G,shape=10,scale=1/2)
#mu_2<- rgamma(n=G,shape=10,scale=1/1.5)


#mu_1<- c(rep(5,(G/2)),rep(10,G/2))
#mu_2<- c(rep(10,(G/2)),rep(5,G/2))

mu_1<- c(rep(5,G))
mu_2<- c(rep(10,G))
MU<- rbind(mu_1,mu_2)
#MU

#S<- 256
#S<- 128
#S<- 50
S<- 3
P_hat<- matrix(NA,ncol=length(P),nrow=S)
#P_hat<- NULL
Phi_hat<- matrix(NA,ncol=length(P),nrow=S)
overdisp_hat<- matrix(NA,nrow=S,ncol=K)
size_hat<- matrix(NA,nrow=S,ncol=K)

#Phi_hat
#phi_hat2<- NULL
mu1_hat<- matrix(NA,ncol = G,nrow=S)
mu2_hat<- matrix(NA,ncol = G,nrow=S)
#lambda3_hat<- matrix(NA,ncol = G,nrow=S)
#lambda4_hat<- matrix(NA,ncol = G,nrow=S)

#MU_1<- MU_2<- matrix(NA,nrow=S,ncol=G)
V_M<- rep(NA,S)
elapsed<- rep(NA,S)

for(s in 1:S){
  
  # mu_1<- rgamma(n=G,shape=10,scale=1/2)
  
  
  #mu_2<- rgamma(n=G,shape=10,scale=1/0.2)
  
  #MU<- rbind(mu_1,mu_2)
  
  
  
  #for(s in 1:S){
  
  
  #for(iteration in 1:r){
  
  #p_hat<- NULL
  #phi_hat1<- NULL
  #phi_hat2<- NULL
  #lambda1_hat<- matrix(0,ncol = G,nrow=r)
  #lambda2_hat<- matrix(0,ncol = G,nrow=r)
  
  #for(iter in 1:r){
  #clust<- rbinom(size=1,n=N,p=0.5)
  clust<- sample(1:2,size=N,replace = TRUE,prob=P)
  #clust
  #length(clust[clust==1])
  
  
  #clust<- sample(c(1,2),size=N,replace = TRUE,prob=pn) 
  #clust
  
  
  cell_matrix<- matrix(rep(0,N*G),nrow=N,ncol=G)
  
  for(n in 1:N){
    for(g in 1:G){
      if(clust[n]==1) cell_matrix[n,g]<- ZINB_r(n=1,prob=Phi[1],theta=MU[1,g],disp_param=size[1])
      else cell_matrix[n,g]<- ZINB_r(n=1,prob=Phi[2],theta=MU[2,g],disp_param =size[2])
    }
  }
  
  
  
  
  
  
  
  #cell_matrix
  #cell_matrix
  
  #plot(hclust(dist(cell_matrix)))
  
  
  #cell_matrix
  #dim(cell_matrix) 
  
  
  #result_test<- ZINB_mix_K_comp(Y=cell_matrix,P=P,Phi=Phi,mu=MU,dispresiopn_parameter=size,eps=1e-5, maxit=1000)
  elapsed[s]<- system.time({result_test<- ZINB_mix_K_comp(Y=cell_matrix,P=P,Phi=Phi,mu=MU,dispresiopn_parameter=size,eps=1e-5, maxit=1000)})["elapsed"]
  P_hat[s,]<- result_test$P
  Phi_hat[s,]<- result_test$Phi
  mu1_hat[s,]<- result_test$MU[1,]
  mu2_hat[s,]<- result_test$MU[2,]
  overdisp_hat[s,]<- result_test$a_hat
  size_hat[s,]<- result_test$size_est
  Z0<- clust
  zhat<- result_test$postprobs
  zhat<- apply(zhat,1,which.max)
  V_M[s]<- v_measure(Z0,zhat)
  
  
  
}



###############################################################
L1<- cbind(colMeans(mu1_hat),mu_1)
colnames(L1)<- c("mu1_hat","mu1")
#print(L1)
L2<- cbind(colMeans(mu2_hat),mu_2)
colnames(L2)<- c("mu2_hat","mu2")
#print(L2)

dif_mu1<- dif_mu2<- matrix(NA,nrow=S,ncol=G)

for(s in 1:S){
  dif_mu1[s,]<- (mu1_hat[s,]-mu_1)^2
  dif_mu2[s,]<- (mu2_hat[s,]-mu_2)^2
  
}

#dif_mu1
#dif_mu2

MSE_mu1<- apply(dif_mu1,2,mean)
#MSE_mu1
#mean(MSE_mu1)

MSE_mu2<- apply(dif_mu2,2,mean)
#MSE_mu2
#mean(MSE_mu2)
#mu1_hat
#boxplot(mu1_hat)
#abline(h=5)
#abline(h=10)

#boxplot(mu2_hat)
#abline(h=5)
#abline(h=10)

#MSE_mu1
#MSE_mu2
MSE_MU<- rbind(mean(MSE_mu1),mean(MSE_mu2))
#MSE_MU
#xtable(MSE_MU)
#mean(MSE_mu1)

#########################################################################################


P_sim2_n120<- write.csv(P_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/P_sim2_n120.csv")
Phi_sim2_n120<- write.csv( Phi_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/Phi_sim2_n120.csv")
overdisp_sim2_n120<- write.csv(overdisp_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/overdisp_sim2_n120.csv")
size_sim2_n120<- write.csv(size_hat,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/size_sim2_n120.csv")
MU1_sim2_n120<- write.csv(MSE_mu1,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/MU1_sim2_n120.csv")
MU2_sim2_n120<- write.csv(MSE_mu2,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/MU2_sim2_n120.csv")
V_M_sim2_n120<- write.csv(V_M,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/v_M_sim2_n120.csv")
elapsed_sim2_n120<- write.csv(elapsed,"E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario1/elapsed_sim2_n120.csv")


```
