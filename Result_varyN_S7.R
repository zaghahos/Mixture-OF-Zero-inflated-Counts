rm(list=ls())

#####################################################################################################################
P_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/P_sim1_n60.csv")
Phi_sim1_n60 <- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Phi_sim1_n60.csv")
overdisp_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/overdisp_sim1_n60.csv")
size_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/size_sim1_n60.csv")
MSE_MU1_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU1_sim1_n60.csv")
MSE_MU2_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU2_sim1_n60.csv")
V_M_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/v_M_sim1_n60.csv")
elapsed_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/elapsed_sim1_n60.csv")
#Med_MU1_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU1_sim1_n60.csv")
#Med_MU2_sim1_n60<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU2_sim1_n60.csv")

######################################################################################################################################
P_sim1_n60<- P_sim1_n60[,-1]
Phi_sim1_n60<- Phi_sim1_n60[,-1]
overdisp_sim1_n60<- overdisp_sim1_n60[,-1]
size_sim1_n60<- size_sim1_n60[,-1]
MSE_MU1_sim1_n60<- MSE_MU1_sim1_n60[,-1]
MSE_MU2_sim1_n60<- MSE_MU2_sim1_n60[,-1]
V_M_sim1_n60<- V_M_sim1_n60[,-1]
elapsed_sim1_n60<- elapsed_sim1_n60[,-1]
#Med_MU1_sim1_n60<- Med_MU1_sim1_n60[,-1]
#Med_MU2_sim1_n60<- Med_MU2_sim1_n60[,-1]
#######################################################################################################################################
P_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/P_sim2_n120.csv")
Phi_sim2_n120 <- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Phi_sim2_n120.csv")
overdisp_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/overdisp_sim2_n120.csv")
size_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/size_sim2_n120.csv")
MSE_MU1_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU1_sim2_n120.csv")
MSE_MU2_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU2_sim2_n120.csv")
V_M_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/v_M_sim2_n120.csv")
elapsed_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/elapsed_sim2_n120.csv")
Med_MU1_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU1_sim2_n120.csv")
Med_MU2_sim2_n120<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU2_sim2_n120.csv")
######################################################################################################################################
P_sim2_n120<- P_sim2_n120[,-1]
Phi_sim2_n120<- Phi_sim2_n120[,-1]
overdisp_sim2_n120<- overdisp_sim2_n120[,-1]
size_sim2_n120<- size_sim2_n120[,-1]
MSE_MU1_sim2_n120<- MSE_MU1_sim2_n120[,-1]
MSE_MU2_sim2_n120<- MSE_MU2_sim2_n120[,-1]
V_M_sim2_n120<- V_M_sim2_n120[,-1]
elapsed_sim2_n120<- elapsed_sim2_n120[,-1]
Med_MU1_sim2_n120<- Med_MU1_sim2_n120[,-1]
Med_MU2_sim2_n120<- Med_MU2_sim2_n120[,-1]
#######################################################################################################################################
P_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/P_sim3_n300.csv")
Phi_sim3_n300 <- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Phi_sim3_n300.csv")
overdisp_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/overdisp_sim3_n300.csv")
size_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/size_sim3_n300.csv")
MSE_MU1_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU1_sim3_n300.csv")
MSE_MU2_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU2_sim3_n300.csv")
V_M_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/v_M_sim3_n300.csv")
elapsed_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/elapsed_sim3_n300.csv")
Med_MU1_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU1_sim3_n300.csv")
Med_MU2_sim3_n300<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU2_sim3_n300.csv")
######################################################################################################################################
P_sim3_n300<- P_sim3_n300[,-1]
Phi_sim3_n300<- Phi_sim3_n300[,-1]
overdisp_sim3_n300<- overdisp_sim3_n300[,-1]
size_sim3_n300<- size_sim3_n300[,-1]
MSE_MU1_sim3_n300<- MSE_MU1_sim3_n300[,-1]
MSE_MU2_sim3_n300<- MSE_MU2_sim3_n300[,-1]
V_M_sim3_n300<- V_M_sim3_n300[,-1]
elapsed_sim3_n300<- elapsed_sim3_n300[,-1]
Med_MU1_sim3_n300<- Med_MU1_sim3_n300[,-1]
Med_MU2_sim3_n300<- Med_MU2_sim3_n300[,-1]

#######################################################################################################################################
#######################################################################################################################################

P_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/P_sim4_n600.csv")
Phi_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Phi_sim4_n600.csv")
overdisp_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/overdisp_sim4_n600.csv")
size_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/size_sim4_n600.csv")
MSE_MU1_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU1_sim4_n600.csv")
MSE_MU2_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU2_sim4_n600.csv")
V_M_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/v_M_sim4_n600.csv")
elapsed_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/elapsed_sim4_n600.csv")
Med_MU1_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU1_sim4_n600.csv")
Med_MU2_sim4_n600<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU2_sim4_n600.csv")

P_sim4_n600<- P_sim4_n600[,-1]
Phi_sim4_n600<- Phi_sim4_n600[,-1]
overdisp_sim4_n600<- overdisp_sim4_n600[,-1]
size_sim4_n600<- size_sim4_n600[,-1]
MSE_MU1_sim4_n600<- MSE_MU1_sim4_n600[,-1]
MSE_MU2_sim4_n600<- MSE_MU2_sim4_n600[,-1]
V_M_sim4_n600<- V_M_sim4_n600[,-1]
elapsed_sim4_n600<- elapsed_sim4_n600[,-1]

Med_MU1_sim4_n600<- Med_MU1_sim4_n600[,-1]
Med_MU2_sim4_n600<- Med_MU2_sim4_n600[,-1]

#############################################################################################################################



P_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/P_sim5_n1200.csv")
Phi_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Phi_sim5_n1200.csv")
overdisp_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/overdisp_sim5_n1200.csv")
size_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/size_sim5_n1200.csv")
MSE_MU1_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU1_sim5_n1200.csv")
MSE_MU2_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/MSE_MU2_sim5_n1200.csv")
V_M_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/v_M_sim5_n1200.csv")
elapsed_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/elapsed_sim5_n1200.csv")
Med_MU1_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU1_sim5_n1200.csv")
Med_MU2_sim5_n1200<- read.csv("E:/Thesis-Feb,March2022/Programming Project 2/ZINB_A/Scenario7/Med_MU2_sim5_n1200.csv")

P_sim5_n1200<- P_sim5_n1200[,-1]
Phi_sim5_n1200<- Phi_sim5_n1200[,-1]
overdisp_sim5_n1200<- overdisp_sim5_n1200[,-1]
size_sim5_n1200<- size_sim5_n1200[,-1]
MSE_MU1_sim5_n1200<- MSE_MU1_sim5_n1200[,-1]
MSE_MU2_sim5_n1200<- MSE_MU2_sim5_n1200[,-1]
V_M_sim5_n1200<- V_M_sim5_n1200[,-1]
elapsed_sim5_n1200<- elapsed_sim5_n1200[,-1]

Med_MU1_sim5_n1200<- Med_MU1_sim5_n1200[,-1]
Med_MU2_sim5_n1200<- Med_MU2_sim5_n1200[,-1]
###########################################################################################################################
#size_sim1_n60
#size_sim2_n120
#size_sim3_n300
#size_sim5_n1200
#############################################################################################################################
par(mfrow=c(1,2))
boxplot(size_sim1_n60[,1],size_sim2_n120[,1],size_sim3_n300[,1],size_sim4_n600[,1],size_sim5_n1200[,1],main=expression(nu[1]),cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
abline(h=5)

boxplot(size_sim1_n60[,2],size_sim2_n120[,2],size_sim3_n300[,2],size_sim4_n600[,2],size_sim5_n1200[,2],main=expression(nu[2]),cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
abline(h=20)

#boxplot(size_sim1_n60[,1],size_sim2_n120[,1],size_sim3_n300[,1],size_sim4_n600[,1],size_sim5_n1200[,1],main=expression(nu[1]),outline =FALSE,cex=0.5,xlab="N_rm_oultiler",names=c(60,120,300,600,1200))
#abline(h=5)

#boxplot(size_sim1_n60[,2],size_sim2_n120[,2],size_sim3_n300[,2],size_sim4_n600[,2],size_sim5_n1200[,2],main=expression(nu[2]),outline = FALSE,cex=0.5,xlab="N_rm_outlier",names=c(60,120,300,600,1200))
#abline(h=20)

#################################################################################################################################
par(mfrow=c(1,2))
boxplot(P_sim1_n60[,1],P_sim2_n120[,1],P_sim3_n300[,1],P_sim4_n600[,1],P_sim5_n1200[,1],main=expression(pi[1]),cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
abline(h=0.5)

boxplot(P_sim1_n60[,2],P_sim2_n120[,2],P_sim3_n300[,2],P_sim4_n600[,2],P_sim5_n1200[,2],main=expression(pi[2]),cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
abline(h=0.5)

############################################################################################################################
par(mfrow=c(1,2))
boxplot(Phi_sim1_n60[,1],Phi_sim2_n120[,1],Phi_sim3_n300[,1],Phi_sim4_n600[,1],Phi_sim5_n1200[,1],main=expression(phi[1]),cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
abline(h=0.1)

boxplot(Phi_sim1_n60[,2],Phi_sim2_n120[,2],Phi_sim3_n300[,2],Phi_sim4_n600[,2],Phi_sim5_n1200[,2],main=expression(phi[2]),cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
abline(h=0.1)

############################################################################################################################
#par(mfrow=c(1,2))
boxplot(V_M_sim1_n60,V_M_sim2_n120,V_M_sim3_n300,V_M_sim4_n600,V_M_sim5_n1200,main="V-Measures",cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
abline(h=1)

boxplot(elapsed_sim1_n60,elapsed_sim2_n120,elapsed_sim3_n300,elapsed_sim4_n600,elapsed_sim5_n1200,main="Computing time",cex.axis=0.5,xlab="N",names=c(60,120,300,600,1200))
#################################################################################################################################
###########################################################################################################################################
M_P1<- c(mean(P_sim1_n60[,1]),mean(P_sim2_n120[,1]),mean(P_sim3_n300[,1]),mean(P_sim4_n600[,1]),mean(P_sim5_n1200[,1]))
M_P1
M_P2<- c(mean(P_sim1_n60[,2]),mean(P_sim2_n120[,2]),mean(P_sim3_n300[,2]),mean(P_sim4_n600[,2]),mean(P_sim5_n1200[,2]))
M_P2

SE_P1<- c(sd(P_sim1_n60[,1]),sd(P_sim2_n120[,1]),sd(P_sim3_n300[,1]),sd(P_sim4_n600[,1]),sd(P_sim5_n1200[,1]))
SE_P1

SE_P2<- c(sd(P_sim1_n60[,2]),sd(P_sim2_n120[,2]),sd(P_sim3_n300[,2]),sd(P_sim4_n600[,2]),sd(P_sim5_n1200[,2]))
SE_P2


K<- c(rep(1,5),rep(2,5))

N<- rep(c(60,120,300,600,1200),2)

P_est_N<- cbind(K,N,rbind(cbind(M_P1,SE_P1),cbind(M_P2,SE_P2)))
P_est_N
xtable(P_est_N,digits = 9)
########################################
### Phi###
M_Phi1<- c(mean(Phi_sim1_n60[,1]),mean(Phi_sim2_n120[,1]),mean(Phi_sim3_n300[,1]),mean(Phi_sim4_n600[,1]),mean(Phi_sim5_n1200[,1]))
M_Phi1
M_Phi2<- c(mean(Phi_sim1_n60[,2]),mean(Phi_sim2_n120[,2]),mean(Phi_sim3_n300[,2]),mean(Phi_sim4_n600[,2]),mean(Phi_sim5_n1200[,2]))
M_Phi2

SE_Phi1<- c(sd(Phi_sim1_n60[,1]),sd(Phi_sim2_n120[,1]),sd(Phi_sim3_n300[,1]),sd(Phi_sim4_n600[,1]),sd(Phi_sim5_n1200[,1]))
SE_Phi1

SE_Phi2<- c(sd(Phi_sim1_n60[,2]),sd(Phi_sim2_n120[,2]),sd(Phi_sim3_n300[,2]),sd(Phi_sim4_n600[,2]),sd(Phi_sim5_n1200[,2]))
SE_Phi2


K<- c(rep(1,5),rep(2,5))

N<- rep(c(60,120,300,600,1200),2)

Phi_est_N<- cbind(K,N,rbind(cbind(M_Phi1,SE_Phi1),cbind(M_Phi2,SE_Phi2)))
Phi_est_N
xtable(Phi_est_N,digits = 9)

# ### Size###
M_size_K1<- c(mean(size_sim1_n60[,1])
              ,mean(size_sim2_n120[,1])
              ,mean(size_sim3_n300[,1])
              ,mean(size_sim4_n600[,1])
              ,mean(size_sim5_n1200[,1]))

M_size_K2<- c(mean(size_sim1_n60[,2])
              ,mean(size_sim2_n120[,2])
              ,mean(size_sim3_n300[,2])
              ,mean(size_sim4_n600[,2])
              ,mean(size_sim5_n1200[,2]))

SE_size_K1<- c(sd(size_sim1_n60[,1])
               ,sd(size_sim2_n120[,1])
               ,sd(size_sim3_n300[,1])
               ,sd(size_sim4_n600[,1])
               ,sd(size_sim5_n1200[,1]))

SE_size_K2<- c(sd(size_sim1_n60[,2])
               ,sd(size_sim2_n120[,2])
               ,sd(size_sim3_n300[,2])
               ,sd(size_sim4_n600[,2])
               ,sd(size_sim5_n1200[,2]))



size_N<- cbind(K,N,rbind(cbind(M_size_K1,SE_size_K1),cbind(M_size_K2,SE_size_K2)))
xtable(size_N)
size_N

#########################################################################################################################
MU_N_K1<- c(mean(MSE_MU1_sim1_n60),mean(MSE_MU1_sim2_n120),mean(MSE_MU1_sim3_n300),mean(MSE_MU1_sim4_n600),mean(MSE_MU1_sim5_n1200)) 
MU_N_K2<- c(mean(MSE_MU2_sim1_n60),mean(MSE_MU2_sim2_n120),mean(MSE_MU2_sim3_n300),mean(MSE_MU2_sim4_n600),mean(MSE_MU2_sim5_n1200)) 
clust<- c(1,2)
MU_N_K1
MU_N_K2

MSE_MU<- rbind(MU_N_K1,MU_N_K2)
MSE_MU


colnames(MSE_MU)<- c(60,120,300,600,1200)
MSE_MU
res_MSE_MU<- cbind(clust,MSE_MU)

res_MSE_MU

xtable(res_MSE_MU,digits = 5)

#MU_G_K1_Med<- c(median(Med_MU1_sim1_G12),median(Med_MU1_sim2_G60),median(Med_MU1_sim3_G120),median(Med_MU1_sim4_G600),median(Med_MU1_sim5_G1500))
#MU_G_K1_Med

#MU_G_K2_Med<- c(median(Med_MU2_sim1_G12),median(Med_MU2_sim2_G60),median(Med_MU2_sim3_G120),median(Med_MU2_sim4_G600),median(Med_MU2_sim5_G1500))
#MU_G_K2_Med


#Med_MU<- rbind(MU_G_K1_Med,MU_G_K2_Med)
#Med_MU


#colnames(Med_MU)<- c(12,60,120,600,1500)
#Med_MU
#res_Med_MU<- cbind(clust,Med_MU)

#res_Med_MU
##############################################################################################################
M_time<- c(mean(elapsed_sim1_n60),mean(elapsed_sim2_n120),mean(elapsed_sim3_n300),mean(elapsed_sim4_n600),mean(elapsed_sim5_n1200))
SE_time<- c(sd(elapsed_sim1_n60),sd(elapsed_sim2_n120),sd(elapsed_sim3_n300),sd(elapsed_sim4_n600),sd(elapsed_sim5_n1200))

N1<- c(60,120,300,600,1200)
comp_time<- cbind(N1,M_time,SE_time)
xtable(comp_time)
