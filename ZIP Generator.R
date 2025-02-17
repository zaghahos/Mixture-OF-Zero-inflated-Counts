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

myvec<- ZIP_r(n=10000,pi=0.2,theta=20)
hist(myvec)


