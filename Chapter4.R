#install.packages(c('tabulizer','copula','pracma','CASdatasets','SurvCorr'))
library(tabulizer)
library(MASS) # for ginv
library(dplyr)
library(copula) # for claytonCopula
library(pracma) # for NewtonRaphson procedure
library(CASdatasets) # Canadian life insurance dataset # nit available for some R versions
library(SurvCorr) # kidney dataset
library(survival) # for Kaplan-Meier function

install.packages('mhazard')
library(mhazard) # for Dabrowska (1988) estimator for bivariate survivor


# Installing copula package in colab may require 'gsl'. To install gsl in colab:
# https://stackoverflow.com/questions/64009256/rcppgsl-install-fails-on-google-colab-r-notebook
# 1. open PyThon notebook in Colab Google
# 2. Run the following set of commands:
# # !sudo apt install libgsl-dev
# # !pip install rpy2
# # %reload_ext rpy2.ipython
# 3. Continue working in the Python notebook, by adding %%R at the beginning of each cell

# # # Functions


# # - - - - - -  - - - - - - - - - - - - - - - 
# # - - - - - -  - - - - - - - - - - - - - - - 



# Dabrowska (1988) bivariate survival estimator, for comparison of my work
# in 4.1.3 and Wang and Wells (2000)

#install.packages(c('mhazard','SurvCorr'))
#update.packages('Rcpp') # need to update 'Rcpp'

library(SurvCorr) # kidney dataset

data(kidney)
model_dab <- KM2(Y1=kidney$TIME1,Y2=kidney$TIME2,Delta1=kidney$STATUS1,Delta2=kidney$STATUS2,estimator="dabrowska")

# Function giving Wang and Wells gamma (tau with ties) estimates from a data set

GammaWangWellsFunc <- function(dataset){
  
  # Call example: GammaWangWellsFunc(kidney)
  
  require(mhazard) # for Dabrowska (1988) estimator for bivariate survivor
  
  # Dabrowska estimator, at grid of the observed data
  
  T1.unique <- unique(c(0,dataset$TIME1))[order(unique(c(0,dataset$TIME1)))]
  T2.unique <- unique(c(0,dataset$TIME2))[order(unique(c(0,dataset$TIME2)))]
  F_dab <- matrix(NA,nrow=length(T1.unique),ncol=length(T2.unique))
  
  model_dab <- KM2(Y1=mydata$TIME1,Y2=mydata$TIME2,Delta1=mydata$STATUS1,Delta2=mydata$STATUS2,estimator="dabrowska")
  
  for(i in 1:length(T1.unique)){
    for(j in 1:length(T2.unique)){
      
      pos.i = 1+min(which(model_dab$T1 >= T1.unique[i]))
      pos.j = 1+min(which(model_dab$T2 >= T1.unique[j]))
      
      if(is.finite(pos.i) & is.finite(pos.j)){
        F_dab[i,j] <- model_dab$Fhat[pos.i,pos.j]
      }else{
        if(is.infinite(pos.i)){pos.i = max(mydata$TIME1)+1}
        if(is.infinite(pos.j)){pos.j = max(mydata$TIME2)+1}
        
        F_dab[i,j] <- KM2(Y1=mydata$TIME1,Y2=mydata$TIME2,Delta1=mydata$STATUS1,Delta2=mydata$STATUS2,newT1=pos.i,newT2=pos.j,estimator="dabrowska")$Fhat_est
        
      }
      #F_dab[i,j] <- model_dab$Fhat[1+min(which(model_dab$T1 >= T1.unique[i])),1+min(which(model_dab$T2 >= T2.unique[j]))]
    }
  }
  
  F_dab[1,1] <- 1 # needed? shouldn't be automatic from the computations?
  
  # Masses, at grid of the observed data
  
  F_dab_mass <- matrix(NA,nrow=dim(F_dab)[1],ncol=dim(F_dab)[2])
  F_dab_mass[1,] <- F_dab[1,]
  F_dab_mass[,1] <- F_dab[,1]
  
  for(i in 2:dim(F_dab)[1]){
    for(j in 2:dim(F_dab)[2]){
      F_dab_mass[i,j] <- F_dab[i,j]-F_dab[i-1,j]-F_dab[i,j-1]+F_dab[i-1,j-1]
    }
  }
  
  # Kendall's tau (tau_tilde of Wang and Wells)
  
  tau_tilde <- 4*sum(F_dab_mass[-1,-1] * F_dab[-1,-1])-1
  tau_tilde
  
  # Kendall's tau modified for ties (gamma of Wang and Wells)
  
  # ties sets
  O1 <- unique(dataset[duplicated(dataset$TIME1) & dataset$STATUS1==1,"TIME1"])
  O2 <- unique(dataset[duplicated(dataset$TIME2) & dataset$STATUS2==1,"TIME2"])
  
  if(length(O1)*length(O2)!=0){# if ties on T1 and T2. Need to rewrite this later.
    
    # Probability of ties
    
    Pr1 <- rep(NA,times=length(O1))
    for(l in 1:length(O1)){
      Pr1[l] <- (F_dab[which(T1.unique==O1[l])-1,1]-F_dab[which(T1.unique==O1[l]),1])**2
    }
    
    Pr2 <- rep(NA,times=length(O2))
    for(m in 1:length(O2)){
      Pr2[m] <- (F_dab[which(T2.unique==O2[m])-1,1]-F_dab[which(T2.unique==O2[m]),1])**2
    }
    
    Pr12 <- matrix(NA,nrow=length(O1),ncol=length(O2))
    for(l in 1:length(O1)){
      for(m in 1:length(O2)){
        Pr12[l,m] <- F_dab_mass[which(T1.unique==O1[l]),which(T2.unique==O2[m])]**2
      }
    }
    
    ProbTie <- sum(Pr1) + sum(Pr1) - sum(Pr12)
    
    gamma <- tau_tilde/(1-ProbTie)
  }else{gamma <- tau_tilde}
  
  return(gamma)
  
}


# # - - - - - -  - - - - - - - - - - - - - - - 
# # - - - - - -  - - - - - - - - - - - - - - - 


FunZDel = function(n,theta,lambda=1,seed,alpha=2,beta=2){
  
  require(copula)
  
  MyCopula <- mvdc(copula=claytonCopula(param=theta), # Clayton copula for (F(X), F(Y))
                   margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                   paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale
  
  set.seed(seed)
  XY <- rMvdc(n=n,MyCopula)
  X <- XY[,1]
  Y <- XY[,2]
  
  # Censoring variable

  Cx = rexp(n=n,rate = lambda)
  Cy = rexp(n=n,rate = lambda)
  
  # Observations
  
  Z1 <- pmin(X,Cx)
  Z2 <- pmin(Y,Cy)
  
  ordre = order(Z1,Z2)
  
  xinf=max(Z1,Z2)+1 # point at infinity
  Z1=c(Z1[ordre],xinf)
  Z2=c(Z2[ordre],xinf)
  
  del1 = as.integer(X <= Cx)
  del2 = as.integer(Y <= Cy)
  
  del1=c(del1,1)
  del2=c(del2,1)
  del=del1*del2
  del=del[ordre]
  
  return(list(Z1,Z2,del))
}


# Log-pseudo-likelihood function
#LogL = function(x){
#  sum(phat*(log(x+1)-(x+1)*log(Fn1*Fn2)-(2+1/x)*log(Fn1**(-x)+Fn2**(-x)-1)))
#}


FunZDel_canlifins = function(size,seed=1){
  
  data(canlifins) # load the dataset
  # 14,889 contracts where one annuitant is male and the other female
  
  canlifins = canlifins[(canlifins$EntryAgeM>=18)|(canlifins$EntryAgeF>=18),]
  
  set.seed(seed)
  Sample = sample(1:dim(canlifins)[1],size = size)
  canlifins = canlifins[Sample,]
  
  #  c(
  #    sum((canlifins$DeathTimeM == 0) & (canlifins$DeathTimeF == 0)), # number of doubly censored couples
  #    sum((canlifins$DeathTimeM > 0) & (canlifins$DeathTimeF == 0)), # number of couples where only the woman is censored
  #    sum((canlifins$DeathTimeM == 0) & (canlifins$DeathTimeF > 0)), # number of couples where only the man is censored
  #    sum((canlifins$DeathTimeM > 0) & (canlifins$DeathTimeF > 0)) # number of doubly uncensored couples
  #  )
  
  #del1 = as.integer((canlifins$DeathTimeM > 0) & (canlifins$DeathTimeF == 0)) 
  #del2 = as.integer((canlifins$DeathTimeM == 0) & (canlifins$DeathTimeF > 0))
  
  #del1=c(del1,1)
  #del2=c(del2,1)
  #del=del1*del2
  
  del = (canlifins$DeathTimeM > 0) * (canlifins$DeathTimeF > 0)

  ordre = order(canlifins$DeathTimeM+canlifins$EntryAgeM,canlifins$DeathTimeF+canlifins$EntryAgeF)
  Z_ordered = cbind(canlifins$DeathTimeM+canlifins$EntryAgeM,canlifins$DeathTimeF+canlifins$EntryAgeF)[ordre,]
  
  del=c(del[ordre],1)
  
  xinf=max(Z_ordered[,1],Z_ordered[,2])+1
  Z1=c(Z_ordered[,1],xinf)
  Z2=c(Z_ordered[,2],xinf)
  
  return(list(Z1,Z2,del))
}

Fun_ThetaHatFn12 <- function(Z1,Z2,del,n,phat,Fbar){
  
  # Starting point of Newton-Raphson
  
  Tau_hat = 4*(t(phat) %*% Fbar)-1
  Theta0 = 2*Tau_hat/(1-Tau_hat) # Clayton
  #Theta0 = 1+4(incgam(3,1.2)/)
  
  #s = 10**(-3)
  
  # # MLE of theta
  
  # Fn1
  
  az1=matrix(rep(Z1,n+1),ncol=n+1)
  A=t(az1)>=az1
  
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  Id = diag(1,n+1)
  
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fn1bar=as.vector(MMinv%*%b) # Fn1 bar
  Fn1 = 1-Fn1bar # Fn1
  #Fn1[Fn1 < s]=s # too small values will yield infinity when inverted
  
  # Fn2
  
  Z2_ordered = Z2[order(Z2)]
  az2=matrix(rep(Z2_ordered,n+1),ncol=n+1)
  A=t(az2)>=az2
  
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  del2 = del[order(Z2)]
  b=c((1-eps)*del2[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fn2bar=as.vector(MMinv%*%b)
  Fn2 = 1-Fn2bar
  #Fn2[Fn2 < s]=s
  
  # MLE estimate of theta
  # theta_hat[i] = optimise(f=LogL,interval=c(0,10**2),maximum = TRUE)$maximum
  
  Fn1=Fn1[-(n+1)]
  Fn2=Fn2[-(n+1)]
  
  # Log-pseudo-likelihood function derivative aka score function
  DiffLogL = function(x){ # Clayton
    #x = 10
    sum(phat*(rep(1/(x+1),times=n)-
                log(Fn1*Fn2)+
                (rep(1/x**2,times=n))*log(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1)+
                (rep(2+1/x,times=n))*(Fn1**(rep(-x,times=n))*log(Fn1)+ 
                # this last term (the line above + the 2 below) has limit log(a)(log(b)+2)+2log(b)
                                             Fn2**(rep(-x,times=n))*log(Fn2))/
                                      (Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1)),na.rm=TRUE)
  
    #term1 = suppressWarnings(log(Fn1*Fn2))
    #term1[is.na(term1)] = -10**6# min(term1[!is.na(term1)])-100
    #term2 = suppressWarnings((rep(1/x**2,times=n))*log(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1)+
    #  (rep(2+1/x,times=n))*(Fn1**(rep(-x,times=n))*log(Fn1)+
    #                                Fn2**(rep(-x,times=n))*log(Fn2))/
    #                      (Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1))
    # term2+term3 has limit log(a)(log(b)+2)+2log(b)
    #term2[is.na(term2)] = -10**6# min(term3[!is.na(term3)])-100
    
    #return(sum(phat*(rep(1/(x+1),times=n)-term1+term2)))
    #return(na.omit(result))
   }
    
  #X=matrix(seq(from=0,to=15,by=0.01)); ForPlot=apply(X=X,MARGIN=1,FUN=DiffLogL);
  #plot(X,ForPlot,type='l',ylim=range(c(ForPlot,0)));abline(h=0);
  
  #  DiffLogL = function(x){ # 4.2.20 Nelsen
  #    (-1/x)*sum(phat*((-exp(Fn1**(-rep(-x,times=n)))*Fn1**(-rep(-x,times=n))*log(Fn1)-exp(Fn2**(-rep(-x,times=n)))*Fn2**(-rep(-x,times=n))*log(Fn2))/
  #      ((exp(Fn1**(-rep(-x,times=n)))+exp(Fn2**(-rep(-x,times=n)))-exp(1))*log(exp(Fn1**(-rep(-x,times=n)))+exp(Fn2**(-rep(-x,times=n)))-exp(1)))))
  #  }
    
  #theta_hat = newtonRaphson(fun=DiffLogL, x0=Theta0)$root
  theta_hat = uniroot(f=DiffLogL, interval=c(0.01,15))$root
  
  # Coding Newton-Raphson:
  # https://rpubs.com/aaronsc32/newton-raphson-method#:~:text=%23%23%20%5B1%5D%203.162278-,Newton%2DRaphson%20Method%20in%20R,rootSolve%20package%20features%20the%20uniroot.
  
  return(list(theta_hat,Fn1,Fn2))
} 

Fun_VarMLE <- function(Z1,Z2,del,theta_hat,Fn1,Fn2,n){
  
  # # Variance of An
  
  # Phi
  
  Phi = rep(1/(theta_hat+1),times=n)-
    log(Fn1*Fn2)+
    rep(1/theta_hat**2,times=n)*log(Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)-
    rep(2+1/theta_hat,times=n)*(-Fn1**rep(-theta_hat,times=n)*log(Fn1)-Fn2**rep(-theta_hat,times=n)*log(Fn2))/(Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)
  
  # VStar
  
  az1=matrix(rep(Z1,n+1),ncol=n+1)
  A=(t(az1)>=az1)
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  Id = diag(1,n+1)
  
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fbar=as.vector(MMinv%*%b)
  
  Fbar=Fbar[-(n+1)]
  b=b[-(n+1)]
  A=A[-(n+1),][,-(n+1)]
  B=B[-(n+1),][,-(n+1)]
  Id=diag(1,n)
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  
  D=(1-A)*(1-t(A))*(A%*%t(A))
  bf=b*Fbar
  BF=diag(bf)
  S=rbind(A%*%BF,t(bf))
  R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
  U=(t(M)%*%R)%*%M
  V_z0=(MMinv%*%U)%*%MMinv
  
  Z2_ordered = Z2[order(Z2)]
  az2=matrix(rep(Z2_ordered,n+1),ncol=n+1)
  A=(t(az2)>=az2)
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  del2 = del[order(Z2)]
  
  b=c((1-eps)*del2[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  Id=diag(1,n+1)
  
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fbar=as.vector(MMinv%*%b)
  
  Fbar=Fbar[-(n+1)]
  b=b[-(n+1)]
  A=A[-(n+1),][,-(n+1)]
  B=B[-(n+1),][,-(n+1)]
  Id=diag(1,n)
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  
  D=(1-A)*(1-t(A))*(A%*%t(A))
  bf=b*Fbar
  BF=diag(bf)
  S=rbind(A%*%BF,t(bf))
  R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
  U=(t(M)%*%R)%*%M
  V_0z=(MMinv%*%U)%*%MMinv
  
  Fun_PhiStar <- function(x,y,theta){
    -1/x-
      rep(1/theta,times=n)*x**rep(-theta-1,times=n)/(x**rep(-theta,times=n)+y**rep(-theta,times=n)-1)+
      rep(2+1/theta,times=n)*x**rep(-theta-1,times=n)*((x**rep(-theta,times=n)+y**rep(-theta,times=n)-1)*(-rep(theta,times=n)*log(x)+1)+
                                                         rep(theta,times=n)*(x**rep(-theta,times=n)*log(x)+y**rep(-theta,times=n)*log(y)))/
      (x**rep(-theta,times=n)+y**rep(-theta,times=n)-1)**2
    }
  
  #PhiStar1  = -1/Fn1-
  #  rep(1/theta_hat,times=n)*Fn1**rep(-theta_hat-1,times=n)/(Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)+
  #  rep(2+1/theta_hat,times=n)*Fn1**rep(-theta_hat-1,times=n)*((Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)*(-rep(theta_hat,times=n)*log(Fn1)+1)+
  #                                                               rep(theta_hat,times=n)*(Fn1**rep(-theta_hat,times=n)*log(Fn1)+Fn2**rep(-theta_hat,times=n)*log(Fn2)))/
  #  (Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)**2
  
  PhiStar1 = Fun_PhiStar(x=Fn1,y=Fn2,theta=theta_hat)
  PhiStar2 = Fun_PhiStar(x=Fn2,y=Fn1,theta=theta_hat)
  
  #PhiStar2  = -1/Fn2-
  #  rep(1/theta_hat,times=n)*Fn2**rep(-theta_hat-1,times=n)/(Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)+
  #  rep(2+1/theta_hat,times=n)*Fn2**rep(-theta_hat-1,times=n)*((Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)*(-rep(theta_hat,times=n)*log(Fn2)+1)+
  #                                                               rep(theta_hat,times=n)*(Fn1**rep(-theta_hat,times=n)*log(Fn1)+Fn2**rep(-theta_hat,times=n)*log(Fn2)))/
  #  (Fn1**rep(-theta_hat,times=n)+Fn2**rep(-theta_hat,times=n)-1)**2
  
  #PhiStar1[is.nan(PhiStar1)] <- 0
  #PhiStar2[is.nan(PhiStar2)] <- 0
  
  PhiStar1_matrix = matrix(rep(PhiStar1,times=n),ncol=n)
  PhiStar2_matrix = matrix(rep(PhiStar2,times=n),ncol=n)
  
  VStar = V_z0*PhiStar1_matrix + V_0z*PhiStar2_matrix
  
  # R_hat
  
  az1=matrix(rep(Z1,n+1),ncol=n+1)
  az2=matrix(rep(Z2,n+1),ncol=n+1)
  A=(t(az1)>=az1)*(t(az2)>=az2)
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  Id=diag(rep(1,n+1))
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fbar=as.vector(MMinv%*%b)
  
  Fbar=Fbar[-(n+1)]
  b=b[-(n+1)]
  A=A[-(n+1),][,-(n+1)]
  B=B[-(n+1),][,-(n+1)]
  Id=diag(1,n)
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  
  D=(1-A)*(1-t(A))*(A%*%t(A))
  
  Matrix1 = rbind(Id - A%*%B,-t(b))
  
  #Matrix3 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%(Phi*PhiStar1))
  Matrix3 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%diag(Fbar)%*%B
  #Matrix4 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%(Phi*PhiStar2))
  
  Matrix4 = ginv(Matrix1)%*%Matrix3
  
  r1_hat = Matrix4%*%(Phi*PhiStar1)
  r2_hat = Matrix4%*%(Phi*PhiStar2)
  
  # e_hat
  
  #Matrix5 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%(Phi*Phi))
  
  e_hat = Matrix4%*%(Phi*Phi)
  
  # V (variance of Fbar)
  
  V=(MMinv%*%U)%*%MMinv
  
  # Var(An)
  
  Var_An = t(PhiStar1)%*%diag(Fbar)%*%B%*%V%*%B%*%diag(Fbar)%*%PhiStar2+
    2*b%*%VStar%*%diag(Fbar)%*%b+
    2*t(r1_hat+r2_hat)%*%diag(Fbar)%*%b+
    t(Phi)%*%B%*%V%*%B%*%Phi+
    2*b%*%e_hat+
    t(Phi)%*%B%*%diag(Fbar)%*%(diag(Fbar)%*%b+B%*%D%*%B%*%diag(Fbar)%*%B%*%Phi)
  
  # # Variance of sqrt(n)(theta_hat - theta)
  
  VarTheta = Var_An/(t(Phi)**2%*%Fbar)
  
  return(VarTheta)
  
}



# # # Kendall's tau:simulated data

n = 500
Max = 1000 # number of samples for MSE 
Tau_hat = rep(NA,times=Max)
del_mean = rep(NA,times=Max)
theta = 3
beta=2
lambda=2

for(m in 1:Max){

  # # Estimator 
  
  ZDel = FunZDel(n=n,theta=theta,seed=m,lambda=lambda,beta=beta)
  Z1 = ZDel[[1]]
  Z2 = ZDel[[2]]
  del = ZDel[[3]]
  
  del_mean[m] = mean(del)
  
  az1=matrix(rep(Z1,n+1),ncol=n+1)
  az2=matrix(rep(Z2,n+1),ncol=n+1)
  A=(t(az1)>=az1)*(t(az2)>=az2)
  # Because of ties, A may not be quite upper triangular. We need to convert some numbers to 0.
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  eps=1/(n+1)
  
  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  Id=diag(rep(1,n+1))
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fbar=as.vector(MMinv%*%b)
  
  # The function 'solve' has some trouble sinverting some A's for n >= 1500
  # (example with seed=4). 'ginv' can do it, but is more time consuming. 
  # A tradeoff is to use ginv where solve fails. 
  
  #if("try-error" %in% class(try(solve(A)))){
  #  phat=(ginv(A)%*%Fbar)
  #} else {
   phat=(solve(A)%*%Fbar)#[-(n+1)]
  #}
  
   phat=phat[-(n+1)]
   Fbar=Fbar[-(n+1)]
   
  # # Kendall tau estimate
  
  Tau_hat[m] = 4*(t(phat) %*% Fbar)-1

}


# Graph for bias in estimating Kendall's tau

Tau = theta/(theta+2)
#plot(Tau_hat,type='l',ylim=range(c(Tau_hat,Tau)));abline(h=Tau);

Xlabels = seq(from=0,to=Max,by=100)
Ylabels = seq(from=0,to=1,by=0.2)
Ylim = c(0,1)
par(mfrow=c(1,1))
plot(Tau_hat,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none",ylim=Ylim)
abline(h=Tau,col='red',lwd = 3)
#abline(h=mean(Tau_hat),col='black',lwd = 3)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "iteration rank", font=1.5,cex=1.25)
legend("bottomright",legend=c(expression(tau), expression(tau[500])),pch=c(NA,1),lwd = 2.5,col=c("red", "black"),lty=c(1,NA),cex=1.25)#,bty="n"


# MSE and its decomposition

Var = mean((Tau_hat-mean(Tau_hat))**2)
Bias2 = (mean(Tau_hat)-Tau)**2

MSE=Var+Bias2

# Fbar Variance

b=b[-(n+1)]
A=A[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf)
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# # Kendall's tau variance

Matrix1 = rbind(Id - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = 4**2*(t(Fbar)%*%B%*%V%*%B%*%Fbar+
  2*t(Fbar)%*%B%*%W_hat+
  t(Fbar)%*%B%*%diag(Fbar)%*%(Id+B%*%D%*%B)%*%diag(Fbar)%*%B%*%Fbar)

# # Graphic for Fbar (not good for big n)

FBarFun = function(x,y){
  sum(phat[((Z1 >= x)*(Z2 >= y))])
}

x = unique(Z1[order(Z1)])#
y = unique(Z2[order(Z2)])#
F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))

persp(x,y,F_bar_grid, theta = 30, phi = 30)



# # # Real data

# Epidemiology Data from 
# McGilchrist, C. A., and C. W. Aisbett. "Regression with Frailty in Survival Analysis."
# Biometrics, vol. 47, no. 2, 1991, pp. 461-466.


#location = 'C:/Users/djimd/OneDrive/Documents/Concordia - PhD/Thesis/McGilchrist_Aisbett-1991.pdf'
# Data: https://rdrr.io/cran/SurvCorr/man/kidney.html

# Extract the table
data(kidney)
# extract_tables(file=location,pages=6)

#data_mc = data.frame(mydata)[-(1:3),2:3]

# strsplit(toString(data.frame(mydata)[-(1:3),2]),",")

#times = matrix(as.numeric(unlist(strsplit(toString(data.frame(mydata)[-(1:3),2]),","))),ncol = 2,byrow =TRUE)

#unlist(strsplit(substring(data.frame(mydata)[-(1:3),3],1,4),","))

#censor = substring(data.frame(mydata)[-(1:3),3],1,4)
#censor[7] = sub(" ",",",censor[7],fixed=TRUE)
#censor[34] = sub("?",",",censor[34],fixed=TRUE)

#censoring = matrix(as.numeric(unlist(strsplit(censor,","))),ncol = 2,byrow =TRUE)

#KidneyInfection = cbind.data.frame(times,censoring)
#colnames(KidneyInfection) = c("T1","T2","uncensored1","uncensored2")

# # Estimator 

n = dim(kidney)[1]

ordre = order(kidney$TIME1,kidney$TIME2)#
#Z_ordered = cbind(kidney$T1,kidney$T2)[ordre,]

del1 = kidney$STATUS1 
del2 = kidney$STATUS2

xinf=max(kidney$TIME1,kidney$TIME2)+1 # CREATING POINT AT INFINITY
Z1=c(kidney$TIME1[ordre],xinf)
Z2=c(kidney$TIME2[ordre],xinf)

del1=c(del1,1)
del2=c(del2,1)
del=del1*del2
del=del[ordre]

az1=matrix(rep(Z1,n+1),ncol=n+1)
az2=matrix(rep(Z2,n+1),ncol=n+1)
A=(t(az1)>=az1)*(t(az2)>=az2)
# Because of ties, A may not be quite upper triangular. We need to convert some numbers to 0.
A[lower.tri(A, diag = FALSE)] = 0
rsumA=apply(A,1,sum)

eps=1/(n+1)

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

Id=diag(rep(1,n+1))
M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=as.vector(MMinv%*%b) ### <---- THIS IS THE \bar{F} VECTOR
phat=(solve(A)%*%Fbar)#[-(n+1)] # weights


# Variance estimator

Fbar=Fbar[-(n+1)]
phat=phat[-(n+1)]
b=b[-(n+1)]
A=A[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf)
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# # Graph of Fbar (good)

FBarFun = function(x,y){
  sum(phat[((Z1 >= x)*(Z2 >= y))])
}

x = unique(Z1)
y = unique(Z2[order(Z2)])
F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))

persp(x,y,F_bar_grid, theta = 30, phi = 30)

# # Kendall's tau estimate

Tau_hat = 4*(t(phat) %*% Fbar)-1

# # Kendall's tau variance

Matrix1 = rbind(Id - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = t(Fbar)%*%B%*%V%*%B%*%Fbar+
        2*t(Fbar)%*%B%*%W_hat+
        t(Fbar)%*%B%*%diag(Fbar)%*%(Id+B%*%D%*%B)%*%diag(Fbar)%*%B%*%Fbar




# # # Insurance data from 
# # Frees, Edward W., et al. "Annuity Valuation with Dependent Mortality." 
# The Journal of Risk and Insurance, vol. 63, no. 2, 1996, pp. 229-261.

# Link: http://cas.uqam.ca/

#library(xts)
#library(sp)
#library(zoo)

#install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/R/", type="source")


# Estimator

n=3000

ZDel = FunZDel_canlifins(size=n)
Z1 = ZDel[[1]]
Z2 = ZDel[[2]]
del = ZDel[[3]]

# Graph of the dependence between the two lifetimes

Xlabels = Ylabels = seq(from=0,to=110,by=20)
Xlim = Ylim = range(c(Z1[-(n+1)],Z2[-(n+1)]))
par(mfrow=c(1,1))
plot(Z1[-(n+1)],Z2[-(n+1)],lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none",xlim=Xlim,ylim=Ylim)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "male lifetime", font=1,cex=1.25)
mtext(side=2, line=2.25, "female lifetime", font=1,cex=1.25)

# Phat and Fbar

az1=matrix(rep(Z1,n+1),ncol=n+1)
az2=matrix(rep(Z2,n+1),ncol=n+1)
A=(t(az1)>=az1)*(t(az2)>=az2)
# Because of ties, A may not be quite upper triangular. We need to convert some numbers to 0.
A[lower.tri(A, diag = FALSE)] = 0
rsumA=apply(A,1,sum)

eps=1/(n+1)

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

Id=diag(rep(1,n+1))
M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=as.vector(MMinv%*%b) ### <---- THIS IS THE \bar{F} VECTOR
phat=(ginv(A)%*%Fbar)#[-(n+1)] # weights

phat = phat[-(n+1)]
Fbar=Fbar[-(n+1)]

# Graph of Fbar: (not good for n=200,500, good for n=100,300)

FBarFun = function(x,y){
  sum(phat[((Z1 >= x)*(Z2 >= y))])
}

x = unique(Z1[order(Z1)])
y = unique(Z2[order(Z2)])
F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))

persp(x,y,F_bar_grid, theta = 30, phi = 30)

# #  Pseudo-likelihood maximization

ThetaFn12 = Fun_ThetaHatFn12(Z1,Z2,del,n,phat,Fbar)
theta_hat = ThetaFn12[[1]]

# Variance estimator

#Fbar=Fbar[-(n+1)]
#phat=phat[-(n+1)]
b=b[-(n+1)]
A=A[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf)
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# # Kendall's tau estimate

Tau_hat = 4*(t(phat) %*% Fbar)-1

# # Kendall's tau variance

Matrix1 = rbind(diag(1,n) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(diag(1,n) + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = t(Fbar)%*%B%*%V%*%B%*%Fbar+
  2*t(Fbar)%*%B%*%W_hat+
  t(Fbar)%*%B%*%diag(Fbar)%*%(diag(1,n+1)+B%*%D%*%B)%*%diag(Fbar)%*%B%*%Fbar

# # MLE and variance

#theta_hat = ThetaFn12[[1]]
Fn1 = ThetaFn12[[2]]
Fn2 = ThetaFn12[[3]]

VarThetaHat = Fun_VarMLE(Z1=Z1,Z2=Z2,del=del,theta_hat=as.vector(theta_hat),Fn1=Fn1,Fn2=Fn2,n=n)

 

# #  # Convergence of the stimator (Simulation)

n_vect = seq(from=100,to=500,by=50)

theta_hat_vect = rep(NA,length(n_vect))

theta = 3

beta = 2

for(i in 1:length(n_vect)){
  
  n = n_vect[i]
  ZDel = FunZDel(n=n,theta=theta,beta=beta,seed=i)
  Z1 = ZDel[[1]]
  Z2 = ZDel[[2]]
  del = ZDel[[3]]
  
  az1=matrix(rep(Z1,n+1),ncol=n+1)
  az2=matrix(rep(Z2,n+1),ncol=n+1)
  A=(t(az1)>=az1)*(t(az2)>=az2)
  # Because of ties, A may not be quite upper triangular. We need to convert some numbers to 0.
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  eps=1/(n+1)
  
  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  Id=diag(rep(1,n+1))
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fbar=as.vector(MMinv%*%b) 
  phat=(solve(A)%*%Fbar)
  
  Fbar = Fbar[-(n+1)]
  phat = phat[-(n+1)]
  
  #FBarFun = function(x,y){
  #  sum(phat[((Z1 >= x)*(Z2 >= y))])
  #}
  
  #x = unique(Z1[order(Z1)])#
  #y = unique(Z2[order(Z2)])#
  #F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun)) 
  
  #persp(x,y,F_bar_grid, theta = 30, phi = 30)# works
  
  theta_hat_vect[i] = Fun_ThetaHatFn12(Z1,Z2,del,n,phat,Fbar)[[1]]
  
}

# plot(theta_hat_vect,type='l',ylim=range(c(theta_hat_vect,theta)));abline(h=theta);


# # MLE and variance (through simulated)

theta_vect =  seq(from=0.5,to=10,by=1.5)# exp(seq(from=0.1,to=1,by=0.1))-1 #  
VarThetaHat = rep(NA,times=length(theta_vect))

beta = 2
n=1000

for(j in 1:length(theta_vect)){
  
  ZDel = FunZDel(n=n,theta=theta_vect[j],lambda=1.5,beta=2,seed=j)
  Z1 = ZDel[[1]]
  Z2 = ZDel[[2]]
  del = ZDel[[3]]
  
  az1=matrix(rep(Z1,n+1),ncol=n+1)
  az2=matrix(rep(Z2,n+1),ncol=n+1)
  A=(t(az1)>=az1)*(t(az2)>=az2)
  # Because of ties, A may not be quite upper triangular. We need to convert some numbers to 0.
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  eps=1/(n+1)
  
  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  Id=diag(rep(1,n+1))
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fbar=as.vector(MMinv%*%b) 
  phat=(solve(A)%*%Fbar)
  
  Fbar=Fbar[-(n+1)]
  phat=phat[-(n+1)]
  
  #FBarFun = function(x,y){
  #  sum(phat[((Z1 >= x)*(Z2 >= y))])
  #}
  
  #x = unique(Z1[order(Z1)])#
  #y = unique(Z2[order(Z2)])#
  #F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))
  
  #persp(x,y,F_bar_grid, theta = 30, phi = 30) # works here
  
  ThetaFn12 = Fun_ThetaHatFn12(Z1,Z2,del,n,phat,Fbar)
  
  theta_hat = as.vector(ThetaFn12[[1]])
  Fn1 = ThetaFn12[[2]]
  Fn2 = ThetaFn12[[3]]
  
  VarThetaHat[j] = Fun_VarMLE(Z1=Z1,Z2=Z2,del=del,theta_hat=theta_hat,Fn1=Fn1,Fn2=Fn2,n=n)
}


# # # MLE with Epidemiology Data from 
# McGilchrist, C. A., and C. W. Aisbett. "Regression with Frailty in Survival Analysis."
# Biometrics, vol. 47, no. 2, 1991, pp. 461-466.

# Assuming a copula structure

# Data: https://rdrr.io/cran/SurvCorr/man/kidney.html

# Extract the table
data(kidney)

# Estimator

n=dim(kidney)[1]

del = (kidney$STATUS1 == 1) * (kidney$STATUS2 == 1)

ordre = order(kidney$TIME1,kidney$TIME2)
Z_ordered = cbind(kidney$TIME1,kidney$TIME2)[ordre,]

del=c(del[ordre],1)

xinf=max(Z_ordered[,1],Z_ordered[,2])+1
Z1=c(Z_ordered[,1],xinf)
Z2=c(Z_ordered[,2],xinf)

# Phat and Fbar

az1=matrix(rep(Z1,n+1),ncol=n+1)
az2=matrix(rep(Z2,n+1),ncol=n+1)
A=(t(az1)>=az1)*(t(az2)>=az2)
# Because of ties, A may not be quite upper triangular. We need to convert some numbers to 0.
A[lower.tri(A, diag = FALSE)] = 0
rsumA=apply(A,1,sum)

eps=1/(n+1)

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

Id=diag(rep(1,n+1))
M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=as.vector(MMinv%*%b) ### <---- THIS IS THE \bar{F} VECTOR
phat=(ginv(A)%*%Fbar)#[-(n+1)] # weights

phat = phat[-(n+1)]
Fbar=Fbar[-(n+1)]

# Graph of Fbar

FBarFun = function(x,y){
  sum(phat[((Z1 >= x)*(Z2 >= y))])
}

x = unique(Z1[order(Z1)])
y = unique(Z2[order(Z2)])
F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))

persp(x,y,F_bar_grid, theta = 30, phi = 30)

# #  Pseudo-likelihood maximization

ThetaFn12 = Fun_ThetaHatFn12(Z1,Z2,del,n,phat,Fbar)
theta_hat = ThetaFn12[[1]]

# Variance estimator

#Fbar=Fbar[-(n+1)]
#phat=phat[-(n+1)]
b=b[-(n+1)]
A=A[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf)
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# # Kendall's tau estimate

Tau_hat = 4*(t(phat) %*% Fbar)-1

# # Kendall's tau variance

Matrix1 = rbind(diag(1,n) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(diag(1,n) + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = t(Fbar)%*%B%*%V%*%B%*%Fbar+
  2*t(Fbar)%*%B%*%W_hat+
  t(Fbar)%*%B%*%diag(Fbar)%*%(diag(1,n)+B%*%D%*%B)%*%diag(Fbar)%*%B%*%Fbar

# # MLE and variance

#theta_hat = ThetaFn12[[1]]
Fn1 = ThetaFn12[[2]]
Fn2 = ThetaFn12[[3]]

VarThetaHat = Fun_VarMLE(Z1=Z1,Z2=Z2,del=del,theta_hat=as.vector(theta_hat),Fn1=Fn1,Fn2=Fn2,n=n)

