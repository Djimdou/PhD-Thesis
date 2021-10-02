
#install.packages("copula")

library(copula)
library(MASS)

# # Functions

GammaWangWellsFunc <- function(dataset){
  
  # Call example: GammaWangWellsFunc(kidney)
  
  require(mhazard) # for Dabrowska (1988) estimator for bivariate survivor
  
  # Dabrowska estimator, at grid of the observed data
  
  T1.unique <- unique(c(0,dataset$TIME1))[order(unique(c(0,dataset$TIME1)))]
  T2.unique <- unique(c(0,dataset$TIME2))[order(unique(c(0,dataset$TIME2)))]
  F_dab <- matrix(NA,nrow=length(T1.unique),ncol=length(T2.unique))
  
  model_dab <- KM2(Y1=dataset$TIME1,Y2=dataset$TIME2,Delta1=dataset$STATUS1,Delta2=dataset$STATUS2,estimator="dabrowska")
  
  for(i in 1:length(T1.unique)){
    for(j in 1:length(T2.unique)){
      
      pos.i = 1+min(which(model_dab$T1 >= T1.unique[i]))
      pos.j = 1+min(which(model_dab$T2 >= T1.unique[j]))
      
      if(is.finite(pos.i) & is.finite(pos.j)){
        F_dab[i,j] <- model_dab$Fhat[pos.i,pos.j]
      }else{
        if(is.infinite(pos.i)){pos.i = max(dataset$TIME1)+1}
        if(is.infinite(pos.j)){pos.j = max(dataset$TIME2)+1}
        
        F_dab[i,j] <- KM2(Y1=dataset$TIME1,Y2=dataset$TIME2,Delta1=dataset$STATUS1,Delta2=dataset$STATUS2,newT1=pos.i,newT2=pos.j,estimator="dabrowska")$Fhat_est
        
      }
      #F_dab[i,j] <- model_dab$Fhat[1+min(which(model_dab$T1 >= T1.unique[i])),1+min(which(model_dab$T2 >= T2.unique[j]))]
    }
  }
  
  F_dab[1,1] <- 1 # needed? shouldn't be automatic from the computations?
  
  # Masses, at grid of the observed data
  
  F_dab_mass <- matrix(NA,nrow=dim(F_dab)[1],ncol=dim(F_dab)[2])
  F_dab_mass[1,] <- F_dab[1,]
  F_dab_mass[,1] <- F_dab[,1]
  
  #for(i in 2:dim(F_dab)[1]){
  #  for(j in 2:dim(F_dab)[2]){
  #    F_dab_mass[i,j] <- F_dab[i,j]-F_dab[i-1,j]-F_dab[i,j-1]+F_dab[i-1,j-1]
  #  }
  #}
  
  F_dab_mass[2:dim(F_dab)[1],2:dim(F_dab)[2]] <- F_dab[2:dim(F_dab)[1],2:dim(F_dab)[2]] -
    F_dab[1:(dim(F_dab)[1]-1),2:dim(F_dab)[2]] -
    F_dab[2:dim(F_dab)[1],1:(dim(F_dab)[2]-1)] +
    F_dab[1:(dim(F_dab)[1]-1),1:(dim(F_dab)[2]-1)] 
  
  # Kendall's tau (tau_tilde of Wang and Wells)
  
  tau_tilde <- 4*sum(F_dab_mass[-1,-1] * F_dab[-1,-1])-1
  #tau_tilde
  
  # Kendall's tau modified for ties (gamma of Wang and Wells)
  
  # ties sets
  O1 <- unique(dataset[duplicated(dataset$TIME1) & dataset$STATUS1==1,"TIME1"])
  O2 <- unique(dataset[duplicated(dataset$TIME2) & dataset$STATUS2==1,"TIME2"])
  
  if(length(O1)*length(O2)!=0){# if ties on both T1 and T2. Need to rewrite this later, to take into account case of ties on one variable
    
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

FunZDel_WangWells = function(n,theta,lambda=1,seed,alpha=2,beta=2){
  
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
  
  #ordre = order(Z1,Z2)
  
  #xinf=max(Z1,Z2)+1 # point at infinity
  #Z1=c(Z1[ordre],xinf)
  #Z2=c(Z2[ordre],xinf)
  
  del1 = as.integer(X <= Cx)
  del2 = as.integer(Y <= Cy)
  
  #del1=c(del1,1)
  #del2=c(del2,1)
  #del=del1*del2
  #del=del[ordre]
  
  return(list(Z1,Z2,del1,del2))
}

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


# # Simulations

n = 500 # 500
#Max = 10 # 100 number of samples for MSE 
Theta = seq(from=0.1,to=5,by=0.2)
Tau = Theta/(Theta+2)
V_tau = Tau_hat = rep(NA,times=length(Theta)) # gamma
alpha=beta=2
lambda=1/10

for(t in 1:length(Theta)){
  
  theta <- Theta[t]
  
  MyCopula <- mvdc(copula=claytonCopula(param=theta), # Ali-Mikail-Haq copula for (F(X), F(Y)), theta should be in [0,1]
                   margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                   paramMargins=list(shape=alpha,scale=beta) # alpha:shape, beta:scale
                  )
  set.seed(t)
  X <- rMvdc(n=n,MyCopula)
  X1 <- X[,1]
  X2 <- X[,2]
  
  # Censoring variable
  
  #if(CopulaName=="Clayton"){
  C1 = rexp(n=n,rate = lambda)
  C2 = rexp(n=n,rate = lambda)
  #}
  
  # Observations
  
  Z1 <- pmin(X1,C1)
  Z2 <- pmin(X2,C2)
  
  ordre = order(Z1,Z2)
  
  xinf=max(Z1,Z2)+1 # point at infinity
  Z1=c(Z1[ordre],xinf)
  Z2=c(Z2[ordre],xinf)
  
  del1 = as.integer(X1 <= C1)
  del2 = as.integer(X2 <= C2)
  
  del1=del1[ordre]
  del2=del2[ordre]
  
  del1=c(del1,1)
  del2=c(del2,1)
  del=del1*del2
  #del=del[ordre]
  
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
  
  D=(1-A)*(1-t(A))*(A%*%t(A))
  
  bf=b*Fbar
  BF=diag(bf)
  S=rbind(A%*%BF,t(bf))
  R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
  U=(t(M)%*%R)%*%M
  V=(MMinv%*%U)%*%MMinv # Fbar variance
  
  # No longer need xinf and its mass
  
  phat=phat[-(n+1)]
  Fbar=Fbar[-(n+1)]
  A=A[-(n+1),-(n+1)]
  B=B[-(n+1),-(n+1)]
  D=D[-(n+1),-(n+1)]
  V=V[-(n+1),-(n+1)]
  b=b[-(n+1)]
  
  Id=diag(rep(1,n))
  
  # # Kendall tau estimate
  
  Tau_hat[t] = 4*(t(phat) %*% Fbar)-1
  
  #if((m %% 50)==0){print(m)}
  
  # Variance estimator
  
  Matrix1 = rbind(Id - A%*%B,-b)
  Matrix2 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%Fbar)
  
  W_hat = ginv(Matrix1)%*%Matrix2
  
  V_tau[t] = 4**2*(t(Fbar)%*%B%*%V%*%B%*%Fbar+
                  2*t(Fbar)%*%B%*%W_hat+
                  t(Fbar)%*%B%*%diag(Fbar)%*%(Id+B%*%D%*%B)%*%diag(Fbar)%*%B%*%Fbar)
  
}

# not good

CI_upper = Tau_hat+1.96*sqrt(V_tau/n) # 
CI_lower = Tau_hat-1.96*sqrt(V_tau/n) # 
Ylim = range(Tau,CI_upper,CI_lower)
plot(Theta,Tau,type="l",ylim=Ylim,col='green')
lines(Theta,Tau_hat,col='blue')
lines(Theta,CI_upper,col='red')
lines(Theta,CI_lower,col='red')