KB

# # Packages installation and loading

install.packages(c('copula','pracma','SurvCorr','survival','MASS'))# 'tabulizer','tryCatchLog','mhazard','CASdatasets'

library(copula) # for claytonCopula
library(survival) # univariate Kaplan-Meier estimators, kidney data
library(pracma) # for NewtonRaphson procedure
#library(tryCatchLog) # errors and warnings management
library(MASS) # for ginv
#library(CASdatasets)


# # # MSE (simulation data)

# same values for MLE if (lambda,theta)=(1/2,5),(1,1),(1,2),(2,2),(2,5),(1/2,5) # seems like similarity appears for censoring >= 50%
# different values if (lambda,theta)=(1/8,1),(1/4,2),(1/4,1),(1/3,1),(1/4,5) # seems like difference appears for censoring <= 50%, due just to a couple of values

# Pseudo-likelihood function

L_MassShift <- function(theta){
  prod(((theta+1)*(Fn1*Fn2)**(-theta-1)*(Fn1**(-theta)+Fn2**(-theta)-1)**(-2-1/theta))**phat) # , na.rm=TRUE
}

L_ShihLouis <- function(theta){
  prod((theta+1)*(Fn1bar*Fn2bar)**(-theta-1)*(Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(del1*del2)*
         (Fn1bar**(-theta-1)*(Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(-1/theta-1))**(del1*(1-del2))*
         (Fn2bar**(-theta-1)*(Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(-1/theta-1))**((1-del1)*del2)*
         ((Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(-1/theta))**((1-del1)*(1-del2))) # , na.rm=TRUE
}

# Pseudo-loglikelihood functions

logL_MassShift <- function(theta){ # Clayton copula
  sum(
    phat*(log(theta+1)-(theta+1)*log(Fn1*Fn2)-(2+1/theta)*log(Fn1**(-theta)+Fn2**(-theta)-1))
    ,na.rm=TRUE) # 
}

logL_MassShift_AMH <- function(theta){ # AMH copula,  # theta in -1,1
  sum(
    phat*(
      3*log(1-theta*(1-Fn1)*(1-Fn2))+log(1-2*theta+theta**2+2*theta*(1-theta)*Fn1+
                                           theta*(1-theta)*Fn2+2*theta**2*Fn1*Fn2+theta**2*Fn1**2-theta**2*Fn1**2*Fn2)
    )
    ,na.rm=TRUE) # 
}

logL_MassShift_Nelsen <- function(theta){ # Nelsen 4.2.20 copula
  sum(
    phat*(
      log(1+theta)-
        (1/theta+1)*log(log(exp(Fn1**(-theta))+exp(Fn2**(-theta))-exp(1)))-
        (theta+1)*log(Fn1*Fn2) + (Fn1**(-theta)+Fn2**(-theta))-
        2*log(exp(Fn1**(-theta))+exp(Fn2**(-theta))-exp(1))
    )
    ,na.rm=TRUE) #
}

logL_ShihLouis <- function(theta){
  sum(del1*del2*log(theta+1)-(theta+1)*(del1*log(Fn1bar)+del2*log(Fn2bar))-
        (1/theta+del1+del2)*log(Fn1bar**(-theta)+Fn2bar**(-theta)-1)
      ,na.rm=TRUE)
}

logL_ShihLouis_AMH <- function(theta){ # theta in -1,1
  sum(
    del1*del2*log(1-2*theta+theta**2+2*theta*(1-theta)*Fn1+
          theta*(1-theta)*Fn2+2*theta**2*Fn1*Fn2+theta**2*Fn1**2-theta**2*Fn1**2*Fn2)
      ,na.rm=TRUE)
}

# Score functions

Score_MassShift <- function(x){ # Clayton
  sum(phat*(rep(1/(x+1),times=n)-
              log(Fn1*Fn2)+
              (rep(1/x**2,times=n))*log(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1)+
              (rep(2+1/x,times=n))*(Fn1**(rep(-x,times=n))*log(Fn1)+ 
                                      Fn2**(rep(-x,times=n))*log(Fn2))/
              (Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1))) # ,na.rm=TRUE
}

Score_ShihLouis <- function(x){ # Clayton
  sum(del1*del2*(rep(1/(x+1),times=n)-log(Fn1*Fn2))-
        (del1*(1-del2)*log(Fn1)+(1-del1)*del2*log(Fn2))-
        (
          -rep(1/x**2,times=n)*(1-del1*del2)*log(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1)+
            (rep(1/x,times=n)*(1-del1*del2)+del1+del2-2*del1*del2)*((-Fn1**(rep(-x,times=n))*log(Fn1)-Fn2**(rep(-x,times=n))*log(Fn2))/(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1))
        )
  ) # ,na.rm=TRUE
}

Max = 1 # seq(from=100,to=200,by=50)
n = 50
lambda = 0.05 # 1/0.5-1

mean_del = theta_hat_MassShift_vect = theta_hat_ShihLouis_vect = rep(NA,Max)

theta = -1/2

alpha = beta = 2

for(m in 1:Max){
  
  #n = n_vect[i]
  
  
  #FunZDel = function(n,theta,lambda=lambda,alpha=2,beta=2){#,seed
  
  #  require(copula)
  #  require(survival)
  
  MyCopula <- mvdc(copula=claytonCopula(param=theta), # Clayton copula for (F(X), F(Y))
                   margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                   paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale
  
  # MyCopula <- mvdc(copula=amhCopula(param=theta), # Ali-Mikail-Haq copula for (F(X), F(Y)), theta should be in [0,1]
  #                  margins=c("exp","exp"), # exponential distribution for margins X and Y
  #                  paramMargins=list(list(rate=1),list(rate=1)))
  
  set.seed(m)
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
  
  del1=del1[ordre]
  del2=del2[ordre]
  
  del1=c(del1,1)
  del2=c(del2,1)
  del=del1*del2
  #del=del[ordre]
  
  #  return(list(Z1,Z2,del1,del2))
  #}
  
  #ZDel = FunZDel(n=n,theta=theta,lambda=lambda,beta=2)#,seed=m
  #Z1 = ZDel[[1]]
  #Z2 = ZDel[[2]]
  #del1 = ZDel[[3]]
  #del2 = ZDel[[4]]
  #del = del1*del2
  
  
  # Phat and Fbar (for Mass-shifting estimator)
  
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
  
  # # Graph of Fbar: (not good for n=200,500, good for n=100,300)
  
  # FBarFun = function(x,y){
  #   sum(phat[((Z1 >= x)*(Z2 >= y))])
  # }
  
  # x = unique(Z1[order(Z1)])
  # y = unique(Z2[order(Z2)])
  # F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))
  
  # persp(x,y,F_bar_grid, theta = 30, phi = 30)
  
  #FBarFun = function(x,y){
  #  sum(phat[((Z1 >= x)*(Z2 >= y))])
  #}
  
  #x = unique(Z1[order(Z1)])#
  #y = unique(Z2[order(Z2)])#
  #F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun)) 
  
  #persp(x,y,F_bar_grid, theta = 30, phi = 30)# works
  
  # Starting point of Newton-Raphson
  
  # Tau_hat = 4*(t(phat) %*% Fbar)-1 # plug-in estimator for Kendall's tau
  # Theta0 = 2*Tau_hat/(1-Tau_hat) # Clayton
  
  # Univariate Kaplan-Meier estimators
  
  Z1 <- Z1[-(n+1)]
  Z2 <- Z2[-(n+1)]
  del1 <- del1[-(n+1)]
  del2 <- del2[-(n+1)]
  
  model1 <- survfit(Surv(Z1, del1) ~ 1)
  T1 <- model1$time
  Fn1bar <- model1$surv
  Fn1bar <- approx(x=T1, y=Fn1bar, xout=Z1, method="constant", ties = mean)$y # predicting the survival rates at the data times
  Fn1 <- 1-Fn1bar
  
  model2 <- survfit(Surv(Z2, del2) ~ 1)#$surv
  T2 <- model2$time
  Fn2bar <- model2$surv
  Fn2bar <- approx(x=T2, y=Fn2bar, xout=Z2, method="constant", ties = mean)$y  # predicting the survival rates at the data times
  Fn2 <- 1-Fn2bar
  # 
  # Fn1[Fn1==0] <- min(Fn1[Fn1!=0])/10**2
  # Fn2[Fn2==0] <- min(Fn2[Fn2!=0])/10**2
  # 
  # Fn1bar[Fn1bar==0] <- min(Fn1bar[Fn1bar!=0])/10**2
  # Fn2bar[Fn2bar==0] <- min(Fn2bar[Fn2bar!=0])/10**2
  
  # MLE
  
  #theta_hat_MassShift_vect[m] <- optimise(f=L_MassShift,interval=c(0,50),maximum = TRUE)$maximum
  theta_hat_MassShift_vect[m] <- optimise(f=logL_MassShift,interval=c(0,50),maximum = TRUE)$maximum
  #theta_hat_MassShift_vect[m] <- newtonRaphson(fun=Score_MassShift, x0=Theta0)$root
  #theta_hat_MassShift_vect[m] <- uniroot(Score_MassShift, interval=c(0.01,50))$root
  #theta_hat_MassShift_vect[m] <- nlm(f=L_MassShift, p=Theta0)$estimate
  #theta_hat_ShihLouis_vect[m] <- optimise(f=L_ShihLouis,interval=c(0,50),maximum = TRUE)$maximum
  theta_hat_ShihLouis_vect[m] <- optimise(f=logL_ShihLouis,interval=c(0,50),maximum = TRUE)$maximum
  
  mean_del[m] = mean(del)
  
}

100*(1-mean(mean_del))

# c(mean(mean_del)*100,mean(theta_hat_MassShift_vect),mean(theta_hat_ShihLouis_vect),theta,mean(VarTheta))


# # Plots for comparing the two log-likelihoods

Theta <- matrix(seq(from=0.01,to=50,by=1),ncol=1)
Y_MassShift <- apply(X=Theta,MARGIN=1,FUN=logL_MassShift) # logL_MassShift
Y_ShihLouis <- apply(X=Theta,MARGIN=1,FUN=logL_ShihLouis)
#c(length(X),length(Y))
Ylim <- range(c(Y_MassShift,Y_ShihLouis))
plot(Theta,Y_MassShift,type='l',ylim=Ylim)
lines(Theta,Y_ShihLouis,col="red",ylim=Ylim)
#abline(h=0,col="red")
#plot(Theta,Y_ShihLouis,type='l')

# MSE and its decomposition

Var_MassShift = mean((theta_hat_MassShift_vect-mean(theta_hat_MassShift_vect))**2)
Bias2_MassShift = (mean(theta_hat_MassShift_vect)-theta)**2

MSE_MassShift=Var_MassShift+Bias2_MassShift

Var_ShihLouis = mean((theta_hat_ShihLouis_vect-mean(theta_hat_ShihLouis_vect))**2)
Bias2_ShihLouis = (mean(theta_hat_ShihLouis_vect)-theta)**2

MSE_ShihLouis=Var_ShihLouis+Bias2_ShihLouis

# # Plots for comparing the two MSEs

Ylim <- range(c(theta_hat_MassShift_vect,theta_hat_ShihLouis_vect))
plot(theta_hat_MassShift_vect,col="black",ylim=Ylim)
points(theta_hat_ShihLouis_vect,col="blue")
abline(h=theta,col="red")

# # # Variance

# # Variance of An

# Fn1[Fn1==0] <- min(Fn1[Fn1!=0])/10 # values of 0 cause numerical issues
# Fn2[Fn2==0] <- min(Fn2[Fn2!=0])/10 # values of 0 cause numerical issues

theta_hat <- log(theta_hat_MassShift_vect[m])
#theta_hat <- theta_hat_MassShift_vect[m]

# Phi

Phi <- rep(1/(theta_hat+1),times=n)-
  log(Fn1[-(n+1)]*Fn2[-(n+1)])+
  rep(1/theta_hat**2,times=n)*log(Fn1[-(n+1)]**rep(-theta_hat,times=n)+Fn2[-(n+1)]**rep(-theta_hat,times=n)-1)+
  rep(2+1/theta_hat,times=n)*(Fn1[-(n+1)]**rep(-theta_hat,times=n)*log(Fn1[-(n+1)])+Fn2[-(n+1)]**rep(-theta_hat,times=n)*log(Fn2[-(n+1)]))/(Fn1[-(n+1)]**rep(-theta_hat,times=n)+Fn2[-(n+1)]**rep(-theta_hat,times=n)-1)

Phi[is.nan(Phi)] <- Phi[!is.nan(Phi)][1]

# VStar

A_z0 <- A_0z <- matrix(rep(NA,times=(n+1)*(n+1)),ncol=n+1)

A_z0[,n+1] <- 1
A_z0[n+1,] <- 1

A_0z[,n+1] <- 1
A_0z[n+1,] <- 1

for(i in 1:n){
  for(j in 1:n){
    A_z0[i,j] <- Z1[j]>=Z1[i]
    A_0z[i,j] <- Z2[j]>=Z2[i]
    #A[i,j] <- (Z1[j]>=Z1[i])*(Z2[j]>=Z2[i])
  }
}

# b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
# B=diag(b)

Id = diag(1,n+1)

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

# Id = diag(1,n+1)

M=rbind(Id-A_z0%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=as.vector(MMinv%*%b)

Fbar=Fbar[-(n+1)]
b=b[-(n+1)]
A_z0=A_z0[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A_z0%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A_z0)*(1-t(A_z0))*(A_z0%*%t(A_z0))
bf=b*Fbar
BF=diag(bf)
S=rbind(A_z0%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V_z0=(MMinv%*%U)%*%MMinv

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

Id = diag(1,n+1)

M=rbind(Id-A_0z%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=as.vector(MMinv%*%b)

Fbar=Fbar[-(n+1)]
b=b[-(n+1)]
A_0z=A_0z[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A_0z%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A_0z)*(1-t(A_0z))*(A_0z%*%t(A_0z))
bf=b*Fbar
BF=diag(bf)
S=rbind(A_0z%*%BF,t(bf))
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

PhiStar1 = Fun_PhiStar(x=Fn1[-(n+1)],y=Fn2[-(n+1)],theta=theta_hat)
PhiStar2 = Fun_PhiStar(x=Fn2[-(n+1)],y=Fn1[-(n+1)],theta=theta_hat)

PhiStar1[is.nan(PhiStar1)] <- PhiStar1[!is.nan(PhiStar1)][1]
PhiStar2[is.nan(PhiStar2)] <- PhiStar2[!is.nan(PhiStar2)][1]

PhiStar1_matrix = matrix(rep(PhiStar1,times=n),ncol=n)
PhiStar2_matrix = matrix(rep(PhiStar2,times=n),ncol=n)

VStar = V_z0*PhiStar1_matrix + V_0z*PhiStar2_matrix

# R_hat

Z1=c(Z1,xinf)
Z2=c(Z2,xinf)

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
Matrix2 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%diag(Fbar)%*%B
Matrix3 = ginv(Matrix1)%*%Matrix2

r1_hat = Matrix3%*%(Phi*PhiStar1)
r2_hat = Matrix3%*%(Phi*PhiStar2)

# e_hat

#Matrix5 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%(Phi*Phi))

e_hat = Matrix3%*%(Phi*Phi)

# V (variance of Fbar)

V=(MMinv%*%U)%*%MMinv

# Var(An)

Var_An = t(PhiStar1)%*%diag(Fbar)%*%B%*%V%*%B%*%diag(Fbar)%*%PhiStar2+
  2*b%*%VStar%*%diag(Fbar)%*%b+
  2*t(r1_hat+r2_hat)%*%diag(Fbar)%*%b+
  t(Phi)%*%B%*%V%*%B%*%Phi+
  2*b%*%e_hat+
  t(Phi)%*%B%*%diag(Fbar)%*%(diag(Fbar)%*%b+B%*%D%*%B%*%diag(Fbar)%*%B%*%Phi)

# t(PhiStar1)*diag(Fbar)

# # Variance of sqrt(n)(theta_hat - theta)

#VarTheta = Var_An/(t(Phi)**2%*%Fbar)
VarTheta <- Var_An/(t(Phi)**2%*%phat)


# # CAS Canada Life Insurance dataset

n=100

canlifins <- get(load('C:/Users/mloudegu/Downloads/CASdatasets_1.0-11/CASdatasets/data/canlifins.rda'))

#head(canlifins)

# 14,889 contracts where one annuitant is male and the other female

canlifins = canlifins[(canlifins$EntryAgeM>=18) & (canlifins$EntryAgeF>=18),]

seed <- 1
set.seed(seed)

Sample <- sample(1:dim(canlifins)[1],size = n)
canlifins <- canlifins[Sample,]

del1 <- canlifins$DeathTimeM > 0
del2 <- canlifins$DeathTimeF > 0

ordre <- order(canlifins$DeathTimeM+canlifins$EntryAgeM,canlifins$DeathTimeF+canlifins$EntryAgeF)
Z_ordered <- cbind(canlifins$DeathTimeM+canlifins$EntryAgeM,canlifins$DeathTimeF+canlifins$EntryAgeF)[ordre,]

del1 <- c(del1[ordre],1)
del2 <- c(del2[ordre],1)

xinf <- max(Z_ordered[,1],Z_ordered[,2])+1
Z1 <- c(Z_ordered[,1],xinf)
Z2 <- c(Z_ordered[,2],xinf)

del <- del1*del2

# Phat and Fbar (for Mass-shifting estimator)

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

# # Graph of Fbar: (not good for n=200,500, good for n=100,300)

# FBarFun = function(x,y){
#   sum(phat[((Z1 >= x)*(Z2 >= y))])
# }

# x = unique(Z1[order(Z1)])
# y = unique(Z2[order(Z2)])
# F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))

# persp(x,y,F_bar_grid, theta = 30, phi = 30)

#FBarFun = function(x,y){
#  sum(phat[((Z1 >= x)*(Z2 >= y))])
#}

#x = unique(Z1[order(Z1)])#
#y = unique(Z2[order(Z2)])#
#F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun)) 

#persp(x,y,F_bar_grid, theta = 30, phi = 30)# works

# Starting point of Newton-Raphson

# Tau_hat = 4*(t(phat) %*% Fbar)-1 # plug-in estimator for Kendall's tau
# Theta0 = 2*Tau_hat/(1-Tau_hat) # Clayton

# Univariate Kaplan-Meier estimators

Z1 <- Z1[-(n+1)]
Z2 <- Z2[-(n+1)]
del1 <- del1[-(n+1)]
del2 <- del2[-(n+1)]

model1 <- survfit(Surv(Z1, del1) ~ 1)
T1 <- model1$time
Fn1bar <- model1$surv
Fn1bar <- approx(x=T1, y=Fn1bar, xout=Z1, method="constant", ties = mean)$y # predicting the survival rates at the data times
Fn1 <- 1-Fn1bar

model2 <- survfit(Surv(Z2, del2) ~ 1)#$surv
T2 <- model2$time
Fn2bar <- model2$surv
Fn2bar <- approx(x=T2, y=Fn2bar, xout=Z2, method="constant", ties = mean)$y  # predicting the survival rates at the data times
Fn2 <- 1-Fn2bar
# 
# Fn1[Fn1==0] <- min(Fn1[Fn1!=0])/10**2
# Fn2[Fn2==0] <- min(Fn2[Fn2!=0])/10**2
# 
# Fn1bar[Fn1bar==0] <- min(Fn1bar[Fn1bar!=0])/10**2
# Fn2bar[Fn2bar==0] <- min(Fn2bar[Fn2bar!=0])/10**2

# MLE

#theta_hat_MassShift_vect[m] <- optimise(f=L_MassShift,interval=c(0,50),maximum = TRUE)$maximum
theta_hat_MassShift <- optimise(f=logL_MassShift,interval=c(0,50),maximum = TRUE)$maximum
#theta_hat_MassShift <- optimise(f=logL_MassShift_Nelsen,interval=c(0,1.5),maximum = TRUE)$maximum
#theta_hat_MassShift_vect[m] <- newtonRaphson(fun=Score_MassShift, x0=Theta0)$root
#theta_hat_MassShift_vect[m] <- uniroot(Score_MassShift, interval=c(0.01,50))$root
#theta_hat_MassShift_vect[m] <- nlm(f=L_MassShift, p=Theta0)$estimate
#theta_hat_ShihLouis_vect[m] <- optimise(f=L_ShihLouis,interval=c(0,50),maximum = TRUE)$maximum
theta_hat_ShihLouis <- optimise(f=logL_ShihLouis,interval=c(0,50),maximum = TRUE)$maximum

# # Plot for the two log-likelihood

Theta <- matrix(seq(from=0.01,to=2,by=0.01),ncol=1)
Y_MassShift <- apply(X=Theta,MARGIN=1,FUN=logL_MassShift_Nelsen) # logL_MassShift
Y_ShihLouis <- apply(X=Theta,MARGIN=1,FUN=logL_ShihLouis)
#c(length(X),length(Y))
Ylim <- range(c(Y_MassShift,Y_ShihLouis))
plot(Theta,Y_MassShift,type='l')
#lines(Theta,Y_ShihLouis,col="red",ylim=Ylim)
#abline(h=0,col="red")
#plot(Theta,Y_ShihLouis,type='l')

# # # Variance

# # Variance of An

# Fn1[Fn1==0] <- min(Fn1[Fn1!=0])/10 # values of 0 cause numerical issues
# Fn2[Fn2==0] <- min(Fn2[Fn2!=0])/10 # values of 0 cause numerical issues

#theta_hat <- log(theta_hat_MassShift)
theta_hat <- theta_hat_MassShift

# Phi

Phi <- rep(1/(theta_hat+1),times=n)-
  log(Fn1[-(n+1)]*Fn2[-(n+1)])+
  rep(1/theta_hat**2,times=n)*log(Fn1[-(n+1)]**rep(-theta_hat,times=n)+Fn2[-(n+1)]**rep(-theta_hat,times=n)-1)+
  rep(2+1/theta_hat,times=n)*(Fn1[-(n+1)]**rep(-theta_hat,times=n)*log(Fn1[-(n+1)])+Fn2[-(n+1)]**rep(-theta_hat,times=n)*log(Fn2[-(n+1)]))/(Fn1[-(n+1)]**rep(-theta_hat,times=n)+Fn2[-(n+1)]**rep(-theta_hat,times=n)-1)

Phi[is.nan(Phi)] <- Phi[!is.nan(Phi)][1]

# VStar

A_z0 <- A_0z <- matrix(rep(NA,times=(n+1)*(n+1)),ncol=n+1)

A_z0[,n+1] <- 1
A_z0[n+1,] <- 1

A_0z[,n+1] <- 1
A_0z[n+1,] <- 1

for(i in 1:n){
  for(j in 1:n){
    A_z0[i,j] <- Z1[j]>=Z1[i]
    A_0z[i,j] <- Z2[j]>=Z2[i]
    #A[i,j] <- (Z1[j]>=Z1[i])*(Z2[j]>=Z2[i])
  }
}

# b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
# B=diag(b)

Id = diag(1,n+1)

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

# Id = diag(1,n+1)

M=rbind(Id-A_z0%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=as.vector(MMinv%*%b)

Fbar=Fbar[-(n+1)]
b=b[-(n+1)]
A_z0=A_z0[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A_z0%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A_z0)*(1-t(A_z0))*(A_z0%*%t(A_z0))
bf=b*Fbar
BF=diag(bf)
S=rbind(A_z0%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V_z0=(MMinv%*%U)%*%MMinv

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

Id = diag(1,n+1)

M=rbind(Id-A_0z%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=as.vector(MMinv%*%b)

Fbar=Fbar[-(n+1)]
b=b[-(n+1)]
A_0z=A_0z[-(n+1),][,-(n+1)]
B=B[-(n+1),][,-(n+1)]
Id=diag(1,n)
M=rbind(Id-A_0z%*%B,-t(b))
MMinv=solve(t(M)%*%M)

D=(1-A_0z)*(1-t(A_0z))*(A_0z%*%t(A_0z))
bf=b*Fbar
BF=diag(bf)
S=rbind(A_0z%*%BF,t(bf))
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

PhiStar1 = Fun_PhiStar(x=Fn1[-(n+1)],y=Fn2[-(n+1)],theta=theta_hat)
PhiStar2 = Fun_PhiStar(x=Fn2[-(n+1)],y=Fn1[-(n+1)],theta=theta_hat)

PhiStar1[is.nan(PhiStar1)] <- PhiStar1[!is.nan(PhiStar1)][1]
PhiStar2[is.nan(PhiStar2)] <- PhiStar2[!is.nan(PhiStar2)][1]

PhiStar1_matrix = matrix(rep(PhiStar1,times=n),ncol=n)
PhiStar2_matrix = matrix(rep(PhiStar2,times=n),ncol=n)

VStar = V_z0*PhiStar1_matrix + V_0z*PhiStar2_matrix

# R_hat

Z1=c(Z1,xinf)
Z2=c(Z2,xinf)

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
Matrix2 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%diag(Fbar)%*%B
Matrix3 = ginv(Matrix1)%*%Matrix2

r1_hat = Matrix3%*%(Phi*PhiStar1)
r2_hat = Matrix3%*%(Phi*PhiStar2)

# e_hat

#Matrix5 = (rbind(A%*%B%*%diag(Fbar),b%*%diag(Fbar)))%*%(Id + B%*%D%*%B)%*%(diag(Fbar)%*%B%*%(Phi*Phi))

e_hat = Matrix3%*%(Phi*Phi)

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

#VarTheta = Var_An/(t(Phi)**2%*%Fbar)
VarTheta[m] <- Var_An/(t(Phi)**2%*%phat)