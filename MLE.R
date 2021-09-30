# # Packages installation and loading

install.packages(c('copula','pracma','survival','MASS'))# ,'SurvCorr','tabulizer','tryCatchLog','mhazard','CASdatasets'

library(copula) # for claytonCopula
library(survival) # univariate Kaplan-Meier estimators, kidney data
library(pracma) # for NewtonRaphson procedure
#library(tryCatchLog) # errors and warnings management
library(MASS) # for ginv
#library(CASdatasets)


# # # MSE (simulation data)

# same values for MLE if (lambda,theta)=(1/2,5),(1,1),(1,2),(2,2),(2,5),(1/2,5) # seems like similarity appears for censoring >= 50%
# different values if (lambda,theta)=(1/8,1),(1/4,2),(1/4,1),(1/3,1),(1/4,5) # seems like difference appears for censoring <= 50%, due just to a couple of values

# # Pseudo-likelihood function
# 
# L_MassShift <- function(theta){
#   prod(((theta+1)*(Fn1*Fn2)**(-theta-1)*(Fn1**(-theta)+Fn2**(-theta)-1)**(-2-1/theta))**phat) # , na.rm=TRUE
# }
# 
# L_ShihLouis <- function(theta){
#   prod((theta+1)*(Fn1bar*Fn2bar)**(-theta-1)*(Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(del1*del2)*
#          (Fn1bar**(-theta-1)*(Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(-1/theta-1))**(del1*(1-del2))*
#          (Fn2bar**(-theta-1)*(Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(-1/theta-1))**((1-del1)*del2)*
#          ((Fn1bar**(-theta)+Fn2bar**(-theta)-1)**(-1/theta))**((1-del1)*(1-del2))) # , na.rm=TRUE
# }

# Pseudo-loglikelihood functions

logL_MassShift <- function(theta,CopulaName){
  
  if(CopulaName == "Clayton"){ # theta > 0
    logl <- log(theta+1)-(theta+1)*log(Fn1*Fn2)-(2+1/theta)*log(Fn1**(-theta)+Fn2**(-theta)-1)
  }
  
  if(CopulaName == "AMH"){ # AMH copula,  # theta in -1,1
    logl <- -3*log(1-theta*(1-Fn1)*(1-Fn2))+
      log(1-2*theta+theta**2+theta*(1-theta)*(Fn1+Fn2)+theta*(theta+1)*Fn1*Fn2)
  }
  
  if(CopulaName == "Nelsen"){ # Nelsen 4.2.20 copula
    logl <- log(1+theta)-
      (1/theta+1)*log(log(exp(Fn1**(-theta))+exp(Fn2**(-theta))-exp(1)))-
      (theta+1)*log(Fn1*Fn2) + (Fn1**(-theta)+Fn2**(-theta))-
      2*log(exp(Fn1**(-theta))+exp(Fn2**(-theta))-exp(1))
  }
  
  logL <- sum(phat*logl,na.rm=TRUE)
  return(logL)
}


logL_ShihLouis <- function(theta,CopulaName){
  
  if(CopulaName == "Clayton"){ # theta > 0
    logl <- del1*del2*log(theta+1)-(theta+1)*(del1*log(Fn1bar)+del2*log(Fn2bar))-
      (1/theta+del1+del2)*log(Fn1bar**(-theta)+Fn2bar**(-theta)-1)
  }
  
  if(CopulaName == "AMH"){ # -1 <= theta <= 1
    logl <- del1*del2*log(1-2*theta+theta**2+theta*(1-theta)*(Fn1bar+Fn2bar)+theta*(theta+1)*Fn1bar*Fn2bar)-
      (1+del1+del2)*log(1-theta*(1-Fn1bar)*(1-Fn2bar))+
      (1-del1)*log(Fn1bar) + (1-del2)*log(Fn2bar) +
      (1-del1)*del2*log(1-theta*(1-Fn1bar)) + del1*(1-del2)*log(1-theta*(1-Fn2bar))
  }
  
  logL <- sum(logl,na.rm=TRUE)
}

# # Score functions
# 
# Score_MassShift <- function(x){ # Clayton
#   sum(phat*(rep(1/(x+1),times=n)-
#               log(Fn1*Fn2)+
#               (rep(1/x**2,times=n))*log(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1)+
#               (rep(2+1/x,times=n))*(Fn1**(rep(-x,times=n))*log(Fn1)+ 
#                                       Fn2**(rep(-x,times=n))*log(Fn2))/
#               (Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1))) # ,na.rm=TRUE
# }
# 
# Score_ShihLouis <- function(x){ # Clayton
#   sum(del1*del2*(rep(1/(x+1),times=n)-log(Fn1*Fn2))-
#         (del1*(1-del2)*log(Fn1)+(1-del1)*del2*log(Fn2))-
#         (
#           -rep(1/x**2,times=n)*(1-del1*del2)*log(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1)+
#             (rep(1/x,times=n)*(1-del1*del2)+del1+del2-2*del1*del2)*((-Fn1**(rep(-x,times=n))*log(Fn1)-Fn2**(rep(-x,times=n))*log(Fn2))/(Fn1**(rep(-x,times=n))+Fn2**(rep(-x,times=n))-1))
#         )
#   ) # ,na.rm=TRUE
# }

Max = 1 # 500
n = 100 # 500
lambda = 0.05 # 0.05
a = 30 # 3,30 # superior limit of the censoring r.v. uniformly distributed

mean_del = theta_hat_MassShift_vect = theta_hat_ShihLouis_vect = rep(NA,Max)

theta = 1 # copula parameter

CopulaName = "Clayton" # in c("Clayton", "AMH")

if(CopulaName == "Clayton"){
  alpha = beta = 2 # alpha:shape, beta:scale of the margins
  MyCopula <- mvdc(copula=claytonCopula(param=theta), # Ali-Mikail-Haq copula for (F(X), F(Y)), theta should be in [0,1]
                   margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                   paramMargins=list(shape=alpha,scale=beta) # alpha:shape, beta:scale
  )
}

if(CopulaName == "AMH"){
  alpha = beta = 1 # alpha:shape, beta:scale of the margins
  MyCopula <- mvdc(copula=amhCopula(param=theta), # Ali-Mikail-Haq copula for (F(X), F(Y)), theta should be in [0,1]
                   margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                   paramMargins=list(shape=alpha,scale=beta) # alpha:shape, beta:scale
  )
}

for(m in 1:Max){
  
  #n = n_vect[i]
  
  
  #FunZDel = function(n,theta,lambda=lambda,alpha=2,beta=2){#,seed
  
  #  require(copula)
  #  require(survival)
  
  # MyCopula <- mvdc(copula=claytonCopula(param=theta), # Clayton copula for (F(X), F(Y))
  #                 margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
  #                 paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale
  
  
  set.seed(m)
  X <- rMvdc(n=n,MyCopula)
  X1 <- X[,1]
  X2 <- X[,2]
  
  # Censoring variable
  
  if(CopulaName=="Clayton"){
    C1 = rexp(n=n,rate = lambda)
    C2 = rexp(n=n,rate = lambda)
  }
  
  if(CopulaName=="AMH"){
    C1 = runif(n=n,min=0,max=a)
    C2 = runif(n=n,min=0,max=a)
  }
  
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
  #Z2 <- Z2[order(Z2)]
  del1 <- del1[-(n+1)]
  del2 <- del2[-(n+1)]
  #del2 <- del2[order(Z2)]
  
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
  if(CopulaName=="Clayton"){
    theta_hat_MassShift_vect[m] <- optimise(f=function(theta){logL_MassShift(theta,CopulaName = "Clayton")},interval=c(0,50),maximum = TRUE)$maximum
    theta_hat_ShihLouis_vect[m] <- optimise(f=function(theta){logL_ShihLouis(theta,CopulaName = "Clayton")},interval=c(0,50),maximum = TRUE)$maximum
  }
  
  if(CopulaName=="AMH"){
    theta_hat_MassShift_vect[m] <- optimise(f=function(theta){logL_MassShift(theta,CopulaName = "AMH")},interval=c(-1,1),maximum = TRUE)$maximum
    theta_hat_ShihLouis_vect[m] <- optimise(f=function(theta){logL_ShihLouis(theta,CopulaName = "AMH")},interval=c(-1,1),maximum = TRUE)$maximum
  }
  
  #theta_hat_MassShift_vect[m] <- optimise(f=logL_MassShift,interval=c(-1,1),maximum = TRUE)$maximum
  #theta_hat_MassShift_vect[m] <- newtonRaphson(fun=Score_MassShift, x0=Theta0)$root
  #theta_hat_MassShift_vect[m] <- uniroot(Score_MassShift, interval=c(0.01,50))$root
  #theta_hat_MassShift_vect[m] <- nlm(f=L_MassShift, p=Theta0)$estimate
  #theta_hat_ShihLouis_vect[m] <- optimise(f=L_ShihLouis,interval=c(0,50),maximum = TRUE)$maximum
  #theta_hat_ShihLouis_vect[m] <- optimise(f=logL_ShihLouis,interval=c(-1,1),maximum = TRUE)$maximum
  
  mean_del[m] = mean(del)
  
}

# 100*(1-mean(mean_del)) # censoring proportion

# c(mean(mean_del)*100,mean(theta_hat_MassShift_vect),mean(theta_hat_ShihLouis_vect),theta,mean(VarTheta))

# # Plot of copulas

#persp(claytonCopula(param=1,dim=2),FUN=dCopula,theta = 30,phi = 30) #
#plot(X1,X2)

Xlim = c(0,ceiling(max(X1)))
Ylim = c(0,ceiling(max(X2)))
Xlabels = seq(from=0,to=ceiling(max(X1)),by=1)
Ylabels = seq(from=0,to=ceiling(max(X2)),by=1)

plot(X1,X2,xlim=Xlim,ylim=Ylim,xlab="",ylab="",xaxt="none",yaxt="none")

axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, expression(x[1]), adj=0.6, font=2,cex=1.5)
axis(2, at=Ylabels,labels=Ylabels,las=0,font=2,hadj=1,padj=0)
mtext(side=2, line=2, expression(x[2]), adj=0.5, font=2,cex=1.5)

# # Plots for comparing the two log-likelihoods

if(CopulaName=="Clayton"){
  MinTheta <- 0
  MaxTheta <- 10
}

if(CopulaName=="AMH"){
  MinTheta <- -1
  MaxTheta <- 1
}

Theta <- matrix(seq(from=MinTheta,to=MaxTheta,by=0.01),ncol=1)

if(CopulaName=="Clayton"){
  Y_MassShift <- apply(X=Theta,MARGIN=1,FUN=function(theta){logL_MassShift(theta,CopulaName = "Clayton")})
  Y_ShihLouis <- apply(X=Theta,MARGIN=1,FUN=function(theta){logL_ShihLouis(theta,CopulaName = "Clayton")})
}

if(CopulaName=="AMH"){
  Y_MassShift <- apply(X=Theta,MARGIN=1,FUN=function(theta){logL_MassShift(theta,CopulaName = "AMH")})
  Y_ShihLouis <- apply(X=Theta,MARGIN=1,FUN=function(theta){logL_ShihLouis(theta,CopulaName = "AMH")})
}

#c(length(X),length(Y))
#Ylim <- range(c(Y_MassShift[is.finite(Y_MassShift)],Y_ShihLouis[is.finite(Y_ShihLouis)]))

Xlim = c(MinTheta,MaxTheta)
Ylim = range(c(Y_ShihLouis[is.finite(Y_ShihLouis)]))
Xlabels = c(seq(from=MinTheta,to=MaxTheta,length.out=6))
Ylabels = round(seq(from=Ylim[1],to=Ylim[2],length.out=6),digits=0)

plot(Theta,Y_ShihLouis,type='l',xlim=Xlim,ylim=Ylim,xlab="",ylab="",xaxt="none",yaxt="none",lwd = 2)
abline(v=c(theta_hat_ShihLouis_vect[m],theta),col=c('green','red'),lwd = 2,lty=c(2,2))

mtext(side=1, line=0.8, expression(hat(theta)), at=theta_hat_ShihLouis_vect[m], font=2,cex=1.5,col='green')
mtext(side=1, line=0.7, expression(theta[0]), at=theta, font=2,cex=1.5,col='red')

axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, expression(theta), adj=0.5, font=2,cex=1.5)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2,hadj=1,padj=0)

legend(x=0,y=-150,legend=c("log-likelihood","MLE","(true) parameter value"),lwd = 3,col=c("black","green", "red"),lty=c(1,2,2),cex=1.25,bty="n")


#plot(Theta,Y_MassShift,type='l') # ,ylim=Ylim
#plot(Theta,Y_ShihLouis,type='l') # ,ylim=Ylim
#lines(Theta,Y_ShihLouis,col="red") # ,ylim=Ylim
#abline(h=0,col="red")
#plot(Theta,Y_ShihLouis,type='l')

# MSE and its decomposition

Var_MassShift = mean((theta_hat_MassShift_vect-mean(theta_hat_MassShift_vect))**2)
Bias2_MassShift = (mean(theta_hat_MassShift_vect)-theta)**2

MSE_MassShift=Var_MassShift+Bias2_MassShift

Var_ShihLouis = mean((theta_hat_ShihLouis_vect-mean(theta_hat_ShihLouis_vect))**2)
Bias2_ShihLouis = (mean(theta_hat_ShihLouis_vect)-theta)**2

MSE_ShihLouis=Var_ShihLouis+Bias2_ShihLouis

MSE_MassShift/MSE_ShihLouis

# # Plots for comparing the two MSEs

Ylim <- range(c(theta_hat_MassShift_vect,theta_hat_ShihLouis_vect))
plot(theta_hat_MassShift_vect,col="black",ylim=Ylim)
points(theta_hat_ShihLouis_vect,col="blue")
abline(h=theta,col="red")

# # # Variance

# # Variance of An

# Fn1[Fn1==0] <- min(Fn1[Fn1!=0])/10 # values of 0 cause numerical issues
# Fn2[Fn2==0] <- min(Fn2[Fn2!=0])/10 # values of 0 cause numerical issues

#theta_hat <- log(theta_hat_MassShift_vect[m])
theta_hat <- theta_hat_MassShift_vect[m]

# Phi

if(CopulaName=="Clayton"){
  Phi <- rep(1/(theta_hat+1),times=n)-
    log(Fn1[-(n+1)]*Fn2[-(n+1)])+
    rep(1/theta_hat**2,times=n)*log(Fn1[-(n+1)]**rep(-theta_hat,times=n)+Fn2[-(n+1)]**rep(-theta_hat,times=n)-1)+
    rep(2+1/theta_hat,times=n)*(Fn1[-(n+1)]**rep(-theta_hat,times=n)*log(Fn1[-(n+1)])+Fn2[-(n+1)]**rep(-theta_hat,times=n)*log(Fn2[-(n+1)]))/(Fn1[-(n+1)]**rep(-theta_hat,times=n)+Fn2[-(n+1)]**rep(-theta_hat,times=n)-1)
}

if(CopulaName=="AMH"){
  Phi <- -3*(1-Fn1)*(1-Fn2)/(1-theta_hat*(1-Fn1)*(1-Fn2))+
    (2*(theta_hat-1)+(1-2*theta_hat)*(Fn1+Fn2)+(2*theta_hat+1)*Fn1*Fn2)/
    (1-2*theta_hat+theta_hat**2+theta_hat*(theta_hat-1)*(Fn1+Fn2)+theta_hat*(theta_hat+1)*Fn1*Fn2)
}

# NaN produced. Better avoid
# Phi[is.nan(Phi)] <- mean(Phi[!is.nan(Phi)]) # not helping. Variance is still negative.

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

Fun_PhiStar <- function(x,y,theta,CopulaName){
  if(CopulaName=="Clayton"){
    PhiStar <- -1/x-
      rep(1/theta,times=n)*x**rep(-theta-1,times=n)/(x**rep(-theta,times=n)+y**rep(-theta,times=n)-1)+
      rep(2+1/theta,times=n)*x**rep(-theta-1,times=n)*((x**rep(-theta,times=n)+y**rep(-theta,times=n)-1)*(-rep(theta,times=n)*log(x)+1)+
                                                         rep(theta,times=n)*(x**rep(-theta,times=n)*log(x)+y**rep(-theta,times=n)*log(y)))/
      (x**rep(-theta,times=n)+y**rep(-theta,times=n)-1)**2
  }
  if(CopulaName=="AMH"){
    PhiStar <- -3*(1-y)/(1-theta*(1-x)*(1-y))**2 +
      (1-2*theta+(2*theta+1)*y)/(1-2*theta+theta**2+theta*(1-theta)*(x+y)+theta*(theta+1)*x*y)-
      (theta*(2*(theta-1)+(1-2*theta)*(x+y)+(2*theta+1)*x*y)*((1-theta)*x+theta*(theta+1)*y))/
      (1-2*theta+theta**2+theta*(1-theta)*(x+y)+theta*(theta+1)*x*y)**2
  }
  return(PhiStar)
}

PhiStar1 = Fun_PhiStar(x=Fn1[-(n+1)],y=Fn2[-(n+1)],theta=theta_hat,CopulaName=CopulaName)
PhiStar2 = Fun_PhiStar(x=Fn2[-(n+1)],y=Fn1[-(n+1)],theta=theta_hat,CopulaName=CopulaName)

# Watch out! NaNs produced. Inputing mean values to the missing ones.
# PhiStar1[is.nan(PhiStar1)] <- mean(PhiStar1[!is.nan(PhiStar1)]) # not helping. Variance is still negative.
# PhiStar2[is.nan(PhiStar2)] <- mean(PhiStar2[!is.nan(PhiStar2)]) # not helping. Variance is still negative.

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