library(tabulizer)
library(MASS) # for ginv
library(dplyr)
library(copula) # for claytonCopula
#install.packages('pracma')
library(pracma) # for NewtonRaphson procedure
library(CASdatasets) # dataset 



FunZDel = function(n,theta,lambda=1,seed,alpha=2,beta){
  
  # Uniform variables with a Clayton copula joint distribution
  
  #clayton <- claytonCopula(param=theta)
  #U <- rCopula(n=n, clayton)[,1]
  #V <- rCopula(n=n, clayton)[,2]
  
  # Weibull distribution for X and Y
  #alpha = 10 # Weibull distribution shape parameter
  #beta = 1.7 # Weibull distribution scale parameter
  
  MyCopula <- mvdc(copula=claytonCopula(param=theta), # Clayton copula for (F(X), F(Y))
                   margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                   paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale
  
  set.seed(seed)
  XY <- rMvdc(n=n,MyCopula)
  X <- XY[,1]
  Y <- XY[,2]
  
  #X = beta*(-log(U))**(1/alpha)
  #Y = beta*(-log(V))**(1/alpha)
  
  # Censoring variable

  Cx = rexp(n=n,rate = lambda)
  Cy = rexp(n=n,rate = lambda)
  
  # Observations
  
  Z1 <- pmin(X,Cx)
  Z2 <- pmin(Y,Cy)
  
  ordre = order(Z1,Z2)
  
  xinf=max(Z1,Z2)+1 # CREATING POINT AT INFINITY
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
LogL = function(x){
  sum(phat*(log(x+1)-(x+1)*log(Fn1*Fn2)-(2+1/x)*log(Fn1**(-x)+Fn2**(-x)-1)))
}

# Log-pseudo-likelihood function derivative
DiffLogL = function(x){
  sum(phat*(rep(1/(x+1),times=n+1)-log(Fn1*Fn2)+
              (rep(1/x**2,times=n+1))*log(Fn1**(rep(-x,times=n+1))+
                                            Fn2**(rep(-x,times=n+1))-1)+
              (rep(2+1/x,times=n+1))*(Fn1**(rep(-x,times=n+1))*
                                        log(Fn1)+Fn2**(rep(-x,times=n+1))*log(Fn2))/
              (Fn1**(rep(-x,times=n+1))+Fn2**(rep(-x,times=n+1))-1)))
}

#DiffLogL2 = function(x){
#  sum(phat*(1/(x+1)-log(Fn1*Fn2)+(1/x**2)*log(Fn1**(-x)+Fn2**(-x)-1)
#            +(2+1/x)*(Fn1**(-x)*log(Fn1) +  Fn2**(-x)*log(Fn2))/(Fn1**(-x)+Fn2**(-x)-1)))
#}


FunZDel_canlifins = function(size,seed=1){
  
  data(canlifins) # load the dataset
  # 14,889 contracts where one annuitant is male and the other female
  
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
  
  del = c((canlifins$DeathTimeM > 0) * (canlifins$DeathTimeF > 0),1)
  
  #canlifins[canlifins$DeathTimeM == 0,"DeathTimeM"] = canlifins[canlifins$DeathTimeM == 0,"AnnuityExpiredM"]
  #canlifins[canlifins$DeathTimeF == 0,"DeathTimeF"] = canlifins[canlifins$DeathTimeF == 0,"AnnuityExpiredM"]
  
  ordre = order(canlifins$DeathTimeM,canlifins$DeathTimeF)
  Z_ordered = cbind(canlifins$DeathTimeM,canlifins$DeathTimeF)[ordre,]
  
  del=del[ordre]
  
  xinf=max(Z_ordered[,1],Z_ordered[,2])+1
  Z1=c(Z_ordered[,1],xinf)
  Z2=c(Z_ordered[,2],xinf)
  
  return(list(Z1,Z2,del))
}


# # # Kendall's tau:simulated data

n = 1000
Max = 1000 # number of samples for MSE 
Tau_hat = rep(NA,times=Max)
del_mean = rep(NA,times=Max)
theta = 2
beta=1.1

for(m in 1:Max){

  # # Estimator 
  
  ZDel = FunZDel(n=n,theta=theta,seed=m,beta=beta)
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
  Fbar=MMinv%*%b
  
  # The function 'solve' has some trouble sinverting some A's for n >= 1500
  # (example with seed=4). 'ginv' can do it, but is more time consuming. 
  # A tradeoff is to use ginv where solve fails. 
  
  if("try-error" %in% class(try(solve(A)))){
    phat=(ginv(A)%*%Fbar)
  } else {
    phat=(solve(A)%*%Fbar)#[-(n+1)]
  }
  
  # # Kendall tau estimate
  
  Tau_hat[m] = 4*(t(phat) %*% Fbar)-1

}


# MSE and its decomposition

Tau = theta/(theta+2)
#plot(Tau_hat,type='l',ylim=range(c(Tau_hat,Tau)));abline(h=Tau);

Xlabels = seq(from=0,to=Max)
Ylabels = seq(from=0,to=1,by=0.2)
Ylim = c(0,1)
par(mfrow=c(1,1))
plot(Tau_hat,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none",ylim=Ylim)
abline(h=Tau,col='red',lwd = 3)
abline(h=mean(Tau_hat),col='black',lwd = 3)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "iteration rank", font=2,cex=1.5)

Var = var(Tau_hat)*((Max-1)/Max)
Bias2 = (mean(Tau_hat)-Tau)**2

MSE=Var+Bias2

# Fbar Variance

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf[1:(n+1)])
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# # Kendall's tau variance

Matrix1 = rbind(diag(1,n+1) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = 4**2*(t(Fbar)%*%B%*%V%*%B%*%Fbar+
  2*t(Fbar)%*%B%*%W_hat+
  t(Fbar)%*%B%*%diag(as.vector(Fbar))%*%(diag(1,n+1)+B%*%D%*%B)%*%diag(as.vector(Fbar))%*%B%*%Fbar)

# # Graphic for Fbar 

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


location = 'C:/Users/djimd/OneDrive/Documents/Concordia - PhD/Thesis/McGilchrist_Aisbett-1991.pdf'
# Data: https://rdrr.io/cran/SurvCorr/man/kidney.html

# Extract the table
mydata <- extract_tables(file=location,pages=6)

data_mc = data.frame(mydata)[-(1:3),2:3]

# strsplit(toString(data.frame(mydata)[-(1:3),2]),",")

times = matrix(as.numeric(unlist(strsplit(toString(data.frame(mydata)[-(1:3),2]),","))),ncol = 2,byrow =TRUE)

#unlist(strsplit(substring(data.frame(mydata)[-(1:3),3],1,4),","))

censor = substring(data.frame(mydata)[-(1:3),3],1,4)
censor[7] = sub(" ",",",censor[7],fixed=TRUE)
censor[34] = sub("?",",",censor[34],fixed=TRUE)

censoring = matrix(as.numeric(unlist(strsplit(censor,","))),ncol = 2,byrow =TRUE)

KidneyInfection = cbind.data.frame(times,censoring)
colnames(KidneyInfection) = c("T1","T2","uncensored1","uncensored2")

# # Estimator 

n = dim(KidneyInfection)[1]

ordre = order(KidneyInfection$T1,KidneyInfection$T2)#
Z_ordered = cbind(KidneyInfection$T1,KidneyInfection$T2)[ordre,]

del1 = KidneyInfection$uncensored1 
del2 = KidneyInfection$uncensored2

xinf=max(Z_ordered[,1],Z_ordered[,2])+1 # CREATING POINT AT INFINITY
Z1=c(Z_ordered[,1],xinf)
Z2=c(Z_ordered[,2],xinf)

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
Fbar=MMinv%*%b ### <---- THIS IS THE \bar{F} VECTOR
phat=(solve(A)%*%Fbar)[-(n+1)] # weights

# Variance estimator

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf[1:(n+1)])
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# # Kendall's tau estimate

FBarFun = function(x,y){
  sum(phat[((Z_ordered[,1] >= x)*(Z_ordered[,2] >= y))])
}

x = unique(Z_ordered[order(Z_ordered[,1]),1])
y = unique(Z_ordered[order(Z_ordered[,2]),2])
F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))

persp(x,y,F_bar_grid, theta = 30, phi = 30)

Tau = 4*(phat %*% Fbar[-(n+1)])-1

# # Kendall's tau variance

Matrix1 = rbind(diag(1,n+1) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = t(Fbar)%*%B%*%V%*%B%*%Fbar+
        2*t(Fbar)%*%B%*%W_hat+
        t(Fbar)%*%B%*%diag(as.vector(Fbar))%*%(diag(1,n+1)+B%*%D%*%B)%*%diag(as.vector(Fbar))%*%B%*%Fbar





# # # Insurance data from 
# # Frees, Edward W., et al. "Annuity Valuation with Dependent Mortality." 
# The Journal of Risk and Insurance, vol. 63, no. 2, 1996, pp. 229-261.

# Link: http://cas.uqam.ca/

#library(xts)
#library(sp)
#library(zoo)

#install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/R/", type="source")


# Estimator

n=10000

ZDel = FunZDel_canlifins(size=n)
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
Fbar=MMinv%*%b ### <---- THIS IS THE \bar{F} VECTOR
phat=(solve(A)%*%Fbar)[-(n+1)] # weights

# Graph of Fbar

FBarFun = function(x,y){
  sum(phat[((Z_ordered[,1] >= x)*(Z_ordered[,2] >= y))])
}

x = unique(Z_ordered[order(Z_ordered[,1]),1])
y = unique(Z_ordered[order(Z_ordered[,2]),2])
F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))

persp(x,y,F_bar_grid, theta = 30, phi = 30)

# Variance estimator

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf[1:(n+1)])
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# # Kendall's tau estimate

Tau_hat = 4*(phat %*% Fbar[-(n+1)])-1

# # Kendall's tau variance

Matrix1 = rbind(diag(1,n+1) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = t(Fbar)%*%B%*%V%*%B%*%Fbar+
  2*t(Fbar)%*%B%*%W_hat+
  t(Fbar)%*%B%*%diag(as.vector(Fbar))%*%(diag(1,n+1)+B%*%D%*%B)%*%diag(as.vector(Fbar))%*%B%*%Fbar


# # #  Pseudo-likelihood maximization

# n = 5

n_vect = seq(from=100,to=1500,by=100)

theta_hat = rep(NA,length(n_vect))

#lambda = 1/4 # check

theta = 2

for(i in 1:length(n_vect)){
  
  n = n_vect[i]
  ZDel = FunZDel(n=n,theta=theta,seed=i)
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
  Fbar=MMinv%*%b 
  phat=(solve(A)%*%Fbar)
  
  # Starting point of Newton-Raphson
  
  Tau = 4*(t(phat[-(n+1)]) %*% Fbar[-(n+1)])-1
  Theta0 = 2*Tau/(1-Tau)
  
  # # MLE of theta
  
  # Fn1
  
  az1=matrix(rep(Z1,n+1),ncol=n+1)
  A=t(az1)>=az1
  
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)

  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)

  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fn1= 1-MMinv%*%b 
  Fn1[Fn1<=10**(-5)]=min(Fn1[Fn1>=10**(-5)])
  
  # Fn2
  
  az2=matrix(rep(Z2,n+1),ncol=n+1)
  A=t(az2)>=az2
  
  A[lower.tri(A, diag = FALSE)] = 0
  rsumA=apply(A,1,sum)
  
  b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
  B=diag(b)
  
  M=rbind(Id-A%*%B,-t(b))
  MMinv=solve(t(M)%*%M)
  Fn2=1-MMinv%*%b
  Fn2[Fn2<=10**(-5)]=min(Fn2[Fn2>=10**(-5)])

  # x=matrix(seq(from=0.01,to=5,by=0.1)); ForPlot=apply(X=x,MARGIN=1,FUN=DiffLogL);
  # plot(x,ForPlot,type='l',ylim=range(c(ForPlot,0)));abline(h=0);
  
  # MLE estimate of theta
  # theta_hat[i] = optimise(f=LogL,interval=c(0,10**2),maximum = TRUE)$maximum
  theta_hat[i] = newtonRaphson(fun=DiffLogL, x0=Theta0)$root
  #theta_hat[i] = uniroot(f=DiffLogL2, interval=c(0.01,10**2))$root
  
  # Coding Newton-Raphson:
  # https://rpubs.com/aaronsc32/newton-raphson-method#:~:text=%23%23%20%5B1%5D%203.162278-,Newton%2DRaphson%20Method%20in%20R,rootSolve%20package%20features%20the%20uniroot.
  
}

# plot(theta_hat,type='l',ylim=range(c(theta_hat,theta)));abline(h=theta);

# # Variance of An

# Phi

Phi = 1/(theta_hat[i]+1)-log(Fn1*Fn2)+(1/theta_hat[i]**2)*log(Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1)-
  (2+1/theta_hat[i])*(-Fn1**(-theta_hat[i])*log(Fn1)-Fn2**(-theta_hat[i])*log(Fn2))/(Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1)

# VStar

az1=matrix(rep(Z1,n+1),ncol=n+1)
A=(t(az1)>=az1)
A[lower.tri(A, diag = FALSE)] = 0
rsumA=apply(A,1,sum)

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=MMinv%*%b
phat=(solve(A)%*%Fbar)[-(n+1)]

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf[1:(n+1)])
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V_z0=(MMinv%*%U)%*%MMinv

az2=matrix(rep(Z2,n+1),ncol=n+1)
A=(t(az2)>=az2)
A[lower.tri(A, diag = FALSE)] = 0
rsumA=apply(A,1,sum)

b=c((1-eps)*del[1:n]/((1-eps)*(rsumA[1:n]-1)+eps*n),1)
B=diag(b)

M=rbind(Id-A%*%B,-t(b))
MMinv=solve(t(M)%*%M)
Fbar=MMinv%*%b ### <---- THIS IS THE \bar{F} VECTOR
phat=(solve(A)%*%Fbar)[-(n+1)] # weights

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf[1:(n+1)])
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V_0z=(MMinv%*%U)%*%MMinv

PhiStar1  = -1/Fn1-(1/theta_hat[i])*(Fn1**(-theta_hat[i]-1)/(Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1))+
  (2+1/theta_hat[i])*Fn1**(-theta_hat[i]-1)*((Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1)*(-theta_hat[i]*log(Fn1)+1)+
  theta_hat[i]*Fn1**(-theta_hat[i])*log(Fn1))/(Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1)**2

PhiStar2  = -1/Fn2-(1/theta_hat[i])*(Fn2**(-theta_hat[i]-1)/(Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1))+
  (2+1/theta_hat[i])*Fn1**(-theta_hat[i]-1)*((Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1)*(-theta_hat[i]*log(Fn2)+1)+
                                               theta_hat[i]*Fn1**(-theta_hat[i])*log(Fn2))/(Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1)**2

PhiStar1_matrix = matrix(rep(PhiStar1,times=n+1),ncol=n+1)
PhiStar2_matrix = matrix(rep(PhiStar2,times=n+1),ncol=n+1)

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
Fbar=MMinv%*%b
phat=(solve(A)%*%Fbar)

D=(1-A)*(1-t(A))*(A%*%t(A))

Matrix1 = rbind(diag(1,n+1) - A%*%B,-b)

Matrix3 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%(Phi*PhiStar1))
Matrix4 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%(Phi*PhiStar2))

r1_hat = ginv(Matrix1)%*%Matrix3
r2_hat = ginv(Matrix1)%*%Matrix4

# e_hat

Matrix5 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%(Phi*Phi))

e_hat = ginv(Matrix1)%*%Matrix5

# V (variance of estimator)

D=(1-A)*(1-t(A))*(A%*%t(A))
bf=b*Fbar
BF=diag(bf[1:(n+1)])
S=rbind(A%*%BF,t(bf))
R=S%*%(Id+((B%*%D)%*%B))%*%t(S)
U=(t(M)%*%R)%*%M
V=(MMinv%*%U)%*%MMinv

# Var(An)

Var_An = t(PhiStar1)%*%diag(as.vector(Fbar))%*%B%*%V%*%B%*%diag(as.vector(Fbar))%*%PhiStar2+
  2*b%*%VStar%*%diag(as.vector(Fbar))%*%b+
  2*t(r1_hat+r2_hat)%*%diag(as.vector(Fbar))%*%b+
  t(Phi)%*%B%*%V%*%B%*%Phi+
  2*b%*%e_hat+
  t(Phi)%*%B%*%diag(as.vector(Fbar))%*%(diag(as.vector(Fbar))%*%b+B%*%D%*%B%*%diag(as.vector(Fbar))%*%B%*%Phi)

# # Variance of sqrt(n)(theta_hat - theta)

VarTheta = Var_An/(t(Phi)**2%*%Fbar)