library(tabulizer)
library(MASS) # for ginv
library(dplyr)
library(copula) # for claytonCopula


FunZDel = function(n,theta=2,lambda = 1/4){
  
  # Uniform variables with a Clayton copula joint distribution
  
  clayton <- claytonCopula(param=theta)
  U <- rCopula(n=n, clayton)[,1]
  V <- rCopula(n=n, clayton)[,2]
  
  # Weibull distribution for X and Y
  alpha = 10 # Weibull distribution shape parameter
  beta = 1.7 # Weibull distribution scale parameter
  
  X = beta*(-log(U))**(1/alpha)
  Y = beta*(-log(V))**(1/alpha)
  
  # Censoring variable

  C = rexp(n=n,rate = lambda)
  
  # Observations
  
  ordre = order(X,Y)
  Z_ordered = cbind(X,Y)[ordre,]
  
  xinf=max(Z_ordered[,1],Z_ordered[,2])+1 # CREATING POINT AT INFINITY
  Z1=c(Z_ordered[,1],xinf)
  Z2=c(Z_ordered[,2],xinf)
  
  del1 = as.integer(X <= C)
  del2 = as.integer(Y <= C)
  
  del1=c(del1,1)
  del2=c(del2,1)
  del=del1*del2
  
  return(list(Z1,Z2,del))
}


# # #  Numerical application

# Function: copula::rCopula

n = 5

# # Estimator 

Z1 = FunZDel(n=n)[[1]]
Z2 = FunZDel(n=n)[[2]]
del = FunZDel(n=n)[[3]]

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





# # Kendall's tau




# # # Real data

# Epidemiology Data from 
# McGilchrist, C. A., and C. W. Aisbett. "Regression with Frailty in Survival Analysis."
# Biometrics, vol. 47, no. 2, 1991, pp. 461-466.

#################################### USE THIS:####################################
# For data (Z1, Z2, del1, del2), sample size = n

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
library(CASdatasets)
data(canlifins) # load the dataset
# 14,889 contracts where one annuitant is male and the other female

#write.csv(canlifins,file="C:/Users/djimd/OneDrive/Documents/Concordia - PhD/Thesis/canlifins.csv",row.names = FALSE)

Sample = sample(1:dim(canlifins)[1],size = 500)
canlifins = canlifins[Sample,]

c(
sum((canlifins$DeathTimeF == 0) & (canlifins$DeathTimeF == 0)), # number of doubly censored couples
sum((canlifins$DeathTimeM > 0) & (canlifins$DeathTimeF == 0)), # number of couples where only the woman is censored
sum((canlifins$DeathTimeM == 0) & (canlifins$DeathTimeF > 0)), # number of couples where only the man is censored
sum((canlifins$DeathTimeM > 0) & (canlifins$DeathTimeF > 0)) # number of doubly uncensored couples
)

# # Estimator 

n = dim(canlifins)[1]

canlifins[canlifins$DeathTimeM == 0,"DeathTimeM"] = max(canlifins$AnnuityExpiredM)
canlifins[canlifins$DeathTimeF == 0,"DeathTimeF"] = max(canlifins$AnnuityExpiredM)

ordre = order(canlifins$DeathTimeM,canlifins$DeathTimeF)
Z_ordered = cbind(canlifins$DeathTimeM,canlifins$DeathTimeF)[ordre,]

del1 = as.integer(canlifins$DeathTimeM > 0) 
del2 = as.integer(canlifins$DeathTimeF > 0)

xinf=max(Z_ordered[,1],Z_ordered[,2])+1 # CREATING POINT AT INFINITY
Z1=c(Z_ordered[,1],xinf)
Z2=c(Z_ordered[,2],xinf)

del1=c(del1,1)
del2=c(del2,1)
del=del1*del2

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


# # #  Pseudo-likelihood maximization: variance

# n = 1000

n_vect = seq(from=100,to=1000,by=10)

theta_hat = rep(NA,length(n_vect))

lambda = 1/4 # check
alpha = 10 # Weibull distribution shape parameter
beta = 1.7 # Weibull distribution scale parameter

for(i in 1:length(n_vect)){
  
  n = n_vect[i]
  Z1 = FunZDel(n=n)[[1]]
  Z2 = FunZDel(n=n)[[2]]
  del = FunZDel(n=n)[[3]]
  
  xinf=max(Z_ordered[,1],Z_ordered[,2])+1 # point at infinity
  Z1=c(Z_ordered[,1],xinf)
  Z2=c(Z_ordered[,2],xinf)
  
  del1=c(del1,1)
  del2=c(del2,1)
  del=del1*del2
  
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
  phat=(solve(A)%*%Fbar)#[-(n+1)] # weights
  
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
  Fn1= 1-MMinv%*%b # 1-round((MMinv%*%b),digits = 10)# avoiding floating points # 1-MMinv%*%b
  #Fn1=Fn1[-(n+1)]
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
  Fn2=1-MMinv%*%b # 1-round((MMinv%*%b),digits = 10)
  #Fn2=Fn2[-(n+1)]
  Fn2[Fn2<=10**(-5)]=min(Fn2[Fn2>=10**(-5)])
  
  # Log-pseudo-likelihood function
  DiffLogL = function(x){
    sum(phat*(1/(x+1)-log(Fn1*Fn2))+(1/x**2)*log(Fn1**(-x)+Fn2**(-x)-1)-(2+1/x)*(-Fn1**(-x)*log(Fn1)-Fn2**(-x)*log(Fn2))/(Fn1**(-x)+Fn2**(-x)-1))
  }
  
  # MLE estimate of theta
  theta_hat[i] = optimise(f=DiffLogL,interval=c(0,10**2),maximum = TRUE)$maximum
}


# # Variance of An

# Psi

Psi = 1/(theta_hat[i]+1)-log(Fn1*Fn2)+(1/theta_hat[i]**2)*log(Fn1**(-theta_hat[i])+Fn2**(-theta_hat[i])-1)-
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

Matrix3 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%(Psi*PhiStar1))
Matrix4 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%(Psi*PhiStar2))

r1_hat = ginv(Matrix1)%*%Matrix3
r2_hat = ginv(Matrix1)%*%Matrix4


