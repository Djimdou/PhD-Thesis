library(tabulizer)
library(MASS) # for ginv
library(dplyr)


# # #  Numerical application

# # Estimator 



# # Kendall's tau




# # # Real data

# Epidemiology Data from 
# McGilchrist, C. A., and C. W. Aisbett. "Regression with Frailty in Survival Analysis."
# Biometrics, vol. 47, no. 2, 1991, pp. 461-466.

#################################### USE THIS:####################################
# For data (Z1, Z2, del1, del2), sample size = n

location = 'C:/Users/djimd/OneDrive/Documents/Concordia - PhD/Thesis/McGilchrist_Aisbett-1991.pdf'

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

#Z1 = Z_ordered[,1]
#Z2 = Z_ordered[,2]
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

#PHatFun = function(x,y){
#  phat[((Z_ordered[,1] >= x)*(Z_ordered[,2] >= y))]
#}

FBarFun = function(x,y){
  sum(phat[((Z_ordered[,1] >= x)*(Z_ordered[,2] >= y))])
}

#Fn_func = function(x,y){
#  sum(phat[((Z_ordered[,1] <= x)*(Z_ordered[,2] <= y))])
#}

x = unique(Z_ordered[order(Z_ordered[,1]),1])
y = unique(Z_ordered[order(Z_ordered[,2]),2])
#F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))
Fn_grid <- outer(X=x,Y=y, FUN=Vectorize(Fn_func))

#P_hat_grid <- outer(X=x,Y=y, FUN=Vectorize(PHatFun))
#P_hat_grid*F_bar_grid

persp(x,y,F_bar_grid, theta = 30, phi = 30)

#Tau = sum(diff(c(as.vector(t(Fn_grid)),1))*as.vector(t(F_bar_grid)))
#Tau = cor(x=KidneyInfection$T1,y=KidneyInfection$T2,method="kendall")
Tau = phat %*% Fbar[-(n+1)]

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

#install.packages('xts')
library(xts)
#install.packages('sp')
library(sp)
#install.packages('zoo')
library(zoo)

install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/R/", type="source")
library(CASdatasets)
data(canlifins) # load the dataset
# 14,889 contracts where one annuitant is male and the other female

#write.csv(canlifins,file="C:/Users/djimd/OneDrive/Documents/Concordia - PhD/Thesis/canlifins.csv",row.names = FALSE)

canlifins$UncensoredM = as.integer((canlifins$AnnuityExpiredM >= canlifins$DeathTimeM) & (canlifins$DeathTimeM > 0))
canlifins$UncensoredF = as.integer((canlifins$AnnuityExpiredM >= canlifins$DeathTimeF) & (canlifins$DeathTimeF > 0))

Sample = sample(1:dim(canlifins)[1],size = 3000)
#AtLeastOneCensored = which((canlifins$AnnuityExpiredM > canlifins$DeathTimeM)|(canlifins$AnnuityExpiredM > canlifins$DeathTimeF))
canlifins = canlifins[Sample,]

c(
sum((canlifins$UncensoredM == 0) & (canlifins$UncensoredF == 0)), # number of doubly censored couples
sum((canlifins$UncensoredM != 0) & (canlifins$UncensoredF == 0)), # number of couples where only the woman is censored
sum((canlifins$UncensoredM == 0) & (canlifins$UncensoredF != 0)), # number of couples where only the man is censored
sum((canlifins$UncensoredM != 0) & (canlifins$UncensoredF != 0)) # number of doubly uncensored couples
)

# # Estimator 

n = dim(canlifins)[1]

#canlifins$UncensoredM = 1-as.integer(canlifins$DeathTimeM==0)
#canlifins$UncensoredF = 1-as.integer(canlifins$DeathTimeF==0)

canlifins[canlifins$DeathTimeM == 0,"DeathTimeM"] = max(canlifins$DeathTimeM)
canlifins[canlifins$DeathTimeF == 0,"DeathTimeF"] = max(canlifins$DeathTimeF)

ordre = order(canlifins$DeathTimeM,canlifins$DeathTimeF)
Z_ordered = cbind(canlifins$DeathTimeM,canlifins$DeathTimeF)[ordre,]

del1 = canlifins$UncensoredM 
del2 = canlifins$UncensoredF

xinf=max(Z_ordered[,1],Z_ordered[,2])+1 # CREATING POINT AT INFINITY
Z1=c(Z_ordered[,1],xinf)
Z2=c(Z_ordered[,2],xinf)

del1=c(del1,1)
del2=c(del2,1)
del=del1*del2

az1=matrix(rep(Z1,n+1),ncol=n+1)
az2=matrix(rep(Z2,n+1),ncol=n+1)
A=(t(az1)>=az1)*(t(az2)>=az2)
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

#PHatFun = function(x,y){
#  phat[((Z_ordered[,1] >= x)*(Z_ordered[,2] >= y))]
#}

FBarFun = function(x,y){
  sum(phat[((Z_ordered[,1] >= x)*(Z_ordered[,2] >= y))])
}

#Fn_func = function(x,y){
#  sum(phat[((Z_ordered[,1] <= x)*(Z_ordered[,2] <= y))])
#}

x = unique(Z_ordered[order(Z_ordered[,1]),1])
y = unique(Z_ordered[order(Z_ordered[,2]),2])
#F_bar_grid <- outer(X=x,Y=y, FUN=Vectorize(FBarFun))
Fn_grid <- outer(X=x,Y=y, FUN=Vectorize(Fn_func))

#P_hat_grid <- outer(X=x,Y=y, FUN=Vectorize(PHatFun))
#P_hat_grid*F_bar_grid

persp(x,y,F_bar_grid, theta = 30, phi = 30)

#Tau = sum(diff(c(as.vector(t(Fn_grid)),1))*as.vector(t(F_bar_grid)))
#Tau = cor(x=KidneyInfection$T1,y=KidneyInfection$T2,method="kendall")
Tau = phat %*% Fbar[-(n+1)]

# # Kendall's tau variance

Matrix1 = rbind(diag(1,n+1) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(as.vector(Fbar)),b%*%diag(as.vector(Fbar))))%*%(diag(1,n+1) + B%*%D%*%B)%*%(diag(as.vector(Fbar))%*%B%*%Fbar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = t(Fbar)%*%B%*%V%*%B%*%Fbar+
  2*t(Fbar)%*%B%*%W_hat+
  t(Fbar)%*%B%*%diag(as.vector(Fbar))%*%(diag(1,n+1)+B%*%D%*%B)%*%diag(as.vector(Fbar))%*%B%*%Fbar

