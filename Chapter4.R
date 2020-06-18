library(tabulizer)
library(MASS) # for ginv

#install.packages("dplyr")
library(dplyr)



# # #  Numerical application

# # Estimator 



# # Kendall's tau




# # # Real data

# Epidemiology Data from 
# McGilchrist, C. A., and C. W. Aisbett. "Regression with Frailty in Survival Analysis."
# Biometrics, vol. 47, no. 2, 1991, pp. 461-466.

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

#Z = apply(X=KidneyInfection[,c("T1","T2")], MARGIN=1, FUN=min)
ordre = order(KidneyInfection$T1,KidneyInfection$T2)
Z_ordered = cbind(KidneyInfection$T1,KidneyInfection$T2)[ordre,]

A = matrix(0,ncol=n+1,nrow=n+1)
A[,n+1] = 1
A[n+1,] = 1

for(i in 1:n){
  for(k in 1:n){
    A[i,k] = ifelse((Z_ordered[k,1] >= Z_ordered[i,1]) & (Z_ordered[k,1] >= Z_ordered[i,1]),1,0)
  }
}

b = (KidneyInfection$uncensored1*KidneyInfection$uncensored2)/(1:n+1)
B = diag(c(b,1))

c = b/(1-b)

p = rep(NA,n+1)
p[n+1] = 1/(n+1)

for(i in n:1){
  p[i] = c[i]*sum(A[i,(i+1):(n+1)]*p[(i+1):(n+1)])
}

F_bar = A%*%p

A = A[-(n+1),-(n+1)]
B = B[-(n+1),-(n+1)]
F_bar = F_bar[-(n+1)]

# Its covariance matrix

Matrix1 = rbind(diag(1,n) - A%*%B,-b)
Matrix2 = cbind(t(diag(1,n) - A%*%B),-b)
Matrix3 = (rbind(A%*%B%*%diag(F_bar),b%*%diag(F_bar)))%*%(diag(1,n) + B%*%D%*%B)%*%cbind(diag(F_bar)%*%B%*%A,diag(F_bar)%*%b)

V = ginv(Matrix1)%*%Matrix3%*%ginv(Matrix2)

# # Kendall's tau variance

D = matrix(NA,nrow=n,ncol=n)

for(i in 1:n){
  for(j in 1:n){
    D[i,j] = (1-A[i,j])*(1-A[j,i])*(A[i,]%*%A[,j])
  }
}

Matrix1 = rbind(diag(1,n) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(F_bar),b%*%diag(F_bar)))%*%(diag(1,n) + B%*%D%*%B)%*%(diag(F_bar)%*%B%*%F_bar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = F_bar%*%B%*%V%*%B%*%F_bar+
        2*F_bar%*%B%*%W_hat+
        F_bar%*%B%*%diag(F_bar)%*%(diag(1,n)+B%*%D%*%B)%*%diag(F_bar)%*%B%*%F_bar


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


# # Estimator 

n = dim(canlifins)[1]

canlifins$UncensoredM = 1-as.integer(canlifins$DeathTimeM==0)
canlifins$UncensoredF = 1-as.integer(canlifins$DeathTimeF==0)

#canlifins$DurationM = canlifins$EntryAgeM + canlifins$DeathTimeM
#canlifins$DurationF = canlifins$EntryAgeF + canlifins$DeathTimeF

ordre = order(canlifins$DeathTimeM,canlifins$DeathTimeF)
Z_ordered = cbind(canlifins$DeathTimeM,canlifins$DeathTimeF)[ordre,]

A = matrix(0,ncol=n+1,nrow=n+1)
A[,n+1] = 1
A[n+1,] = 1

for(i in 1:n){
  for(k in 1:n){
    A[i,k] = ifelse((Z_ordered[k,1] >= Z_ordered[i,1]) & (Z_ordered[k,1] >= Z_ordered[i,1]),1,0)
  }
}

b = (canlifins$UncensoredM*canlifins$UncensoredF)/(1:n+1)
B = diag(c(b,1))

c = b/(1-b)

p = rep(NA,n+1)
p[n+1] = 1/(n+1)

for(i in n:1){
  p[i] = c[i]*sum(A[i,(i+1):(n+1)]*p[(i+1):(n+1)])
}

F_bar = A%*%p

A = A[-(n+1),-(n+1)]
B = B[-(n+1),-(n+1)]
F_bar = F_bar[-(n+1)]

# Its covariance matrix

Matrix1 = rbind(diag(1,n) - A%*%B,-b)
Matrix2 = cbind(t(diag(1,n) - A%*%B),-b)
Matrix3 = (rbind(A%*%B%*%diag(F_bar),b%*%diag(F_bar)))%*%(diag(1,n) + B%*%D%*%B)%*%cbind(diag(F_bar)%*%B%*%A,diag(F_bar)%*%b)

V = ginv(Matrix1)%*%Matrix3%*%ginv(Matrix2)

# # Kendall's tau variance

D = matrix(NA,nrow=n,ncol=n)

for(i in 1:n){
  for(j in 1:n){
    D[i,j] = (1-A[i,j])*(1-A[j,i])*(A[i,]%*%A[,j])
  }
}

Matrix1 = rbind(diag(1,n) - A%*%B,-b)
Matrix2 = (rbind(A%*%B%*%diag(F_bar),b%*%diag(F_bar)))%*%(diag(1,n) + B%*%D%*%B)%*%(diag(F_bar)%*%B%*%F_bar)

W_hat = ginv(Matrix1)%*%Matrix2

V_tau = F_bar%*%B%*%V%*%B%*%F_bar+
  2*F_bar%*%B%*%W_hat+
  F_bar%*%B%*%diag(F_bar)%*%(diag(1,n)+B%*%D%*%B)%*%diag(F_bar)%*%B%*%F_bar