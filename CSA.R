# # # # CSA paper: Distribution of the joint survival function of an Archimedean copula

# # # Preliminary

#install.packages("GoFKernel")
library(GoFKernel) # for inverse
#install.packages("copula")
library(copula) # for claytonCopula

# # Goodness-of-fit?

# Simulating values

# m = 1000

# W = 1:m/m # grid for the x-axis

Theta0 = 4

Theta = seq(from=0.1,to=5,by=0.1) # c(0.1,1,3,5)  # 7 seems not good a value for the numerical approximations later

A = seq(from=0.01,to=1/2,by=0.01)
Z = seq(from=0.01, to=0.99,by=0.01)
M = 10**15

cop <- claytonCopula(param=Theta0,dim=2)

# Sampling X and Y

# n = 1000 # size of sample

alpha = 2 # Weibull distribution shape parameter
beta = 2 # Weibull distribution scale parameter

MyCopula <- mvdc(copula=cop, # copula for (F(X), F(Y))
                 margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                 paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale

#set.seed(1)
XY <- rMvdc(n=length(Z),MyCopula)
X <- XY[,1]
Y <- XY[,2]

U1 = pweibull(X, shape=alpha, scale = beta, lower.tail = TRUE, log.p = FALSE)
U2 = pweibull(Y, shape=alpha, scale = beta, lower.tail = TRUE, log.p = FALSE)

Cn_star = rep(NA,times=length(Z))

for(i in 1:length(Z)){
  Cn_star[i] = (sum((U1 <= U1[i])*(U2 <= U2[i])))/length(Z)
}

# Empirical distribution of Cn_star

Kn_star = ecdf(Cn_star)

# # Distance

# Cn_star = Cn_star[order(Cn_star)]

# Kstar

inv_a_z = integrand = diff_a_z = matrix(NA,ncol=length(Z),nrow=length(A)) 

k_est = matrix(NA,nrow=length(Theta),ncol=length(Z)) # density of C_star

for(t in 1:length(Theta)){
  for(i in 1:length(A)){
    psi_inverse = inverse(function (z) {1-(Theta[t]*A[i]*z+1)**(-1/Theta[t])-(Theta[t]*(1-A[i])*z+1)**(-1/Theta[t])+(Theta[t]*z+1)**(-1/Theta[t])}, 0, M)
    for(j in 1:length(Z)){
      inv_a_z[i,j] = unlist(psi_inverse(Z[j]))
      diff_a_z[i,j] = 1/(A[i]*(Theta[t]*A[i]*inv_a_z[i,j]+1)**(-1/Theta[t]-1)+(1-A[i])*(Theta[t]*(1-A[i])*inv_a_z[i,j]+1)**(-1/Theta[t]-1)-(Theta[t]*inv_a_z[i,j]+1)**(-1/Theta[t]-1))
      integrand[i,j] = (Theta[t]+1)*inv_a_z[i,j]*(Theta[t]*inv_a_z[i,j]+1)**(-2-1/Theta[t]) * diff_a_z[i,j]
    }
  }
  
  for(j in 1:length(Z)){
    k_est[t,j] = 2*integrate(approxfun(x=A,y=integrand[,j]), lower = min(A), upper = max(A))[[1]]
  }
}

K_star = matrix(NA,nrow=length(Theta),ncol=length(Z))  # distribution of C_star
K_star[,1] = rep(0,times=length(Theta))

for(t in 1:length(Theta)){
  for(j in 2:length(Z)){
    K_star[t,j] = integrate(approxfun(x=Z,y=k_est[t,]), lower = min(Z), upper = Z[j])[[1]]
  }
}

# plot(Z,K_star[40,],type="l")

dn = rep(NA,times=length(Theta))

for(t in 1:length(Theta)){
  dn[t] = max(abs(Kn_star(Z) - K_star[t,]))
}

# plot(Z,Kn_star(Z),type="l",col="blue");lines(Z,K_star[10,],col="green")

plot(Theta,dn,type="l")
abline(v=Theta0,col='red')