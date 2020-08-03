# # # # CHAPTER 2: ESTIMATION OF THE GENERATOR


# # Preps

library(copula) # for claytonCopula

n = 100

W = seq(from=1/n,to=1,by=0.01)

# # Simulating values from Archimedean copulas

# Independence copula

indep <- indepCopula(dim=2)
U <- rCopula(n=n, indep)[,1]
V <- rCopula(n=n, indep)[,2]

Phi_true = -log(W)

# Clayton copula

theta = 0.25
clayton <- claytonCopula(param=theta,dim=2)
U <- rCopula(n=n, clayton)[,1]
V <- rCopula(n=n, clayton)[,2]

Phi_true = (W**(-theta)-1)/theta

# AMH copula

theta = 0.5
AMH <- amhCopula(param=theta,dim=2)
U <- rCopula(n=n, AMH)[,1]
V <- rCopula(n=n, AMH)[,2]

Phi_true = log((1-theta*(1-W))/W)/(1-theta)

# Frank

theta = 1
frank <- frankCopula(param=theta,dim=2)
U <- rCopula(n=n, frank)[,1]
V <- rCopula(n=n, frank)[,2]

Phi_true = -((exp(theta)-1)/theta)*log((exp(-theta*W)-1)/(exp(-theta)-1))

# Gumbel (counter-example)

theta = 10
gumbel <- gumbelCopula(param=theta,dim=2)
U <- rCopula(n=n, gumbel)[,1]
V <- rCopula(n=n, gumbel)[,2]

Phi_true = (-log(W))**theta

# # Weibull distribution for X and Y

alpha = 10 # Weibull distribution shape parameter
beta = 1.7 # Weibull distribution scale parameter

X = beta*(-log(U))**(1/alpha)
Y = beta*(-log(V))**(1/alpha)

# # Pseudo-sample

V = rep(NA,times=n)

for(i in 1:n){
  V[i] = (sum((X <= X[i])*(Y <= Y[i])))/(n-1)
}

# # Empirical distribution

Kn = ecdf(V)
#Vn = unique(V[order(V)])
Vn = V[order(V)]

# # # Estimator 1: differential equation

Phi_est = rep(NA,times=length(W))

for(i in 1:length(W)){
  
  if(W[i]< 1-1/n){
    T = unique(c(W[i],Vn[which((W[i] <= Vn) & (Vn <= 1-1/n))],1-1/n))
    ProdTerm = rep(NA,times=length(T)-1)
    for(j in 1:(length(T)-1)){
      if(Kn(T[j])!=T[j]){
        ProdTerm[j] = (Kn(T[j])-T[j])/(Kn(T[j])-T[j+1])
      }else{
        ProdTerm[j] = 1
      }
    }
    Phi_est[i] = abs(prod(ProdTerm))/n
    }else{
      Phi_est[i] = 0
  }
}

# # Graph

Ylim = range(c(Phi_est,Phi_true),na.rm = TRUE)

plot(W,Phi_est,type='l',col='red',ylim=Ylim) #
lines(W,Phi_true,col='blue')

# # # Estimator 2: integral equation

h = rep(NA,times=n)
h[n] = 1 

for(i in (n-1):1){
  if(Kn(i/n)-i/n > 0){
  h[i] = (1/(n*Kn(i/n)-i))*sum(h[(i+1):n])
  }else{
    h[i] = h[i+1]
  }
}

Phi_est = rep(NA,times=length(W))

for(i in 1:length(W)){
  j0 = which.max((1:n/n <= W[i]) & (W[i] < (1:n+1)/n))
  Phi_est[i] = h[j0]*(1-W[i])+((n/2)*(1-W[i]**2)-i*(1-W[i]))*(h[j0+1]-h[j0])
}

Ylim = range(c(Phi_est,Phi_true),na.rm = TRUE)

plot(W,Phi_est,type='l',col='red',ylim=Ylim)
lines(W,Phi_true,col='blue')