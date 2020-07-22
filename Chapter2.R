# # CHAPTER 2: ESTIMATION OF THE GENERATOR


# Librairies

library(copula) # for claytonCopula


n = 200

# Clayton copula

clayton <- claytonCopula(param=2)
U <- rCopula(n=n, clayton)[,1]
V <- rCopula(n=n, clayton)[,2]

# Weibull distribution for X and Y
alpha = 10 # Weibull distribution shape parameter
beta = 1.7 # Weibull distribution scale parameter

X = beta*(-log(U))**(1/alpha)
Y = beta*(-log(V))**(1/alpha)

# Pseudo-sample

V = rep(NA,times=n)

for(i in 1:n){
  V[i] = (sum((X <= X[i])*(Y <= Y[i])))/(n-1)
}

# Empirical distribution

Kn = ecdf(V)
Vn = unique(V[order(V)])

# # Estimator 1: differential equation

W = seq(from=0,to=1,by=0.01)
Phi = rep(NA,times=length(W))

for(i in 1:length(W)){
  
  if(W[i]<=1-1/n){
    T = unique(c(W[i],Vn[which((W[i] <= Vn) & (Vn <= 1-1/n))],1-1/n))
    ProdTerm = rep(NA,times=length(T)-1)
    for(j in 1:(length(T)-1)){
      if(Kn(T[j])!=T[j]){
        ProdTerm[j] = (Kn(T[j])-T[j])/(Kn(T[j])-T[j+1])
      }else{
        ProdTerm[j] = 1
      }
    }
    Phi[i] = abs(prod(ProdTerm))/n
    }else{
    Phi[i] = 0
  }
}

plot(W,Phi,type='l')

# # Estimator 2: integral equation