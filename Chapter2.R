# # # # CHAPTER 2: ESTIMATION OF THE GENERATOR


# # Preps: loading libraries

#install.packages("copula")
library(copula) # for claytonCopula


# # Simulating values from Archimedean copulas

m = 1000

W = 1:m/m #seq(from=1/n,to=1,by=0.01) # grid for the x-axis

# Independence copula

cop <- indepCopula(dim=2)
Phi_true = -log(W)

# Clayton copula

theta = 2
cop <- claytonCopula(param=theta,dim=2)
Phi_true = (W**(-theta)-1)/theta

# AMH copula

theta = 0.5
cop <- amhCopula(param=theta,dim=2)
Phi_true = log((1-theta*(1-W))/W)/(1-theta)

# Frank

theta = -1
cop <- frankCopula(param=theta,dim=2)
Phi_true = -((exp(theta)-1)/theta)*log((exp(-theta*W)-1)/(exp(-theta)-1))

# Gumbel (counter-example)

theta = 10
cop <- gumbelCopula(param=theta,dim=2)
Phi_true = (-log(W))**theta

# # Sampling X and Y

n = 750 # size of sample

alpha = 2 # Weibull distribution shape parameter
beta = 2 # Weibull distribution scale parameter

MyCopula <- mvdc(copula=cop, # copula for (F(X), F(Y))
                 margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                 paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale

#set.seed(1)
XY <- rMvdc(n=n,MyCopula)
X <- XY[,1]
Y <- XY[,2]

# # Pseudo-sample

V = rep(NA,times=n)

for(i in 1:n){
  V[i] = (sum((X <= X[i])*(Y <= Y[i])))/(n-1)
}

# # Empirical distribution

Kn = ecdf(V)
Vn = V[order(V)]

# # # Estimator 1: differential equation

Phi_est = rep(0,times=n)

for(i in 1:n){
  if(Vn[i] < 1-1/n){
    T = unique(c(Vn[(Vn[i] <= Vn) & (Vn <= 1-1/n)],1-1/n))
    ProdTerm = rep(NA,times=length(T)-1)
    ProdTerm = (Kn(T[-length(T)])-T[-length(T)])/(Kn(T[-length(T)])-T[-1])
    Phi_est[i] = (1/n)*abs(prod(ProdTerm))
  }
}

Phi_est_diff = approx(x=c(0,Vn,1), y=c(Phi_est[1],Phi_est,0), xout=W, method="constant", ties = "ordered")$y 


# # # Estimator 2: integral equation

h = rep(NA,times=n)
h[n] = 1 

for(i in (n-1):1){
  h[i] = (1/(n*Kn(i/n)-i))*sum(h[(i+1):n])
}

# Step interpolation

cond = h>=0 & is.finite(h) # may have 0's, when Kn(i/n)-i/n <= 0

h.interp = approx(x=c(0,Vn[cond],1), y=c(h[cond][1],h[cond],1), xout=W, method="constant", ties = mean)$y 

Phi_est_int = (sum(h.interp) - cumsum(c(0,h.interp[-length(h.interp)])))/n


# # # Graph: comparison Phi_est_int, Phi_est_diff, Phi_true

Ylim = range(c(Phi_est_int,Phi_est_diff,Phi_true),na.rm = TRUE)

Xlabels = seq(from=0,to=1,by=0.2)
#Ylabels = format(seq(from=0,to=Ylim[2],length.out=5),scientific=TRUE,digits=2)
Ylabels = format(seq(from=0,to=Ylim[2],length.out=5),digits=2)#scientific=TRUE,

plot(W,Phi_true,type='l',col='blue',ylim=Ylim,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none")

lines(W,Phi_est_int,col='red',lwd = 2)
lines(W,Phi_est_diff,col='orange',lwd = 2)

axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "t", font=2,cex=1.5)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2,hadj=1,padj=0)

legend("topright",legend=c("true density","integral estimator","differential estimator"),lwd = 4,col=c("blue","red", "orange"),lty=1,cex=1.25,bty="n")


# # # Graph (for Statistics Papers article): comparision to empirical copula process

Kn_hat <- W+Phi_est_int/h.interp
Ylim <- range(c(Kn_hat,W),na.rm = TRUE)

plot(W,Kn(W),type="l",ylim=Ylim,xlab="t",ylab="Kn")
lines(W,Kn_hat,col="red")

# Write to csv (for LaTex)

write.csv(x=data.frame(cbind(W,Kn_hat,Kn(W))),
          file="C:/Users/djimd/OneDrive/Documents/Concordia - PhD/Thesis/Written Articles/Statistical Papers/clayton_empirical_kendall_750.csv",
          row.names = FALSE)

# # # # # # # # #  Limiting behavior of Estimator 1 (differential) # # # # # # # # # # # # # 


# # Simulating values from the copulas

SimuCopula <- function(CopulaName,W){
  
  # # Description
  # Function to simulate values of Archimedean copula and of the generator. 
  # Limited to the independence, Clayton, Frank and AMH copulas
  
  # # Inputs
  # CopulaName: name of the copula. One of 'Indep','Clayton','AMH' or 'Frank'
  # W: vector of values in (0,1)
  
  # # Outputs
  # A list with two elements: 
  # * cop: copula 
  # * gen_true: generator (phi)
  
  if(CopulaName == 'Indep'){
    # Independence copula
    
    cop <- indepCopula(dim=2)
    Phi_true = -log(W)
  }
  
  if(CopulaName == 'Clayton'){
    # Clayton copula
    
    theta = 3
    cop <- claytonCopula(param=theta,dim=2)
    Phi_true = (W**(-theta)-1)/theta
  }
  
  if(CopulaName == 'AMH'){
    # AMH copula
    
    theta = 0.5
    cop <- amhCopula(param=theta,dim=2)
    Phi_true = log((1-theta*(1-W))/W)/(1-theta)
  }
  
  if(CopulaName == 'Frank'){
    # Frank copula
    
    theta = -1
    cop <- frankCopula(param=theta,dim=2)
    Phi_true = -((exp(theta)-1)/theta)*log((exp(-theta*W)-1)/(exp(-theta)-1))
  }
  return(list(cop=cop,gen_true=Phi_true))
}

# Grid for the x-axis
W = seq(from=0.2,to=0.9,by=0.2) 

# CopulaNames <- c('Indep','Clayton','AMH','Frank')

CopulaName <- 'Frank' # select the copula
cop = SimuCopula(CopulaName,W)$cop
Phi_true = SimuCopula(CopulaName,W)$gen_true

# # Sampling X and Y

N <- seq(from=0,to=1000,by=50) # the 0 will be ignored. N should have at least 2 elements.

ConvergenceMAtrix <- matrix(rep(0,times=length(W)*(length(N)-1)),ncol=length(W)) # colnames = w, rownames = n
colnames(ConvergenceMAtrix) <- paste('w',1:length(W),sep='')

for(n.index in 2:length(N)){ # ignoring the value 0 value in N
  
  #n.index <- 1
  n <- N[n.index]
  
  alpha = 2 # Weibull distribution shape parameter
  beta = 2 # Weibull distribution scale parameter
  
  MyCopula <- mvdc(copula=cop, # copula for (F(X), F(Y))
                   margins=c("weibull","weibull"), # Weibull distribution for margins X and Y
                   paramMargins=list(shape=alpha,scale=beta)) # alpha:shape, beta:scale
  
  
  set.seed(5) # 5 for Clayton, AMH, indep
  XY <- rMvdc(n=n,MyCopula)
  X <- XY[,1]
  Y <- XY[,2]
  
  # # Pseudo-sample
  
  V = rep(NA,times=n)
  
  for(i in 1:n){
    V[i] = (sum((X <= X[i])*(Y <= Y[i])))/(n-1)
  }
  
  # # Empirical distribution
  
  Kn = ecdf(V)
  Vn = V[order(V)]
  
  # # # Estimator 1: differential equation
  
  Phi_est = rep(0,times=n)
  
  for(i in 1:n){
    if(Vn[i] < 1-1/n){
      T = unique(c(Vn[(Vn[i] <= Vn) & (Vn <= 1-1/n)],1-1/n))
      ProdTerm = rep(NA,times=length(T)-1)
      ProdTerm = (Kn(T[-length(T)])-T[-length(T)])/(Kn(T[-length(T)])-T[-1])
      Phi_est[i] = (1/n)*abs(prod(ProdTerm))
    }
  }
  
  Phi_est_diff = approx(x=c(0,Vn,1), y=c(Phi_est[1],Phi_est,0), xout=W, method="constant", ties = "ordered")$y 
  ConvergenceMAtrix[n.index-1,] <- (sqrt(n)/log(1/n, base = exp(1)))*(Phi_est_diff-Phi_true)/Phi_true
}

# ConvergenceMAtrix may contain infinite values. Replacing such values.
for(w.index in 1:length(W)){
  if(sum(is.infinite(ConvergenceMAtrix[,w.index])) != 0){
    ConvergenceMAtrix[,w.index][is.infinite(ConvergenceMAtrix[,w.index])] <- min(ConvergenceMAtrix[,w.index][is.finite(ConvergenceMAtrix[,w.index])])-0.25
  }
}

#apply(X=ConvergenceMAtrix,
#      MARGIN=2,
#      FUN=function(vecteur){min(vecteur[is.finite(vecteur)])-0.25})

# colors to be used to graph
MyColors <- c('aquamarine4','bisque1','blue','blueviolet','brown1','burlywood3','darkgreen')

# axes ranges
#Xlim = c(0,500)
Ylim = range(ConvergenceMAtrix,na.rm = TRUE)

# axis labels
Xlabels = N[-1] #seq(from=0,to=500,by=50)
Ylabels = format(seq(from=Ylim[1],to=Ylim[2],length.out=5),digits=2)

plot(N[-1],ConvergenceMAtrix[,1],type='l',col=MyColors[1],ylim=Ylim,lwd = 2,xlab="",ylab="",xaxt="none",yaxt="none")

for(w.index in 2:length(W)){
  lines(N[-1],ConvergenceMAtrix[,w.index],col=MyColors[w.index],lwd = 2)
}

axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=3, "sample size n", font=2,cex=1.25)#
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2,hadj=1,padj=0)

legend("bottomright",legend=paste('w=',W,sep=''),lwd = 3,col=MyColors[1:length(W)],lty=1,cex=1.25,bty="n")







