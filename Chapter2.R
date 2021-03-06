# # # # CHAPTER 2: ESTIMATION OF THE GENERATOR


# # Preps

library(copula) # for claytonCopula

m = 1000

W = 1:m/m #seq(from=1/n,to=1,by=0.01) # grid for the x-axis

# # Simulating values from Archimedean copulas

# Independence copula

cop <- indepCopula(dim=2)
Phi_true = -log(W)

# Clayton copula

theta = 3
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

n = 500 # size of sample

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
    #T = unique(c(W[i],Vn[which((W[i] <= Vn) & (Vn <= 1-1/n))],1-1/n))
    #T = unique(c(Vn[which((Vn[i] <= Vn) & (Vn <= 1-1/n))],1-1/n))
    T = unique(c(Vn[(Vn[i] <= Vn) & (Vn <= 1-1/n)],1-1/n))
    ProdTerm = rep(NA,times=length(T)-1)
    #for(j in 1:(length(T)-1)){
      #if((Kn(T[j])!=T[j])){# & (Kn(T[j])!=T[j+1])
      #  ProdTerm[j] = (Kn(T[j])-T[j])/(Kn(T[j])-T[j+1])
      #}else{
      #  ProdTerm[j] = 1
      #}
    #}
    ProdTerm = (Kn(T[-length(T)])-T[-length(T)])/(Kn(T[-length(T)])-T[-1])
    #ProdTerm = ProdTerm[is.finite(ProdTerm) & (ProdTerm != 0)]
    Phi_est[i] = (1/n)*abs(prod(ProdTerm))
  #  }else{
  #    Phi_est[i] = 0
  }
}

#Phi_est_diff = approx(x=c(0,Vn,1), y=c(Phi_est[1],Phi_est,0), xout=W, method="constant", ties = "ordered")$y 
Phi_est_diff = approx(x=c(0,Vn,1), y=c(Phi_est[1],Phi_est,0), xout=W, method="constant", ties = "ordered")$y 



# # # Estimator 2: integral equation

h = rep(NA,times=n)
h[n] = 1 

for(i in (n-1):1){
#  if(Kn(i/n)-i/n > 0){
  h[i] = (1/(n*Kn(i/n)-i))*sum(h[(i+1):n])
#  }else{
#    h[i] = h[i+1]
#  }
}

# Step interpolation
#Phi_est= rep(NA,times=length(W))

cond = h>=0 & is.finite(h) # may have 0's, when Kn(i/n)-i/n <= 0

h.interp = approx(x=c(0,Vn[cond],1), y=c(h[cond][1],h[cond],1), xout=W, method="constant", ties = mean)$y 

Phi_est_int = (sum(h.interp) - cumsum(c(0,h.interp[-n])))/n


# # # Graph

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



#Ylim = range(c(Phi_est_int,Phi_true),na.rm = TRUE)

#plot(W,Phi_est,type='l',col='red',ylim=Ylim)
#lines(W,Phi_true,col='blue')