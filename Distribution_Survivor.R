# # Kendall function for survival distribution

install.packages("GoFKernel")
library(GoFKernel)

#inverse = function (f, lower, upper) {
#   function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
#}

A = seq(from=0.01,to=1/2,by=0.01)
B = seq(from=0,to=500,by=0.001)
Z = seq(from=0.01, to=0.99,by=0.001)
M = 10**15
N = 10**4


# # Independence copula

# Known distribution of independence survival distribution

k_true  = -log(Z)

# Estimation through our method

inv_a_z = integrand = diff_a_z = matrix(NA,ncol=length(Z),nrow=length(A)) 

for(i in 1:length(A)){
  psi_inverse = inverse(function (z) {1-exp(-A[i]*z)-exp(-(1-A[i])*z)+exp(-z)}, 0, M)
  for(j in 1:length(Z)){
    inv_a_z[i,j] = unlist(psi_inverse(Z[j]))
    diff_a_z[i,j] = 1/(A[i]*exp(-A[i]*inv_a_z[i,j])+(1-A[i])*exp(-(1-A[i])*inv_a_z[i,j])-exp(-inv_a_z[i,j]))
    integrand[i,j] = inv_a_z[i,j] * exp(-inv_a_z[i,j]) * diff_a_z[i,j]
  }
}
k_est = rep(NA,length(Z)) #
for(j in 1:length(Z)){
  k_est[j] = 2*integrate(approxfun(x=A,y=integrand[,j]), lower = min(A), upper = max(A))[[1]]
}


# Graphical comparison for the densities

Ylim = range(c(k_true,k_est))
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(c(k_true,k_est))),to=ceiling(max(c(k_true,k_est))))
par(mfrow=c(1,1))
plot(Z,k_est,type="l",ylim = Ylim,lwd = 4,col="red",xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
lines(Z,k_true,col="green",lwd = 4)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2, "density", font=2,cex=1.5)
legend(x=0.6,y=4,legend=c("estimated density", "true density"),lwd = 4,col=c("red", "green"),lty=1,cex=1.25,bty="n")

# Adding simulation through McNeil and Neslehova method

#W = runif(n=N)
R = -log(runif(n=N))
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-exp(-R*Y1/(Y1+Y2))-exp(-R*Y2/(Y1+Y2))+(R*Y1/(Y1+Y2))*(R*Y2/(Y1+Y2))

Ylim = range(c(k_true,k_est,density(C_star,from=0,to=1,n=length(Z))$y))
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(c(k_true,k_est,density(C_star,from=0,to=1,n=length(Z))$y))),to=ceiling(max(c(k_true,k_est,density(C_star,from=0,to=1,n=length(Z))$y))))
par(mfrow=c(1,1))
plot(Z,k_est,type="l",ylim = Ylim,lwd = 4,col="red",xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
lines(Z,k_true,col="green",lwd = 4)
lines(density(C_star,from=0,to=1,n=length(Z)),col="blue",lwd = 4)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2, "density", font=2,cex=1.5)
legend(x=0.5,y=10,legend=c("our estimation", "true density","from simulated values"),lwd = 4,col=c("red", "green","blue"),lty=1,cex=1.25,bty="n")


# Density of the Kendall distribution

Ylim = range(k_true)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(k_true)),to=ceiling(max(k_true)))
par(mfrow=c(1,1))
plot(Z,k_true,type="l",ylim = Ylim,lwd = 4,col="green",xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2, "density", font=2,cex=1.5)

# Graphics for psi and its inverse

psi = function(a,b){1-exp(-a*b)-exp(-(1-a)*b)+exp(-b)}
Xlabels = seq(from=0,to=50,by=5)
Ylabels = seq(from=0,to=1,by=0.2)
par(mfrow=c(1,1))
plot(B[B<=50],psi(A[11],B[B<=50]),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)
plot(B[B<=50],psi(A[25],B[B<=50]),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)
plot(B[B<=50],psi(A[45],B[B<=50]),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)
Ylim = range(c(inv_a_z[11,],inv_a_z[25,],inv_a_z[45,]))
par(mfrow=c(1,1))
plot(Z,inv_a_z[11,],type="l",col="red",lwd = 4,xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
plot(Z,inv_a_z[25,],type="l",col="red",lwd = 4,xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
plot(Z,inv_a_z[45,],type="l",col="red",lwd = 4,xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)

# # # Clayton copula

Theta = c(0.1,1,3,5)
inv_a_z = integrand = diff_a_z = matrix(NA,ncol=length(Z),nrow=length(A)) 
k_est = matrix(NA,nrow=length(Theta),ncol=length(Z))
#i0 = j0 = 0
for(t in 1:length(Theta)){
  for(i in 1:length(A)){
    #i0 = i0 + 1
    psi_inverse = inverse(function (z) {1-(Theta[t]*A[i]*z+1)**(-1/Theta[t])-(Theta[t]*(1-A[i])*z+1)**(-1/Theta[t])+(Theta[t]*z+1)**(-1/Theta[t])}, 0, M)
    for(j in 1:length(Z)){
      #j0 = j0 + 1
      inv_a_z[i,j] = unlist(psi_inverse(Z[j]))
      diff_a_z[i,j] = 1/(A[i]*(Theta[t]*A[i]*inv_a_z[i,j]+1)**(-1/Theta[t]-1)+(1-A[i])*(Theta[t]*(1-A[i])*inv_a_z[i,j]+1)**(-1/Theta[t]-1)-(Theta[t]*inv_a_z[i,j]+1)**(-1/Theta[t]-1))
      integrand[i,j] = (Theta[t]+1)*inv_a_z[i,j]*(Theta[t]*inv_a_z[i,j]+1)**(-2-1/Theta[t]) * diff_a_z[i,j]
    }
  }
  
  for(j in 1:length(Z)){
    k_est[t,j] = 2*integrate(approxfun(x=A,y=integrand[,j]), lower = min(A), upper = max(A))[[1]]
  }
}

# Graphics for the estimated density (for survival Kendall function)

Ylim = range(k_est)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=0,to=4)
par(mfrow=c(1,1))

plot(Z,k_est[1,],type="l",lwd = 4,col="green",xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
lines(Z,k_est[2,],lwd = 4,col="blue")
lines(Z,k_est[3,],lwd = 4,col="red")
lines(Z,k_est[4,],lwd = 4,col="grey")
lines(Z,k_est[5,],lwd = 4,col="black")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.8,y=3.5,legend=expression(paste(theta,"=0.1"),paste(theta,"=1"),paste(theta,"=3"),paste(theta,"=5")),lwd=4,col=c("green","blue","red","grey","black"),lty=1,bty="n",cex=1.25)

# Comparison with simulations through McNeil and Neslehova method

t=2

W = runif(n=N)
R = (W**(-Theta[t])-1)/Theta[t]
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-(Theta[t]*R*Y1/(Y1+Y2)+1)**(-1/Theta[t])-(Theta[t]*R*Y2/(Y1+Y2)+1)**(-1/Theta[t])+W

Ylim = range(k_est,density(C_star,from=0,to=1,n=length(Z))$y)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=0,to=5)
par(mfrow=c(1,1))

plot(Z,k_est[t,],type="l",lwd = 4,col="red",xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
lines(density(C_star,from=0,to=1,n=length(Z)),lwd = 4,col="blue")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.5,y=4.5,legend=c("our estimation","from simulated values"),lwd=4,col=c("red","blue"),lty=1,bty="n",cex=1.25)


# Graphics for the density of Kendall function

kendall_copula = matrix(NA,nrow=length(Theta),ncol=length(Z))
for(t in 1:length(Theta)){
  kendall_copula[t,] = (1+1/Theta[t])*(1-Z**Theta[t])
}

Ylim = range(kendall_copula)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(kendall_copula)),to=ceiling(max(kendall_copula)))
par(mfrow=c(1,1))

plot(Z,kendall_copula[1,],type = "l",xlab="",ylab="",lwd = 4,col="green",ylim=Ylim,xaxt="none",yaxt="none")
lines(Z,kendall_copula[2,],type = "l",lwd = 4,col="blue")
lines(Z,kendall_copula[3,],type = "l",lwd = 4,col="red")
lines(Z,kendall_copula[4,],type = "l",lwd = 4,col="grey")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.8,y=4,legend=expression(paste(theta,"=0.1"),paste(theta,"=1"),paste(theta,"=3"),paste(theta,"=5")),lwd=3,col=c("green","blue","red","grey"),lty=1,bty="n")

# Graphics for psi and its inverse

psi = function(a,b){1-(Theta[2]*a*b+1)**(-1/Theta[2])-(Theta[2]*(1-a)*b+1)**(-1/Theta[2])+(Theta[2]*b+1)**(-1/Theta[2])}

Xlabels = seq(from=0,to=500,by=50)
Ylabels = seq(from=0,to=1,by=0.2)
Xlim = range(Xlabels)
Ylim = range(Ylabels)

par(mfrow=c(1,1))
plot(B,psi(A[11],B),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none",xlim=Xlim,ylim=Ylim)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)

plot(B,psi(A[25],B),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none",xlim=Xlim,ylim=Ylim)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)

plot(B,psi(A[45],B),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none",xlim=Xlim,ylim=Ylim)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)

Xlim = c(0,1)
Ylim = range(c(inv_a_z[11,inv_a_z[11,]<=500],inv_a_z[25,inv_a_z[11,]<=500],inv_a_z[45,inv_a_z[11,]<=500]))

par(mfrow=c(1,1))
plot(Z[inv_a_z[11,]<=500],inv_a_z[11,inv_a_z[11,]<=500],type="l",col="red",lwd = 4,xlab="",ylab="",xlim=Xlim,ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)

plot(Z[inv_a_z[25,]<=500],inv_a_z[25,inv_a_z[25,]<=500],type="l",col="red",lwd = 4,xlab="",ylab="",xlim=Xlim,ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)

plot(Z[inv_a_z[45,]<=500],inv_a_z[45,inv_a_z[45,]<=500],type="l",col="red",lwd = 4,xlab="",ylab="",xlim=Xlim,ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)



# # # Gumbel copula


Theta = c(1,1.5,3,5,10)  

inv_a_z = integrand = diff_a_z = matrix(NA,ncol=length(Z),nrow=length(A)) 
k_est = matrix(NA,nrow=length(Theta),ncol=length(Z))

for(t in 1:length(Theta)){
  for(i in 1:length(A)){
    psi_inverse = inverse(function (z) {1-exp(-(A[i]*z)**(1/Theta[t]))-exp(-((1-A[i])*z)**(1/Theta[t]))+exp(-z**(1/Theta[t]))}, 0, M)
    for(j in 1:length(Z)){
      inv_a_z[i,j] = unlist(psi_inverse(Z[j]))
      diff_a_z[i,j] = 1/((1/Theta[t])*A[i]**(1/Theta[t])*inv_a_z[i,j]**(1/Theta[t]-1)*exp(-(A[i]*inv_a_z[i,j])**(1/Theta[t]))+(1/Theta[t])*(1-A[i])**(1/Theta[t])*inv_a_z[i,j]**(1/Theta[t]-1)*exp(-((1-A[i])*inv_a_z[i,j])**(1/Theta[t]))-(1/Theta[t])*inv_a_z[i,j]**(1/Theta[t]-1)*exp(-inv_a_z[i,j]**(1/Theta[t])))
      integrand[i,j] = Theta[t]**(-2)*inv_a_z[i,j]**(1/Theta[t]-1)*exp(-inv_a_z[i,j]**(1/Theta[t]))*(inv_a_z[i,j]**(1/Theta[t])+Theta[t]-1) * diff_a_z[i,j]
    }
  }
  
  for(j in 1:length(Z)){
    if(sum(!is.nan(integrand[,j]))==length(A)){
      # make sure the vector integrand[,j] has at least 2 values, for approxfun
      k_est[t,j] = 2*integrate(approxfun(x=A,y=integrand[,j]), lower = min(A), upper = max(A))[[1]]
    }
  }
}

# Graphics for the estimated density of the survival copula

Ylim = range(k_est,na.rm = TRUE)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(k_est,na.rm = TRUE)),to=ceiling(max(k_est,na.rm = TRUE)))
par(mfrow=c(1,1))

plot(Z,k_est[1,],type="l",lwd = 3,col="green",xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
lines(Z,k_est[2,],type="l",lwd = 3,col="blue")
lines(Z,k_est[3,],type="l",lwd = 3,col="red")
lines(Z,k_est[4,],type="l",lwd = 3,col="grey")
lines(Z,k_est[5,],type="l",lwd = 3,col="black")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.6,y=4,legend=expression(paste(theta,"=1"),paste(theta,"=1.5"),paste(theta,"=3"),paste(theta,"=5"),paste(theta,"=10")),lwd=4,col=c("green","blue","red","grey","black"),lty=1,bty="n",cex=1.25)

# Comparison with simulations through McNeil and Neslehova method

t=2
W = runif(n=N)
R = (-log(W))**Theta[t]
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-exp(-(R*Y1/(Y1+Y2))**(1/Theta[t]))-exp(-(R*Y2/(Y1+Y2))**(1/Theta[t]))+W

Ylim = range(k_est,density(C_star,from=0,to=1,n=length(Z))$y,na.rm = TRUE)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(k_est,density(C_star,from=0,to=1,n=length(Z))$y,na.rm = TRUE)),to=ceiling(max(k_est,density(C_star,from=0,to=1,n=length(Z))$y,na.rm = TRUE)))
par(mfrow=c(1,1))

plot(Z,k_est[t,],type="l",lwd = 3,col="red",xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
lines(density(C_star,from=0,to=1,n=length(Z)),lwd = 3,col="blue")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.5,y=4,legend=c("our estimation","from simulated values"),lwd=4,col=c("red","blue"),lty=1,bty="n",cex=1.25)

# Graphics for the Kendall density (of the copula)

kendall_copula = matrix(NA,nrow=length(Theta),ncol=length(Z))
for(t in 1:length(Theta)){
  kendall_copula[t,] = (1/Theta[t])*(-log(Z) + Theta[t]-1)
}
Ylim = range(kendall_copula)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(kendall_copula)),to=ceiling(max(kendall_copula)))
par(mfrow=c(1,1))
plot(Z,kendall_copula[1,],type = "l",xlab="",ylab="",lwd = 4,col="green",ylim=Ylim,xaxt="none",yaxt="none")
lines(Z,kendall_copula[2,],type = "l",lwd = 4,col="blue")
lines(Z,kendall_copula[3,],type = "l",lwd = 4,col="red")
lines(Z,kendall_copula[4,],type = "l",lwd = 4,col="grey")
lines(Z,kendall_copula[5,],type = "l",lwd = 4,col="black")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.6,y=4,legend=expression(paste(theta,"=1"),paste(theta,"=1.5"),paste(theta,"=3"),paste(theta,"=5"),paste(theta,"=10")),lwd=4,col=c("green","blue","red","grey","black"),lty=1,bty="n")
# Graphics for psi and its inverse
psi = function(a,b){1-exp(-(a*b)**(1/Theta[3]))-exp(-((1-a)*b)**(1/Theta[3]))+exp(-b**(1/Theta[3]))}

Xlabels = seq(from=0,to=500,by=50)
Ylabels = seq(from=0,to=1,by=0.2)
Xlim = range(Xlabels)
Ylim = range(Ylabels)

par(mfrow=c(1,1))
plot(B,psi(A[11],B),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none",xlim=Xlim,ylim=Ylim)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)

plot(B,psi(A[25],B),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none",xlim=Xlim,ylim=Ylim)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)

plot(B,psi(A[45],B),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none",xlim=Xlim,ylim=Ylim)
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)

Xlim = c(0,1)
Ylim = range(c(inv_a_z[11,inv_a_z[11,]<=500],inv_a_z[25,inv_a_z[11,]<=500],inv_a_z[45,inv_a_z[11,]<=500]))

par(mfrow=c(1,1))
plot(Z[inv_a_z[11,]<=500],inv_a_z[11,inv_a_z[11,]<=500],type="l",col="red",lwd = 4,xlab="",ylab="",xlim=Xlim,ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)

plot(Z[inv_a_z[25,]<=500],inv_a_z[25,inv_a_z[25,]<=500],type="l",col="red",lwd = 4,xlab="",ylab="",xlim=Xlim,ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)

plot(Z[inv_a_z[45,]<=500],inv_a_z[45,inv_a_z[45,]<=500],type="l",col="red",lwd = 4,xlab="",ylab="",xlim=Xlim,ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Ylabels,labels=Ylabels,las=1,font=2)
axis(2, at=Xlabels,labels=Xlabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)


# # # Ali-Mikhail-Haq 

Theta = c(-1,-1/2,0,1/2,0.99)  # avoid value 1

inv_a_z = integrand = diff_a_z = matrix(NA,ncol=length(Z),nrow=length(A)) 
k_est = matrix(NA,nrow=length(Theta),ncol=length(Z))

for(t in 1:length(Theta)){
  for(i in 1:length(A)){
    psi_inverse = inverse(function (z) {1-(1-Theta[t])/(exp(A[i]*z)-Theta[t])-(1-Theta[t])/(exp((1-A[i])*z)-Theta[t])+(1-Theta[t])/(exp(z)-Theta[t])}, 0, M)
    for(j in 1:length(Z)){
      inv_a_z[i,j] = unlist(psi_inverse(Z[j]))
      diff_a_z[i,j] = 1/((1-Theta[t])*A[i]*exp(A[i]*inv_a_z[i,j])/((exp(A[i]*inv_a_z[i,j])-Theta[t])**2) + (1-Theta[t])*(1-A[i])*exp((1-A[i])*inv_a_z[i,j])/((exp((1-A[i])*inv_a_z[i,j])-Theta[t])**2)-(1-Theta[t])*exp(inv_a_z[i,j])/((exp(inv_a_z[i,j])-Theta[t])**2))
      integrand[i,j] = ((1-Theta[t])*inv_a_z[i,j]*(exp(2*inv_a_z[i,j])-Theta[t])*exp(inv_a_z[i,j]))/(exp(inv_a_z[i,j])-Theta[t])**4 * diff_a_z[i,j]
    }
  }
  
  for(j in 1:length(Z)){
    if(sum(!is.nan(integrand[,j]) & is.finite(integrand[,j]))==length(A)){
      # make sure the vector integrand[,j] has at least 2 finite values, for approxfun
      k_est[t,j] = 2*integrate(approxfun(x=A,y=integrand[,j]), lower = min(A), upper = max(A))[[1]]
    }
  }
}

# Graphics for the estimated density of the survival density

Ylim = range(k_est,na.rm = TRUE)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(k_est,na.rm = TRUE)),to=ceiling(max(k_est,na.rm = TRUE)))
par(mfrow=c(1,1))

plot(Z,k_est[1,],type="l",lwd = 4,col="red",xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
lines(Z,k_est[2,],type="l",lwd = 4,col="blue")
lines(Z,k_est[3,],type="l",lwd = 4,col="green")
lines(Z,k_est[4,],type="l",lwd = 4,col="grey")
lines(Z,k_est[5,],type="l",lwd = 4,col="black")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.6,y=7,legend=expression(paste(theta,"=-1"),paste(theta,"=-0.5"),paste(theta,"=0"),paste(theta,"=0.5"),paste(theta,"=0.99")),lwd=4,col=c("red","blue","green","grey","black"),lty=1,bty="n",cex=1.25)

# Comparison with simulations through McNeil and Neslehova method

t=2
N = 10000
W = runif(n=N)
R = log((1-Theta[t]*(1-W))/W)
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-(1-Theta[t])/(exp(R*Y1/(Y1+Y2))-Theta[t])-(1-Theta[t])/(exp(R*Y2/(Y1+Y2))-Theta[t])+W

Ylim = range(k_est,density(C_star,from=0,to=1,n=length(Z))$y,na.rm = TRUE)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(k_est,density(C_star,from=0,to=1,n=length(Z))$y,na.rm = TRUE)),to=ceiling(max(k_est,density(C_star,from=0,to=1,n=length(Z))$y,na.rm = TRUE)),by=2)
par(mfrow=c(1,1))

plot(Z,k_est[t,],type="l",lwd = 4,col="red",xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
lines(density(C_star,from=0,to=1,n=length(Z)),type="l",lwd = 4,col="blue")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.5,y=12,legend=c("our estimation","from simulated values"),lwd=4,col=c("red","blue","green","grey","black"),lty=1,bty="n",cex=1.25)

# Graphics for the Kendall density (of the copula)

kendall_copula = matrix(NA,nrow=length(Theta),ncol=length(Z))
for(t in 1:length(Theta)){
  kendall_copula[t,] = (1-Theta[t]+Theta[t]*Z*(2-Z))*log((1-Theta[t]*(1-Z))/Z)/(1-Theta[t])
}

Ylim = range(kendall_copula)
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=floor(min(kendall_copula)),to=ceiling(max(kendall_copula)))
par(mfrow=c(1,1))

plot(Z,kendall_copula[1,],type = "l",xlab="",ylab="",lwd = 4,col="red",ylim=Ylim,xaxt="none",yaxt="none")
lines(Z,kendall_copula[2,],type = "l",lwd = 4,col="blue")
lines(Z,kendall_copula[3,],type = "l",lwd = 4,col="green")
lines(Z,kendall_copula[4,],type = "l",lwd = 4,col="grey")
lines(Z,kendall_copula[5,],type = "l",lwd = 4,col="black")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
mtext(side=2, line=2.5, "density", font=2,cex=1.5)
legend(x=0.6,y=5,legend=expression(paste(theta,"=-1"),paste(theta,"=-0.5"),paste(theta,"=0"),paste(theta,"=0.5"),paste(theta,"=0.99")),lwd=4,col=c("red","blue","green","grey","black"),lty=1,bty="n")

# Graphics for psi and its inverse

psi = function(a,b){1-(1-Theta[4])/(exp(a*b)-Theta[4])-(1-Theta[4])/(exp((1-a)*b)-Theta[4])+(1-Theta[4])/(exp(b)-Theta[4])}
Xlabels = seq(from=0,to=50,by=5)
Ylabels = seq(from=0,to=1,by=0.2)
par(mfrow=c(1,1))
plot(B[B<=50],psi(A[11],B[B<=50]),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)
plot(B[B<=50],psi(A[25],B[B<=50]),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)
plot(B[B<=50],psi(A[45],B[B<=50]),type="l",col="blue",lwd = 4,xlab="",ylab="",xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2.25, "b", font=2,cex=1.5)
Ylim = range(c(inv_a_z[11,],inv_a_z[25,],inv_a_z[45,]))
Xlabels = seq(from=0,to=1,by=0.2)
Ylabels = seq(from=0,to=6,by=1)
par(mfrow=c(1,1))
plot(Z,inv_a_z[11,],type="l",col="red",lwd = 4,xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
plot(Z,inv_a_z[25,],type="l",col="red",lwd = 4,xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)
plot(Z,inv_a_z[45,],type="l",col="red",lwd = 4,xlab="",ylab="",ylim=Ylim,xaxt="none",yaxt="none")
axis(1, at=Xlabels,labels=Xlabels,las=1,font=2)
axis(2, at=Ylabels,labels=Ylabels,las=1,font=2)
mtext(side=1, line=2, "z", font=2,cex=1.5)


# # Simulation (McNeil and Neslehova): histogrammes

# Independence

N = 10000
W = runif(n=N)
R = -log(W)
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-exp(-R*Y1/(Y1+Y2))-exp(-R*Y2/(Y1+Y2))+W
hist(C_star,breaks=100)
lines(density(C_star))

# Clayton

Theta = c(0.1,1,3,5)
t = 2
N = 10000
W = runif(n=N)
R = (W**(-Theta[t])-1)/Theta[t]
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-(Theta[t]*R*Y1/(Y1+Y2)+1)**(-1/Theta[t])-(Theta[t]*R*Y2/(Y1+Y2)+1)**(-1/Theta[t])+W
hist(C_star,breaks=100)

# Gumbel

Theta = c(1,1.5,3,5,10)
t = 2
N = 10000
W = runif(n=N)
R = (-log(W))**Theta[t]
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-exp(-(R*Y1/(Y1+Y2))**(1/Theta[t]))-exp(-(R*Y2/(Y1+Y2))**(1/Theta[t]))+W
hist(C_star,breaks=100)

# Ali-Mikhail-Haq

Theta = c(-1,-1/2,0,1/2,0.99) 
t = 2
N = 10000
W = runif(n=N)
R = log((1-Theta[t]*(1-W))/W)
Y1 = rexp(n=N)
Y2 = rexp(n=N)
C_star = 1-(1-Theta[t])/(exp(R*Y1/(Y1+Y2))-Theta[t])-(1-Theta[t])/(exp(R*Y2/(Y1+Y2))-Theta[t])+W
hist(C_star,breaks=100)
