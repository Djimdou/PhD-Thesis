# Data from 
# McGilchrist, C. A., and C. W. Aisbett. "Regression with Frailty in Survival Analysis."
# Biometrics, vol. 47, no. 2, 1991, pp. 461-466.


# # #  Numerical application

# # Estimator 



# # Kendall's tau




# # # Real data

# install.packages("tabulizer")
#install.packages("rJava")
#library(rJava)
library(tabulizer)

#install.packages("dplyr")
library(dplyr)

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


# # Kendall's tau









# # Police shooting in US

wash_data = read.csv2(file="https://raw.githubusercontent.com/washingtonpost/data-police-shootings/master/fatal-police-shootings-data.csv",header = TRUE,sep = ",")

table(wash_data[,'race'])
year = lubridate::year(wash_data[,'date'])

table(lubridate::year(wash_data[,'date']),wash_data[,'race'])

armed = wash_data[wash_data[,'armed']=='gun',]

table(year = lubridate::year(armed[,'date']),armed[,'race'])
