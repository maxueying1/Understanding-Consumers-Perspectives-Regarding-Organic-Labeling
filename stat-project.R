setwd("D:/RA/STAT PROJECT/v3 this one")
data <- read.csv("data.csv")
data <- data[seq(1, nrow(data), by = 6), ]

data$label <- ifelse(data$label_nonGMO == 1 | data$label_local == 1 | data$label_healthy == 1|
                     data$label_domestic == 1 | data$label_farmer == 1 | data$label_eco == 1|
                     data$label_pesticide == 1 | data$label_sustainable == 1 | data$label_private == 1,1, 0)


data <- data[, c("label_organic", "age", "female","college", 
                 "chronic", "main_habits", "label")]

##stepAIC
library(MASS)
#model=glm(label_organic~.,data=data)
#step_model <- stepAIC(model, direction = "both")

model_1<-glm(label_organic ~ age + female + college + 
               + chronic + main_habits + label,data=data)
summary(model_1)

# generate predicted values for the model using the predict() function
predicted <- predict(model_1, type="response")
# convert predicted values to binary format
predicted_binary <- ifelse(predicted > 0.5, 1, 0)
# assume that the actual values are stored in a variable called "actual"
# calculate the accuracy score
accuracy <- sum(predicted_binary == data$label_organic) / length(data$label_organic)
# print the accuracy score
print(paste0("Accuracy score: ", accuracy))


library(car)
model<-glm(label_organic~.,data=data)
vif(model)

data=data.matrix(data)

c=ncol(data)   ##number of variables
n=nrow(data)  ##number of observations
summary(data)

##desperate data into repose variables and independent variables
y=data[,1]
X=data[,2:7]
x=cbind(1,X)
q=dim(x)[2]


## plot
par(mfrow=c(2,3))
plot(age,label_organic) 
plot(female,label_organic) 
plot(college,label_organic) 
plot(chronic,label_organic) 
plot(main_habits,label_organic) 
plot(label,label_organic) 
pairs(~label_organic+age+female+college+chronic+main_habits+label)



##create correlation plot
library(corrplot)
corr_matrix <- cor(X)
corr_matrix
corrplot(corr_matrix, method = "color",tl.cex=0.75)

#------------------------------------------------------------------------------#
#------------------------------ the logit model -------------------------------#
#------------------------------------------------------------------------------#
set.seed(12345)
invlogit=function(beta,x){
  pr=exp(x%*%beta)/(1+exp(x%*%beta))
  return(pr)
}
oldbeta=rep(0,q) ##c is # of columns
tol=1e-5  
err=1
maxits=3000
its=0
while(err>tol & its<maxits){
  old_pr=invlogit(oldbeta,x)
  s=t(x)%*%(y-old_pr)
  v=matrix(rep(0,n*n),n,n)
  diag(v)=old_pr*(1-old_pr)
  h=t(x)%*%v%*%x
  newbeta=oldbeta+(solve(h))%*%s
  err=max(abs(newbeta-oldbeta))
  its=its+1
  oldbeta=newbeta
}
newbeta
##estimate value using newbeta
temp=x%*%newbeta
pr=exp(temp)/(1+exp(temp))
yhat=rbinom(n,1,pr)
head(yhat)
##calculate the accuracy score
accuracy=function(y,yhat){
  mistake=0
  for (i in 1:n){
    if (yhat[i] !=y[i]){mistake=mistake+1}
  }
  score=(n-mistake)/n
  return(score)
}
score=function(y,beta_est,n,x){
  temp=x%*%beta_est
  pr=exp(temp)/(1+exp(temp))
  yhat=rbinom(n,1,pr)
  score=accuracy(y,yhat)
  return(score)
}
score=accuracy(y,yhat)
score

#####use only 13 variables to model
X2=X[,-1] ##delete the age variable from X
x2=cbind(1,X2)
q=dim(x2)[2]
########the logit model
oldbeta=rep(0,q)  ##c is # of columns
tol=1e-5  
err=1
maxits=3000
its=0
while(err>tol & its<maxits){
  old_pr=invlogit(oldbeta,x2)
  s=t(x2)%*%(y-old_pr)
  v=matrix(rep(0,n*n),n,n)
  diag(v)=old_pr*(1-old_pr)
  h=t(x2)%*%v%*%x2
  newbeta=oldbeta+(solve(h))%*%s
  err=max(abs(newbeta-oldbeta))
  its=its+1
  oldbeta=newbeta
}
newbeta
##estimate value using newbeta
temp=x2%*%newbeta
pr=exp(temp)/(1+exp(temp))
yhat=rbinom(n,1,pr)
head(yhat)
##calculate the accuracy score
score=accuracy(y,yhat)
score

#####use only 12 variables to model
X3=X2[,-4] ##delete the female variable from X
x3=cbind(1,X3)
q=dim(x3)[2]
########the logit model
oldbeta=rep(0,q)  ##c is # of columns
tol=1e-5  
err=1
maxits=3000
its=0
while(err>tol & its<maxits){
  old_pr=invlogit(oldbeta,x3)
  s=t(x3)%*%(y-old_pr)
  v=matrix(rep(0,n*n),n,n)
  diag(v)=old_pr*(1-old_pr)
  h=t(x3)%*%v%*%x3
  newbeta=oldbeta+(solve(h))%*%s
  err=max(abs(newbeta-oldbeta))
  its=its+1
  oldbeta=newbeta
}
newbeta
##estimate value using newbeta
temp=x3%*%newbeta
pr=exp(temp)/(1+exp(temp))
yhat=rbinom(n,1,pr)
head(yhat)
##calculate the accuracy score
score=accuracy(y,yhat)
score

######################################## 
## will use 5 independent variables


#------------------------------------------------------------------------------#
#--------------------------------- MCMC - MH ----------------------------------#
#------------------------------------------------------------------------------#
y=data[,1]
x=x2  ##x is n*6
N=2000

##initialize
beta0=c() ##this is the constant term
beta1=c()
beta2=c()
beta3=c()
beta4=c()
beta5=c()
beta0[1]=-1
beta1[1]=-1
beta2[1]=1
beta3[1]=-1
beta4[1]=1
beta5[1]=1

beta=c(beta0,beta1,beta2,beta3,beta4,beta5)  #beta is 6*1
ex=function(x, beta){
  y=exp(x%*%beta)
  return(y)
}
likelihood=function(beta){
  l=1
  for (i in 1:n){
    l=l*(((ex(x[i,], beta))/(1+ex(x[i,],beta)))^(y[i]))*((1-((ex(x[i,],beta))/(1+ex(x[i,],beta))))^(1-y[i]))*3.27
  }
  return(l)
}
##define the conditional distribution
cond=function(beta,betak){     
  likelihood(beta)*(exp(-(betak^2)/50))
}
qd1=function(x){
  y=exp(-((x+0)^2)/50)
  return(y)
}
qd2=function(x){
  y=exp(-((x-3)^2)/50)
  return(y)
}
qd3=function(x){
  y=exp(-((x+3)^2)/50)
  return(y)
}

#likelihood(beta)

###MH 
for (j in 2:N){
  ##the proposals
  z0=rnorm(1,-3,5)
  z1=rnorm(1,0,5)
  z2=rnorm(1,0,5)
  z3=rnorm(1,0,5)
  z4=rnorm(1,0,5)
  z5=rnorm(1,3,5)

  den0=cond(beta,beta0[j-1])*qd3(z0)
  beta_xi0=rbind(z0,beta[2],beta[3],beta[4],beta[5],beta[6])
  num0=cond(beta_xi0,z0)*qd3(beta0[j-1])
  rho0=num0/den0
  
  den1=cond(beta,beta1[j-1])*qd1(z1)
  beta_xi1=rbind(beta[1],z1,beta[3],beta[4],beta[5],beta[6])
  num1=cond(beta_xi1,z1)*qd1(beta1[j-1])
  rho1=num1/den1
  
  den2=cond(beta,beta2[j-1])*qd1(z2)
  beta_xi2=rbind(beta[1],beta[2],z2,beta[4],beta[5],beta[6])
  num2=cond(beta_xi2,z2)*qd1(beta2[j-1])
  rho2=num2/den2
  
  den3=cond(beta,beta3[j-1])*qd1(z3)
  beta_xi3=rbind(beta[1],beta[2],beta[3],z3,beta[5],beta[6])
  num3=cond(beta_xi3,z3)*qd1(beta3[j-1])
  rho3=num3/den3
  
  den4=cond(beta,beta4[j-1])*qd1(z4)
  beta_xi4=rbind(beta[1],beta[2],beta[3],beta[4],z4,beta[6])
  num4=cond(beta_xi4,z4)*qd1(beta4[j-1])
  rho4=num4/den4
  
  den5=cond(beta,beta5[j-1])*qd2(z5)
  beta_xi5=rbind(beta[1],beta[2],beta[3],beta[4],beta[5],z5)
  num5=cond(beta_xi5,z5)*qd2(beta5[j-1])
  rho5=num5/den5
  
  if(is.na(rho0)==T){rho0=1}
  if(is.na(rho1)==T){rho1=1}
  if(is.na(rho2)==T){rho2=1}
  if(is.na(rho3)==T){rho3=1}
  if(is.na(rho4)==T){rho4=1}
  if(is.na(rho5)==T){rho5=1}
  
  if(rho0>=1){beta0[j]=z0}else{coin0=rbinom(1,1,rho0)
  if(coin0==1){beta0[j]=z0}else{beta0[j]=beta0[j-1]}}
  
  if(rho1>=1){beta1[j]=z1}else{coin1=rbinom(1,1,rho1)
  if(coin1==1){beta1[j]=z1}else{beta1[j]=beta1[j-1]}}
  
  if(rho2>=1){beta2[j]=z2}else{coin2=rbinom(1,1,rho2)
  if(coin2==1){beta2[j]=z2}else{beta2[j]=beta2[j-1]}}
  
  if(rho3>=1){beta3[j]=z3}else{coin3=rbinom(1,1,rho3)
  if(coin3==1){beta3[j]=z3}else{beta3[j]=beta3[j-1]}}
  
  if(rho4>=1){beta4[j]=z4}else{coin4=rbinom(1,1,rho4)
  if(coin4==1){beta4[j]=z4}else{beta4[j]=beta4[j-1]}}
  
  if(rho5>=1){beta5[j]=z5}else{coin5=rbinom(1,1,rho5)
  if(coin5==1){beta5[j]=z5}else{beta5[j]=beta5[j-1]}}
  
  beta=c(beta0[j],beta1[j],beta2[j],beta3[j],beta4[j],beta5[j])
}

newbeta=cbind(beta0,beta1,beta2,beta3,beta4,beta5)
betab=colMeans(newbeta)
betab

# create a 2x3 grid of plots
par(mfrow = c(2, 3))

plot(beta0,type="l", main = "constant term")
plot(beta1,type="l", main = "female")
plot(beta2,type="l", main = "college")
plot(beta3,type="l", main = "chronic")
plot(beta4,type="l", main = "main_habits")
plot(beta5,type="l", main = "label")

newbeta_burn=newbeta[0:500,]
betab_burn=colMeans(newbeta_burn)
betab_burn
##estimate value using newbeta
temp=x%*%betab_burn
pr=exp(temp)/(1+exp(temp))
yhat=rbinom(n,1,pr)
head(yhat)
##calculate the accuracy score
score=accuracy(y,yhat)
score


#------------------------------------------------------------------------------#
#---------------- 5 times  ------ MCMC - MH ----------------------------------#
#------------------------------------------------------------------------------#

N=2000

mcmc=function(a,b,c,d,e,f){
  beta0=c()
  beta1=c()
  beta2=c()
  beta3=c()
  beta4=c()
  beta5=c()
  beta0[1]=a
  beta1[1]=b
  beta2[1]=c
  beta3[1]=d
  beta4[1]=e
  beta5[1]=f
  beta=c(beta0,beta1,beta2,beta3,beta4,beta5)  #beta is 6*1
  ex=function(x, beta){
    y=exp(x%*%beta)
    return(y)
  }
  likelihood=function(beta){
    l=1
    for (i in 1:n){
      l=l*(((ex(x[i,], beta))/(1+ex(x[i,],beta)))^(y[i]))*((1-((ex(x[i,],beta))/(1+ex(x[i,],beta))))^(1-y[i]))*3.27
    }
    return(l)
  }
  ##define the conditional distribution
  cond=function(beta,betak){     
    likelihood(beta)*(exp(-(betak^2)/50))
  }
  qd1=function(x){
    y=exp(-((x+1)^2)/50)
    return(y)
  }
  qd2=function(x){
    y=exp(-((x-1)^2)/50)
    return(y)
  }
  qd3=function(x){
    y=exp(-((x+10)^2)/50)
    return(y)
  }
  
  #likelihood(beta)
  
  ###MH 
  for (j in 2:N){
    ##the proposals
    z0=rnorm(1,-3,5)
    z1=rnorm(1,0,5)
    z2=rnorm(1,0,5)
    z3=rnorm(1,0,5)
    z4=rnorm(1,0,5)
    z5=rnorm(1,3,5)
    
    den0=cond(beta,beta0[j-1])*qd3(z0)
    beta_xi0=rbind(z0,beta[2],beta[3],beta[4],beta[5],beta[6])
    num0=cond(beta_xi0,z0)*qd3(beta0[j-1])
    rho0=num0/den0
    
    den1=cond(beta,beta1[j-1])*qd1(z1)
    beta_xi1=rbind(beta[1],z1,beta[3],beta[4],beta[5],beta[6])
    num1=cond(beta_xi1,z1)*qd1(beta1[j-1])
    rho1=num1/den1
    
    den2=cond(beta,beta2[j-1])*qd1(z2)
    beta_xi2=rbind(beta[1],beta[2],z2,beta[4],beta[5],beta[6])
    num2=cond(beta_xi2,z2)*qd1(beta2[j-1])
    rho2=num2/den2
    
    den3=cond(beta,beta3[j-1])*qd1(z3)
    beta_xi3=rbind(beta[1],beta[2],beta[3],z3,beta[5],beta[6])
    num3=cond(beta_xi3,z3)*qd1(beta3[j-1])
    rho3=num3/den3
    
    den4=cond(beta,beta4[j-1])*qd1(z4)
    beta_xi4=rbind(beta[1],beta[2],beta[3],beta[4],z4,beta[6])
    num4=cond(beta_xi4,z4)*qd1(beta4[j-1])
    rho4=num4/den4
    
    den5=cond(beta,beta5[j-1])*qd2(z5)
    beta_xi5=rbind(beta[1],beta[2],beta[3],beta[4],beta[5],z5)
    num5=cond(beta_xi5,z5)*qd2(beta5[j-1])
    rho5=num5/den5
    
    if(is.na(rho0)==T){rho0=1}
    if(is.na(rho1)==T){rho1=1}
    if(is.na(rho2)==T){rho2=1}
    if(is.na(rho3)==T){rho3=1}
    if(is.na(rho4)==T){rho4=1}
    if(is.na(rho5)==T){rho5=1}
    
    if(rho0>=1){beta0[j]=z0}else{coin0=rbinom(1,1,rho0)
    if(coin0==1){beta0[j]=z0}else{beta0[j]=beta0[j-1]}}
    
    if(rho1>=1){beta1[j]=z1}else{coin1=rbinom(1,1,rho1)
    if(coin1==1){beta1[j]=z1}else{beta1[j]=beta1[j-1]}}
    
    if(rho2>=1){beta2[j]=z2}else{coin2=rbinom(1,1,rho2)
    if(coin2==1){beta2[j]=z2}else{beta2[j]=beta2[j-1]}}
    
    if(rho3>=1){beta3[j]=z3}else{coin3=rbinom(1,1,rho3)
    if(coin3==1){beta3[j]=z3}else{beta3[j]=beta3[j-1]}}
    
    if(rho4>=1){beta4[j]=z4}else{coin4=rbinom(1,1,rho4)
    if(coin4==1){beta4[j]=z4}else{beta4[j]=beta4[j-1]}}
    
    if(rho5>=1){beta5[j]=z5}else{coin5=rbinom(1,1,rho5)
    if(coin5==1){beta5[j]=z5}else{beta5[j]=beta5[j-1]}}
    
    beta=c(beta0[j],beta1[j],beta2[j],beta3[j],beta4[j],beta5[j])
  }
  newbeta=cbind(beta0,beta1,beta2,beta3,beta4,beta5)
  return(newbeta)
}

mcmc1=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc2=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc3=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc4=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc5=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc6=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))

beta_mcmc1=colMeans(mcmc1)
beta_mcmc2=colMeans(mcmc2)
beta_mcmc3=colMeans(mcmc3)
beta_mcmc4=colMeans(mcmc4)
beta_mcmc5=colMeans(mcmc5)

beta_mcmc=colMeans(rbind(beta_mcmc1,beta_mcmc2,beta_mcmc3,beta_mcmc4,beta_mcmc5))
beta_mcmc

newbeta_burn=newbeta[0:500,]
betab_burn=colMeans(newbeta_burn)
betab_burn

score=function(y,beta_est,n,x){
  temp=x%*%beta_est
  pr=exp(temp)/(1+exp(temp))
  yhat=rbinom(n,1,pr)
  score=accuracy(y,yhat)
  return(score)
}
score(y,beta_mcmc,n,x)


####burn no difference on score
beta_mcmc1=colMeans(mcmc1[500:N,])
beta_mcmc2=colMeans(mcmc2[500:N,])
beta_mcmc3=colMeans(mcmc3[500:N,])
beta_mcmc4=colMeans(mcmc4[500:N,])
beta_mcmc5=colMeans(mcmc5[500:N,])

beta_mcmc=colMeans(rbind(beta_mcmc1,beta_mcmc2,beta_mcmc3,beta_mcmc4,beta_mcmc5))


par(mfrow = c(2, 3))

plot(beta0,type="l", main = "constant term")
plot(beta1,type="l", main = "female")
plot(beta2,type="l", main = "college")
plot(beta3,type="l", main = "chronic")
plot(beta4,type="l", main = "main_habits")
plot(beta5,type="l", main = "label")

##beta0
plot(mcmc1[,1],type="l",col="red", main = "constant term")
lines(mcmc2[,1],col="blue")
lines(mcmc3[,1],col="green")
lines(mcmc4[,1],col="black")
lines(mcmc5[,1],col="yellow")


##beta1
plot(mcmc1[,2],type="l",col="red", main = "female")
lines(mcmc2[,2],col="blue")
lines(mcmc3[,2],col="green")
lines(mcmc4[,2],col="black")
lines(mcmc5[,2],col="yellow")


##beta2
plot(mcmc1[,3],type="l",col="red", main = "college")
lines(mcmc2[,3],col="blue")
lines(mcmc3[,3],col="green")
lines(mcmc4[,3],col="black")
lines(mcmc5[,3],col="yellow")


##beta3
plot(mcmc1[,4],type="l",col="red", main = "chronic")
lines(mcmc2[,4],col="blue")
lines(mcmc3[,4],col="green")
lines(mcmc4[,4],col="black")
lines(mcmc5[,4],col="yellow")


##beta4
plot(mcmc1[,5],type="l",col="red", main = "main_habits")
lines(mcmc2[,5],col="blue")
lines(mcmc3[,5],col="green")
lines(mcmc4[,5],col="black")
lines(mcmc5[,5],col="yellow")


##beta5
plot(mcmc1[,6],type="l",col="red", main = "label")
lines(mcmc2[,6],col="blue")
lines(mcmc3[,6],col="green")
lines(mcmc4[,6],col="black")
lines(mcmc5[,6],col="yellow")


#------------------------------------------------------------------------------#
#---------------- 5 times random proposals ------ MCMC - MH -------------------#
#------------------------------------------------------------------------------#
N=2000
var0 = runif(1,1,10)
var1 = runif(1,1,10)
var2 = runif(1,1,10)
var3 = runif(1,1,10)
var4 = runif(1,1,10)
var5 = runif(1,1,10)

mcmc=function(a,b,c,d,e,f){
  beta0=c()
  beta1=c()
  beta2=c()
  beta3=c()
  beta4=c()
  beta5=c()
  beta0[1]=a
  beta1[1]=b
  beta2[1]=c
  beta3[1]=d
  beta4[1]=e
  beta5[1]=f
  beta=c(beta0,beta1,beta2,beta3,beta4,beta5)  #beta is 6*1
  ex=function(x, beta){
    y=exp(x%*%beta)
    return(y)
  }
  likelihood=function(beta){
    l=1
    for (i in 1:n){
      l=l*(((ex(x[i,], beta))/(1+ex(x[i,],beta)))^(y[i]))*((1-((ex(x[i,],beta))/(1+ex(x[i,],beta))))^(1-y[i]))*3.27
    }
    return(l)
  }
  ##define the conditional distribution
  cond=function(beta,betak,var){     
    likelihood(beta)*(exp(-(betak^2)/var^2))
  }
  qd=function(x,avg,var){
    y=exp(-((x-avg)^2)/var^2)
    return(y)
  }
  
  
  #likelihood(beta)
  
  ###MH 
  for (j in 2:N){
    ##the proposals
    z0=rnorm(1,-3,var0)
    z1=rnorm(1,0,var1)
    z2=rnorm(1,0,var2)
    z3=rnorm(1,0,var3)
    z4=rnorm(1,0,var4)
    z5=rnorm(1,3,var5)
    
    den0=cond(beta,beta0[j-1],var0)*qd(z0,-3,var0)
    beta_xi0=rbind(z0,beta[2],beta[3],beta[4],beta[5],beta[6])
    num0=cond(beta_xi0,z0,var0)*qd(beta0[j-1],-3,var0)
    rho0=num0/den0
    
    den1=cond(beta,beta1[j-1],var1)*qd(z1,0,var1)
    beta_xi1=rbind(beta[1],z1,beta[3],beta[4],beta[5],beta[6])
    num1=cond(beta_xi1,z1,var1)*qd(beta1[j-1],0,var1)
    rho1=num1/den1
    
    den2=cond(beta,beta2[j-1],var2)*qd(z2,0,var2)
    beta_xi2=rbind(beta[1],beta[2],z2,beta[4],beta[5],beta[6])
    num2=cond(beta_xi2,z2,var2)*qd(beta2[j-1],0,var2)
    rho2=num2/den2
    
    den3=cond(beta,beta3[j-1],var3)*qd(z3,0,var3)
    beta_xi3=rbind(beta[1],beta[2],beta[3],z3,beta[5],beta[6])
    num3=cond(beta_xi3,z3,var3)*qd(beta3[j-1],0,var3)
    rho3=num3/den3
    
    den4=cond(beta,beta4[j-1],var4)*qd(z4,0,var4)
    beta_xi4=rbind(beta[1],beta[2],beta[3],beta[4],z4,beta[6])
    num4=cond(beta_xi4,z4,var4)*qd(beta4[j-1],0,var4)
    rho4=num4/den4
    
    den5=cond(beta,beta5[j-1],var5)*qd(z5,3,var5)
    beta_xi5=rbind(beta[1],beta[2],beta[3],beta[4],beta[5],z5)
    num5=cond(beta_xi5,z5,var5)*qd(beta5[j-1],3,var5)
    rho5=num5/den5
    
    if(is.na(rho0)==T){rho0=1}
    if(is.na(rho1)==T){rho1=1}
    if(is.na(rho2)==T){rho2=1}
    if(is.na(rho3)==T){rho3=1}
    if(is.na(rho4)==T){rho4=1}
    if(is.na(rho5)==T){rho5=1}
    
    if(rho0>=1){beta0[j]=z0}else{coin0=rbinom(1,1,rho0)
    if(coin0==1){beta0[j]=z0}else{beta0[j]=beta0[j-1]}}
    
    if(rho1>=1){beta1[j]=z1}else{coin1=rbinom(1,1,rho1)
    if(coin1==1){beta1[j]=z1}else{beta1[j]=beta1[j-1]}}
    
    if(rho2>=1){beta2[j]=z2}else{coin2=rbinom(1,1,rho2)
    if(coin2==1){beta2[j]=z2}else{beta2[j]=beta2[j-1]}}
    
    if(rho3>=1){beta3[j]=z3}else{coin3=rbinom(1,1,rho3)
    if(coin3==1){beta3[j]=z3}else{beta3[j]=beta3[j-1]}}
    
    if(rho4>=1){beta4[j]=z4}else{coin4=rbinom(1,1,rho4)
    if(coin4==1){beta4[j]=z4}else{beta4[j]=beta4[j-1]}}
    
    if(rho5>=1){beta5[j]=z5}else{coin5=rbinom(1,1,rho5)
    if(coin5==1){beta5[j]=z5}else{beta5[j]=beta5[j-1]}}
    
    beta=c(beta0[j],beta1[j],beta2[j],beta3[j],beta4[j],beta5[j])
  }
  newbeta=cbind(beta0,beta1,beta2,beta3,beta4,beta5)
  return(newbeta)
}

mcmc1=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc2=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc3=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc4=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc5=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))
mcmc6=mcmc(rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25),rnorm(1,0,25))

beta_mcmc1=colMeans(mcmc1)
beta_mcmc2=colMeans(mcmc2)
beta_mcmc3=colMeans(mcmc3)
beta_mcmc4=colMeans(mcmc4)
beta_mcmc5=colMeans(mcmc5)

beta_mcmc=colMeans(rbind(beta_mcmc1,beta_mcmc2,beta_mcmc3,beta_mcmc4,beta_mcmc5))
beta_mcmc

newbeta_burn=newbeta[0:500,]
betab_burn=colMeans(newbeta_burn)
betab_burn

score(y,beta_mcmc,n,x)


####burn no difference on score
beta_mcmc1=colMeans(mcmc1[500:N,])
beta_mcmc2=colMeans(mcmc2[500:N,])
beta_mcmc3=colMeans(mcmc3[500:N,])
beta_mcmc4=colMeans(mcmc4[500:N,])
beta_mcmc5=colMeans(mcmc5[500:N,])

beta_mcmc=colMeans(rbind(beta_mcmc1,beta_mcmc2,beta_mcmc3,beta_mcmc4,beta_mcmc5))
score(y,beta_mcmc,n,x)

par(mfrow = c(2, 3))

plot(beta0,type="l", main = "constant term")
plot(beta1,type="l", main = "female")
plot(beta2,type="l", main = "college")
plot(beta3,type="l", main = "chronic")
plot(beta4,type="l", main = "main_habits")
plot(beta5,type="l", main = "label")

##beta0
plot(mcmc1[,1],type="l",col="red", main = "constant term")
lines(mcmc2[,1],col="blue")
lines(mcmc3[,1],col="green")
lines(mcmc4[,1],col="black")
lines(mcmc5[,1],col="yellow")


##beta1
plot(mcmc1[,2],type="l",col="red", main = "female")
lines(mcmc2[,2],col="blue")
lines(mcmc3[,2],col="green")
lines(mcmc4[,2],col="black")
lines(mcmc5[,2],col="yellow")


##beta2
plot(mcmc1[,3],type="l",col="red", main = "college")
lines(mcmc2[,3],col="blue")
lines(mcmc3[,3],col="green")
lines(mcmc4[,3],col="black")
lines(mcmc5[,3],col="yellow")


##beta3
plot(mcmc1[,4],type="l",col="red", main = "chronic")
lines(mcmc2[,4],col="blue")
lines(mcmc3[,4],col="green")
lines(mcmc4[,4],col="black")
lines(mcmc5[,4],col="yellow")


##beta4
plot(mcmc1[,5],type="l",col="red", main = "main_habits")
lines(mcmc2[,5],col="blue")
lines(mcmc3[,5],col="green")
lines(mcmc4[,5],col="black")
lines(mcmc5[,5],col="yellow")


##beta5
plot(mcmc1[,6],type="l",col="red", main = "label")
lines(mcmc2[,6],col="blue")
lines(mcmc3[,6],col="green")
lines(mcmc4[,6],col="black")
lines(mcmc5[,6],col="yellow")

hist(mcmc1[500:2000,6], freq=FALSE, col="grey", bor="darkgrey")
