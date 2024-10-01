rm(list=ls())
library(fields)
library(splines)
library(Matrix)
library(matrixcalc)
library(MASS)
library(RobustGaSP)
library(parallel)
library(invgamma)

load('emulation_ns_500_5.Rdata')

#define color palette for plotting
Temp.palette <- colorRampPalette(c("blue","grey","#DDDDDD"),space = "Lab")
resid.palette <- colorRampPalette(c("blue","#FFFFFF","red"),space = "Lab")

#x y coordinates
x=1:nrow(thickness.obs)
y=1:ncol(thickness.obs)

#define parameters to estimate
parameter.to.estimate <- 1:10

theta=param

p=nrow(mat.Y)
n=ncol(mat.Y)

a.sigma.nu=2
b.sigma.nu=1


nx=length(x)
ny=length(y)
xy=cbind(rep(x,times=ny),rep(y,each=nx))

library(plyr)
knots=as.matrix(expand.grid(seq(1,11,3),seq(1,11,3)))

J.d=nrow(knots) 

K.d=exp(-as.matrix(rdist(xy,knots))/4)

image(x,y,matrix(K.d[,2],nx,ny))

K.y=K.eta1

mean.Y=mean.Y1

png("principal_components_binary.png",width=1000,height=1000)
par(mfrow=c(3,4),mar=c(6,6,3,3),oma=c(2,0,0,2))
for(i in 1:J.y)
{
  K.y.draw=matrix(K.y[,i],length(x),length(y))
  image.plot(x,y,K.y.draw,col=tim.colors(50),main=paste("PC",i),xlab="x",ylab="y",cex.main=1.8,cex.lab=1.5,cex.axis=1.3)
}
dev.off()


thickness.obs.draw=thickness.obs
thickness.obs.draw[thickness.obs.draw==0]=NA

model.draw1=thickness[,,243]
model.draw1[model.draw1==0]=NA

model.draw2=thickness[,,394]
model.draw2[model.draw2==0]=NA

model.draw3=thickness[,,71]
model.draw3[model.draw3==0]=NA


image(thickness.obs,col=grey(seq(0,1,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',main='Observation')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "blue")
image(thickness.obs.draw,col=grey(seq(0,1,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)

points.draw=matrix(NA,nrow(thickness.obs.draw),ncol(thickness.obs.draw))
library(plyr)
x.draw=as.matrix(expand.grid(seq(1,11,3),seq(1,11,3)))
points.draw[x.draw]=2000
image(points.draw,col='red',xaxt='n',yaxt='n',add=TRUE)


#svd function for residuals
svd.resid=function(n,p,mat.Y,J.eta)
{
  svd.Y=svd(t(mat.Y), LINPACK = FALSE) #apply SVD
  K.eta=-svd.Y$u%*%diag(svd.Y$d)[,1:J.eta] #orthonormalized Principal component
  explained.var=sum(svd.Y$d[1:J.eta])/sum(svd.Y$d) #the proportion of explained variance
  list(K.eta=K.eta,explained.var=explained.var) #returen the outcomes
}

nx=length(x)
ny=length(y)


J.eta2=ncol(K.eta2)

J.d2=nrow(knots) #number of knots

K.d2=exp(-as.matrix(rdist(xy,knots))/4)

sigma2.Z=1000
rho.v=matrix(0.01,ncol(K.d2),ncol(K.d))




sigma.nu=1
sigma2.v=1000

a.sigma2.v=50
b.sigma2.v=500^2*(a.sigma2.v)/max(diag(K.d2%*%t(K.d2)))

plot(sqrt(seq(1,150000,500)),dinvgamma(seq(1,150000,500),a.sigma2.v,b.sigma2.v))

#parameter specifications
#non-informative prior for the observational error variance(sigma2.Z) and the partial sill(kappa): Inverse-Gamma (2,0.0002)
a.sigma2.Z=50
b.sigma2.Z=100^2*a.sigma2.Z

plot(sqrt(1:130000),dinvgamma(1:130000,a.sigma2.Z,b.sigma2.Z))

kappa.mat2=rep(NA,J.eta2)

for(i in 1:J.eta2)
{
  kappa.mat2[i]=gasp.result_height[[i]]@sigma2_hat
}

kappa.mat.update=kappa.mat2
scale.kappa=sqrt(kappa.mat2)

a.kappa=50
b.kappa.mat=kappa.mat2[1:J.eta2]*(a.kappa+1)


var.kappa=b.kappa.mat^2/((a.kappa-1)^2*(a.kappa-2))


mean.theta1=0.75
sigma.theta1=0.1
  
mean.theta2=0.3
sigma.theta2=0.15
  
mean.theta3=0.8
sigma.theta3=0.1
  
mean.theta5=0.4
sigma.theta5=0.15





##############################################################################
# initialization
##############################################################################

theta.star=unlist(theta[393,])

iter2=2000000

eta=mat.Y.red1[288,]
nu=rnorm(J.d,0,1)
sigma.nu=1

condtional1=function(theta.star)
{
  cond.mean=NULL
  cond.sd=NULL
  
  for(chain.no in 1:J.y)
  {
    predicted=predict.rgasp(gasp.result_binary[[chain.no]],matrix(theta.star,nrow=1,ncol=length(theta.star)))
    
    cond.mean[chain.no]=predicted$mean
    
    cond.sd[chain.no]=predicted$sd
  }
  
  return(list(cond.mean=cond.mean,cond.sd=cond.sd))
}

predicted=condtional1(theta.star)

cond.mean=predicted$cond.mean
cond.sd=predicted$cond.sd

log.CL.cond=function(cond.mean,cond.sd,sigma.nu,eta,nu)
{
  
  pi.eta=dnorm(eta,cond.mean,cond.sd,log=TRUE)
  
  k.y.eta=K.eta1%*%eta+mean.Y1
  
  k.d.nu=K.d%*%nu
  
  gamma=k.y.eta+k.d.nu
  
  sum(Z.binary*(gamma)-log(1+exp(gamma)))+sum(dnorm(nu,0,sigma.nu,log=TRUE))+sum(pi.eta)-(a.sigma.nu+1)*2*log(sigma.nu)-b.sigma.nu/sigma.nu^2
}


condtional2=function(theta.star)
{
  predicted.term=rep(NA,J.eta2)
  
  predicted.sigma=rep(NA,J.eta2)
  
  for(chain.no in 1:J.eta2)
  {
    Y.training=mat.Y.red2[,chain.no];
    
    X.new=cbind(1,predict(ns.list[[chain.no]][[1]],eta[1]))
    for(jj in 2:J.y) X.new=cbind(X.new,predict(ns.list[[chain.no]][[jj]],eta[jj])) 
    
    predicted=predict.rgasp(gasp.result_height[[chain.no]],matrix(theta.star,nrow=1,ncol=length(theta.star)),testing_trend=X.new)
    
    cond.mean=predicted$mean
    
    cond.sd=predicted$sd
    
    predicted.sigma[chain.no]=kappa.mat.update[chain.no]/kappa.mat2[chain.no]*cond.sd^2
    predicted.term[chain.no]=cond.mean
  }
  
  return(list(predicted.term=predicted.term,predicted.sigma=predicted.sigma))
}

predicted=condtional2(theta.star)

predicted.term=predicted$predicted.term
predicted.sigma=predicted$predicted.sigma


log.posterior2=function(kappa.mat.update,sigma2.Z,sigma2.v,predicted.sigma,predicted.term,rho.v,nu) 
{
  cond.mean.v=sqrt(sigma2.v)/sigma.nu*rho.v%*%nu
  sigma2.cond.v=sigma2.v*(diag(1,ncol(K.d2))-rho.v%*%t(rho.v))
  
  mu.height=K.eta2%*%predicted.term+matrix(mean.Y,n,1)+K.d2%*%cond.mean.v
  
  Sigma.eta.nu=as.matrix(bdiag(diag(c(predicted.sigma)),sigma2.cond.v))
  Inv.Part=solve(Sigma.eta.nu)+K2tK2/sigma2.Z
  Sigma.Z.inv=diag(1/sigma2.Z,nrow(K2))-1/sigma2.Z*K2%*%solve(Inv.Part,t(K2))/sigma2.Z
  
  Sigma.Z.log.det=determinant(solve(Sigma.eta.nu)+K2tK2/sigma2.Z,logarithm=TRUE)$modulus[1]+determinant(Sigma.eta.nu,logarithm=TRUE)$modulus[1]+log(sigma2.Z)*nrow(K2)
  
  Z.minus.mu=Z[Z>0]-mu.height[Z>0]
  
  log.lik=-0.5*Sigma.Z.log.det-0.5*Z.minus.mu%*%Sigma.Z.inv%*%Z.minus.mu
  
  log.prior=-(a.sigma2.v-1)*log(sigma2.v)-b.sigma2.v/sigma2.v-(a.sigma2.Z-1)*log(sigma2.Z)-b.sigma2.Z/sigma2.Z-(a.kappa+1)*sum(log(kappa.mat.update))-sum(b.kappa.mat/kappa.mat.update)
  
  log.lik+log.prior
}


Z=c(thickness.obs)
Z.binary=c(thickness.obs>0)


K2=cbind(K.eta2,K.d2)

K2=K2[Z>0,]

K2tK2=t(K2)%*%K2


#the matrices to store the result
mat.theta.star=matrix(NA,iter2,d)
mat.sigma.nu=matrix(NA,iter2,1)
mat.eta=matrix(NA,iter2,J.y)
mat.nu=matrix(NA,iter2,J.d)

mat.posterior1=matrix(NA,iter2,1)
mat.posterior2=matrix(NA,iter2,1)

theta.lower=apply(theta,2,min)
theta.upper=apply(theta,2,max)


#the matrices to store the result
mat.sigma2.v=matrix(NA,iter2,1)
mat.sigma2.Z=matrix(NA,iter2,1)
mat.kappa.update=matrix(NA,iter2,J.eta2)

par.to.estimate=1:length(theta)


init.time=proc.time() #store the starting time
posterior1=log.CL.cond(cond.mean,cond.sd,sigma.nu,eta,nu) #the posterior density at the initial values
posterior2=log.posterior2(kappa.mat.update,sigma2.Z,sigma2.v,predicted.sigma,predicted.term,rho.v,nu) #the posterior density at the initial values

#running MCMC
for(k in 1:iter2)
{
  eta.new <- eta+rnorm(J.y,0,0.01) #propose new theta.star
  posterior.new=log.CL.cond(cond.mean,cond.sd,sigma.nu,eta.new,nu)
  #compute the posterior a the new theta.star value
  ratio=exp(posterior.new-posterior1)  #the acceptance ratio
  if(ratio>runif(1))  #with the probability of the acceptance ratio...
  {
    eta=eta.new   #update theta
    posterior1=posterior.new     #update the current posterior density value
  }
  
  posterior.sum=posterior1+posterior2
  sigma.nu.new <- sigma.nu+rnorm(1,0,1)
  nu.new <- nu+rnorm(J.d,0,1)
  if(sigma.nu.new>0)  #check if proposed sigma2.Z is greater than 0
  {
    posterior1.new=log.CL.cond(cond.mean,cond.sd,sigma.nu.new,eta,nu.new)
    posterior2.new=log.posterior2(kappa.mat.update,sigma2.Z,sigma2.v,predicted.sigma,predicted.term,rho.v,nu.new)
    posterior.sum.new=posterior1.new+posterior2.new #compute the posterior a the new theta.star value
    ratio=exp(posterior.sum.new-posterior.sum)  #the acceptance ratio
    if(ratio>runif(1))  #with the probability of the acceptance ratio...
    {
      sigma.nu=sigma.nu.new
      nu=nu.new   #update theta
      posterior1=posterior1.new
      posterior2=posterior2.new     #update the current posterior density value
    }
  }
  
  mat.sigma.nu[k]=sigma.nu
  mat.nu[k,]=nu
  mat.eta[k,]=eta
  mat.posterior1[k]=posterior1
  mat.posterior2[k]=posterior2
  
  
  #sigma2.v
  sigma2.Z.new <- sigma2.Z+rnorm(1,0,1000) #propose new sigma2.Z
  sigma2.v.new <- sigma2.v+rnorm(1,0,1000) #propose new sigma2.Z
  if(prod(sigma2.Z.new>0)&prod(sigma2.v.new>0))  #check if proposed sigma2.Z is greater than 0
  {
    posterior.new=log.posterior2(kappa.mat.update,sigma2.Z.new,sigma2.v.new,predicted.sigma,predicted.term,rho.v,nu)  #compute the posterior a the new sigma2.Z value
    ratio=exp(posterior.new-posterior2) #the acceptance ratio
    if(ratio>runif(1))  #with the probability of the acceptance ratio...
    {
      sigma2.Z=sigma2.Z.new
      sigma2.v=sigma2.v.new    #update sigma2.Z
      posterior2=posterior.new   #update the current posterior density value
    }
  }
  
  
  rho.v.new <- rho.v+matrix(rnorm(length(rho.v),0,0.05),nrow(rho.v),ncol(rho.v)) #propose new sigma2.Z
  if(prod(rho.v.new>-1)&prod(rho.v.new<1)&is.positive.definite(round(diag(1,ncol(K.d2))-rho.v.new%*%t(rho.v.new),5),tol=1e-36))  #check if proposed sigma2.Z is greater than 0
  {
    posterior.new=log.posterior2(kappa.mat.update,sigma2.Z,sigma2.v,predicted.sigma,predicted.term,rho.v.new,nu)  #compute the posterior a the new sigma2.Z value
    ratio=exp(posterior.new-posterior2) #the acceptance ratio
    if(ratio>runif(1))  #with the probability of the acceptance ratio...
    {
      rho.v=rho.v.new;
      posterior2=posterior.new   #update the current posterior density value
    }
  }
  
  
  #kappa
  kappa.mat.update.new  <- kappa.mat.update
  kappa.mat.update.new  <- kappa.mat.update+mvrnorm(1,rep(0,J.eta2),diag(var.kappa))/10 #propose new sigma2.Z
  if(prod(kappa.mat.update.new>0))  #check if proposed sigma2.Z is greater than 0
  {
    posterior.new=log.posterior2(kappa.mat.update.new,sigma2.Z,sigma2.v,predicted.sigma,predicted.term,rho.v,nu)
    #compute the posterior a the new sigma2.Z value
    ratio=exp(posterior.new-posterior2) #the acceptance ratio
    if(ratio>runif(1))  #with the probability of the acceptance ratio...
    {
      kappa.mat.update=kappa.mat.update.new    #update sigma2.Z
      posterior2=posterior.new   #update the current posterior density value
    }
  }
  
  
  posterior.sum=posterior1+posterior2
  #theta
  theta.star.new=theta.star
  theta.star.new[par.to.estimate] <- theta.star[par.to.estimate]+rnorm(d,0,0.02)[par.to.estimate] #propose new theta.star #estimate Kv only. Other values are fix at the initial(and the best) values
  if(prod(apply(theta-0.02,2,min)<theta.star.new&theta.star.new<apply(theta+0.02,2,max)))  #check if theta.star is within the prior range
  {
    predicted.new=condtional1(theta.star.new)
    cond.mean.new=predicted.new$cond.mean
    cond.sd.new=predicted.new$cond.sd
    
    predicted.new=condtional2(theta.star.new)
    predicted.sigma.new=predicted.new$predicted.sigma
    predicted.term.new=predicted.new$predicted.term
    
    posterior.sum.prior=posterior.sum+dnorm(theta.star[1],mean.theta1,sigma.theta1,log=TRUE)+dnorm(theta.star[2],mean.theta2,sigma.theta2,log=TRUE)+dnorm(theta.star[3],mean.theta3,sigma.theta3,log=TRUE)+dnorm(theta.star[5],mean.theta5,sigma.theta5,log=TRUE)
    posterior1.new=log.CL.cond(cond.mean.new,cond.sd.new,sigma.nu,eta,nu)
    posterior2.new=log.posterior2(kappa.mat.update,sigma2.Z,sigma2.v,predicted.sigma.new,predicted.term.new,rho.v,nu)
    posterior.sum.new=posterior1.new+posterior2.new+dnorm(theta.star.new[1],mean.theta1,sigma.theta1,log=TRUE)+dnorm(theta.star.new[2],mean.theta2,sigma.theta2,log=TRUE)+dnorm(theta.star.new[3],mean.theta3,sigma.theta3,log=TRUE)+dnorm(theta.star.new[5],mean.theta5,sigma.theta5,log=TRUE) #compute the posterior a the new theta.star value
    ratio=exp(posterior.sum.new-posterior.sum.prior)  #the acceptance ratio
    if(ratio>runif(1))  #with the probability of the acceptance ratio...
    {
      theta.star=theta.star.new   #update theta
      
      cond.mean=cond.mean.new
      cond.sd=cond.sd.new
      
      predicted.sigma=predicted.sigma.new
      predicted.term=predicted.term.new
      
      posterior1=posterior1.new
      posterior2=posterior2.new     #update the current posterior density value
    }
  }
  
  mat.sigma2.v[k] = sigma2.v
  mat.sigma2.Z[k] = sigma2.Z  #store the current value of sigma2.Z
  mat.theta.star[k,]=theta.star #store the current value of theta.star
  mat.kappa.update[k,]=kappa.mat.update
  mat.posterior1[k]=posterior1
  mat.posterior2[k]=posterior2
  
  if((k%%30000)==0)
  {
    save.image('calibration_results_obs_scale_5.Rdata')
  }
}

elap.time=proc.time()-init.time #store the elapsed time
elap.time

save.image('calibration_results_obs_scale_5.Rdata')

palette <- colorRampPalette(c("blue","yellow","red"),space = "Lab")

example.settings=c(167,212)
mat.theta.star=mat.theta.star[1:k,]
theta.scale=rep(1,10)
set.seed(2015)
library(fields)
library(KernSmooth)
pdf("parameter_densities_obs_gray.pdf",width=12,height=12)
par( oma=c(3,5,3,1),mfrow=c(9,9),mar=c(1,1,0,0))
param.names=c("OCFAC",
              "CALV",
              "CRH",
              "CRHASE",
              "TAU",
              "LITH",
              "GEO",
              "SUBPIN",
              "USCH",
              "LAPSE")

param.units=rep("",10)

for(i in 1:9)
{
  for(j in 1:9)
  {
    if(i<(11-j))
    {
      mat.theta.draw=cbind(mat.theta.star[,11-j],mat.theta.star[,i])
      est <- bkde2D(mat.theta.draw, bandwidth=c((max(mat.theta.star[,11-j])-min(mat.theta.star[,11-j]))/6,(max(mat.theta.star[,i]*theta.scale[i])-min(mat.theta.star[,i]*theta.scale[i]))/6),gridsize=c(200,200),range.x=list(c(apply(theta,2,min)[11-j],apply(theta,2,max)[11-j]),c(apply(theta,2,min)[i]*theta.scale[i],apply(theta,2,max)[i]*theta.scale[i])))
      contour(est$x1,est$x2,est$fhat,labcex=0.9,xlim = c(0,1),ylim=c(0,1),xaxt='n',yaxt='n')
    } else
    {
      plot.new()
    }
    if( (i+j)==10 ){
      if(i==1) axis( 1, cex.axis=1.5, at=seq(0,1,0.2), labels = c(0.1,NA,round(10^(0.4*(log10(3.16)-log10(0.316))+log10(0.316)),1),NA,round(10^(0.8*(log10(3.16)-log10(0.316))+log10(0.316)),1),NA)) 
      if(i==2) axis( 1, cex.axis=1.1, at=seq(0,1,,5), labels = c(NA,"10^-8",NA,"10^-6",NA)) 
      if(i==3) axis( 1, cex.axis=1.1, at=seq(0,1,,4), labels = c(NA,"10^-2",NA,"1")) 
      if(i==4) axis( 1, cex.axis=1.5, at=seq(0,1,0.2), labels = 1:6) 
      if(i==5) axis( 1, cex.axis=1, at=seq(0,1,,3), labels = c('10^23','10^24','10^25')) 
      if(i==6) axis( 1, cex.axis=1.5, at=seq(0,1,,5), labels = c(NA,60,70,80,NA)) 
      if(i==7) axis( 1, cex.axis=1.5, at=seq(0,1,,5), labels = c(0,1,2,3,4)) 
      if(i==8) axis( 1, cex.axis=1.5, at=seq(0,1,0.2), labels = c(0.1,NA,round(10^(0.4*(log10(5)-log10(0.5))+log10(0.5)),1),NA,round(10^(0.8*(log10(5)-log10(0.5))+log10(0.5)),1),NA)) 
      if(i==9) axis( 1, cex.axis=0.9, at=seq(0,1,0.2), labels = c(-0.005,NA,-0.007,NA,-0.009,NA)) 
      mtext( line=3, side=1,param.names[11-j] )
      mtext( line=3, side=1,param.units[11-j] )}
    if( (11-j)==10 ){
      if(i==1) axis( 2, cex.axis=1.5, at=seq(0,1,0.2), labels = c(0.1,NA,round(10^(0.4*(log10(10)-log10(0.1))+log10(0.1)),1),NA,round(10^(0.8*(log10(10)-log10(0.1))+log10(0.1)),1),NA)) 
      if(i==2) axis( 2, cex.axis=1.5, at=seq(0,1,0.2), labels = c(0.1,NA,round(10^(0.4*(log10(3.16)-log10(0.316))+log10(0.316)),1),NA,round(10^(0.8*(log10(3.16)-log10(0.316))+log10(0.316)),1),NA)) 
      if(i==3) axis( 2, cex.axis=1.1, at=seq(0,1,,5), labels = c(NA,"10^-8",NA,"10^-6",NA)) 
      if(i==4) axis( 2, cex.axis=1.1, at=seq(0,1,,4), labels = c(NA,"10^-2",NA,"1")) 
      if(i==5) axis( 2, cex.axis=1.5, at=seq(0,1,0.2), labels = 1:6) 
      if(i==6) axis( 2, cex.axis=1, at=seq(0,1,,3), labels = c('10^23','10^24','10^25')) 
      if(i==7) axis( 2, cex.axis=1.5, at=seq(0,1,,5), labels = c(NA,60,70,80,NA)) 
      if(i==8) axis( 2, cex.axis=1.5, at=seq(0,1,,5), labels = c(0,1,2,3,4)) 
      if(i==9) axis( 2, cex.axis=1.5, at=seq(0,1,0.2), labels = c(0.1,NA,round(10^(0.4*(log10(5)-log10(0.5))+log10(0.5)),1),NA,round(10^(0.8*(log10(5)-log10(0.5))+log10(0.5)),1),NA)) 
      mtext( line=3, side=2,param.names[i])
      mtext( line=3, side=2,param.units[i])}
    
    if(i==1&j==1) mtext("2-D Posterior Densities for Input Parameters ", side=3, line=1, outer=TRUE, cex=1.2)
  }
}
dev.off()

