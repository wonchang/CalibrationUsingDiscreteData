##############################################################
# Load model output and obs data
##############################################################

rm(list=ls())

load('ground_ice_thickness_Ross_500.Rdata')

J.eta=10

#need to remove a model run with all NA output
param=param[-which(apply(is.na(mat.Y),1,sum)>0),]
mat.Y=mat.Y[-which(apply(is.na(mat.Y),1,sum)>0),]

#########################################################################
# Defining the dimensionality of matrices                               #
# n : The number of spatial observations in UVic outputs                #
# p : The number of parameter settings in UVic outputs                  #
# d : number of input parameters
#########################################################################

n=ncol(mat.Y) #number of spatial locations
p=nrow(mat.Y) #number of model runs
d=ncol(param)

###################################################################
# compute the distance matrix for each parameter                  #
# C : squared distance matrices for covariance computation        #
###################################################################

C=array(NA,dim=c(p,p,d))

for(i in 1:d)
{
  SS=param[,i]
  AA = matrix(SS,length(SS),length(SS))
  BB = t(AA)
  C[,,i]=abs(AA-BB)
}

######################################################################
# Dimeionsl Reduction for Binary Outcome                             #
######################################################################

mat.Y.binary=mat.Y

mat.Y.binary[mat.Y>0]=1

J.y=10

source("SparseLogisticPCA.r")

pca.result.b=sparse.logistic.pca(mat.Y.binary,start.mu=rnorm(n),start.A=matrix(rnorm(J.y*p),p,J.y),start.B=matrix(rnorm(n*J.y),n,J.y))

mat.Y.red1=pca.result.b$A
K.eta1=pca.result.b$B
mean.Y1=pca.result.b$mu


######################################################################
# Emulation for Ice Thickness                                        #
######################################################################
mat.Y.red=mat.Y.red1

source('PCAemulator_500.R')

library(RobustGaSP)

mean.Y.red1 <- matrix(NA,nrow(mat.Y.red1),ncol(mat.Y.red1))

for(j in 1:J.y) 
{
  for(i in 1:10)
  {
    if(i<10)
    {
      leave.out=50*(i-1)+1:50
    } else
    {
      leave.out=451:499
    }
    
    q=length(leave.out)
    
    m.ppgasp_full=gasp.result_binary[[j]]
    
    zeta.hat=m.ppgasp_full@nugget.est
    
    mu.hat=m.ppgasp_full@theta_hat
    
    kappa.hat=m.ppgasp_full@sigma2_hat
    
    phi.hat=m.ppgasp_full@beta_hat
    
    D.array=array(NA,dim=c(p-q,p-q,d))
    
    D.array1=array(NA,dim=c(q,p-q,d))
    
    Cov11=matrix(1,p-q,p-q)
    
    for(k in 1:d) D.array[,,k]=abs(matrix(param[-leave.out,k],p-q,p-q)-t(matrix(param[-leave.out,k],p-q,p-q)))
    
    for(k in 1:d) D.array1[,,k]=abs(matrix(param[-leave.out,k],q,p-q,byrow=TRUE)-matrix(param[leave.out,k],q,p-q))
    
    sum11=0
    for(k in 1:d) sum11=sum11+D.array[,,k]*phi.hat[k]
    
    Cov11=kappa.hat*(exp(-sum11)+diag(zeta.hat,p-q))
    
    Cov12=matrix(1,q,p-q)
    sum12=0
    for(k in 1:d) sum12=sum12+D.array1[,,k]*phi.hat[k]
    
    Cov12=kappa.hat*exp(-sum12)
    
    Cond.mean=rep(mu.hat,q)+Cov12%*%solve(Cov11,mat.Y.red1[-leave.out,j]-rep(mu.hat,p-q))
    
    mean.Y.red1[leave.out,j]=Cond.mean
    #sd.Y.red1[leave.out,j]=m_pred1$sd
  }
}

non.zero.location <- apply(mat.Y==0,2,sum) != p  

predicted.Y.binary=t(K.eta1%*%t(mean.Y.red1)+mean.Y1)>0

#computing accuracy
mean(abs( (predicted.Y.binary[,non.zero.location]>0) == (mat.Y[,non.zero.location]>0)))

#true positive
true.positive=sum((predicted.Y.binary[,non.zero.location]==TRUE)&(mat.Y[,non.zero.location]>0))/sum(mat.Y[,non.zero.location]>0)
round(true.positive,3)
#true negative
true.negative=sum((predicted.Y.binary[,non.zero.location]==FALSE)&(mat.Y[,non.zero.location]==0))/sum(mat.Y[,non.zero.location]==0)
round(true.negative,3)


mae=rep(NA,25)
count=0
p.original=p

param.original=param

mat.Y.red1.original=mat.Y.red1

init.time=proc.time()
count=count+1
mean.Y.red2=matrix(NA,p.original,J.eta)

predicted.heights=matrix(NA,p.original,n)


for(i in 1:10)
{
  if(i<10)
  {
    leave.out=50*(i-1)+1:50
  } else
  {
    leave.out=451:499
  }
  
  q=length(leave.out)
  
  p=p.original-q
  
  ######################################################################
  # Dimension Reduction for Ice Thickness                              #
  ######################################################################
  
  library(pcaMethods)
  
  mat.Y.log=mat.Y
  
  mat.Y.log[mat.Y<1]=log(mat.Y[mat.Y<1])+1
  
  mat.Y.log[mat.Y.log==-Inf]=NA
  
  log.results=ppca(mat.Y.log[-leave.out,],nPcs=J.eta,threshold = 1e-03, maxIterations = 100000,seed=2018)
  
  mat.Y.red2=scores(log.results)
  K.eta2=loadings(log.results)
  
  param=param.original[-leave.out,]
  
  mat.Y.red1=mat.Y.red1.original[-leave.out,]
  
  mat.Y.red=mat.Y.red2
  
  source('PCAemulatorConditional_ns_500_df=1.R')
  
  for(j in 1:J.eta) 
  {
    q=length(leave.out)
    
    m.ppgasp_full=gasp.result_height[[j]]
    
    zeta.hat=m.ppgasp_full@nugget.est
    
    mu.hat=m.ppgasp_full@theta_hat
    
    kappa.hat=m.ppgasp_full@sigma2_hat
    
    phi.hat=m.ppgasp_full@beta_hat
    
    D.array=array(NA,dim=c(p,p,d))
    
    D.array1=array(NA,dim=c(q,p,d))
    
    Cov11=matrix(1,p,p)
    
    param.leave.out=param.original[leave.out,]
    
    for(k in 1:d) D.array[,,k]=abs(matrix(param[,k],p,p)-t(matrix(param[,k],p,p)))
    
    for(k in 1:d) D.array1[,,k]=abs(matrix(param[,k],q,p,byrow=TRUE)-matrix(param.leave.out[,k],q,p))
    
    sum11=0
    for(k in 1:d) sum11=sum11+D.array[,,k]*phi.hat[k]
    
    Cov11=kappa.hat*(exp(-sum11)+diag(zeta.hat,p))
    
    Cov12=matrix(1,q,p)
    sum12=0
    for(k in 1:d) sum12=sum12+D.array1[,,k]*phi.hat[k]
    
    Cov12=kappa.hat*exp(-sum12)
    
    X.new=cbind(1,predict(ns.list[[j]][[1]],mean.Y.red1[leave.out,1]))
    for(jj in 2:J.y) X.new=cbind(X.new,predict(ns.list[[j]][[jj]],mean.Y.red1[leave.out,jj]))  
    
    X=rep(1,p)
    for(jj in 1:J.y) X=cbind(X,splines::ns(mat.Y.red1[,jj],df=df.ns.matrix[j]))
    
    Cond.mean=X.new%*%mu.hat+Cov12%*%solve(Cov11,mat.Y.red2[,j]-X%*%mu.hat)
 
    mean.Y.red2[leave.out,j]=Cond.mean
  }
  predicted.heights[leave.out,]<-t(K.eta2%*%t(mean.Y.red2[leave.out,]))
}

non.zero.location <- apply(mat.Y==0,2,sum) < 499 

predicted.heights[predicted.heights<1]=exp(predicted.heights[predicted.heights<1]-1)

predicted.Y=predicted.heights

mat.combined = predicted.Y*predicted.Y.binary

mae[count]=mean(abs(mat.Y[,non.zero.location]-mat.combined[,non.zero.location]))
print(proc.time()-init.time)


pdf('leave1out_2Dresult_all.pdf',width=9,height=4.5)

non.zero.location <- apply(mat.Y==0,2,sum) < p.original

predicted.heights[predicted.heights<1]=exp(predicted.heights[predicted.heights<1]-1)

mat.combined = predicted.Y*predicted.Y.binary

mean(abs(mat.Y[,non.zero.location]-mat.combined[,non.zero.location]))

library(fields)

for(leave.out in 1:p.original)
{
  true <- mat.Y[leave.out,]
  
  predicted.combined <- predicted.Y[leave.out,]*predicted.Y.binary[leave.out,]
  
  true.draw=matrix(true,dim(thickness)[1],dim(thickness)[2])
  true.draw[true.draw==0]=NA
  
  predicted.combined.draw=matrix(predicted.combined,dim(thickness)[1],dim(thickness)[2])
  predicted.combined.draw[predicted.combined.draw==0]=NA
  
  par(mfrow=c(1,2),mar=c(1,1,3,1),oma=c(0,0,0,5))
  image(true.draw,col=grey(seq(0,1,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',main='Original')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "blue")
  image(true.draw,col=grey(seq(0,1,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)
  
  image(predicted.combined.draw,col=grey(seq(0,1,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',main='Emulated')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "blue")
  image(predicted.combined.draw,col=grey(seq(0,1,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)
  
  par(oma=c(0,0,0,1.5))# reset margin to be much smaller.
  image.plot( legend.only=TRUE,zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),col=grey(seq(0,1,,50)),legend.width = 2) 
}

dev.off()



pdf('leave1out_2Dresult_three.pdf',width=9,height=4.5)

non.zero.location <- apply(mat.Y==0,2,sum) < p.original

predicted.heights[predicted.heights<1]=exp(predicted.heights[predicted.heights<1]-1)

mat.combined = predicted.Y*predicted.Y.binary

mean(abs(mat.Y[,non.zero.location]-mat.combined[,non.zero.location]))

library(fields)

for(leave.out in c(313,314,322))
{
  true <- mat.Y[leave.out,]
  
  predicted.combined <- predicted.Y[leave.out,]*predicted.Y.binary[leave.out,]
  
  true.draw=matrix(true,dim(thickness)[1],dim(thickness)[2])
  true.draw[true.draw==0]=NA
  
  predicted.combined.draw=matrix(predicted.combined,dim(thickness)[1],dim(thickness)[2])
  predicted.combined.draw[predicted.combined.draw==0]=NA
  
  par(mfrow=c(1,2),mar=c(1,1,3,1),oma=c(0,0,0,5))
  image(true.draw,col=grey(seq(0,0.8,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',main='Original')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  image(true.draw,col=grey(seq(0,0.8,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)
  
  image(predicted.combined.draw,col=grey(seq(0,0.8,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',main='Emulated')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  image(predicted.combined.draw,col=grey(seq(0,0.8,,50)),zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)
  
  par(oma=c(0,0,0,1.5))# reset margin to be much smaller.
  image.plot( legend.only=TRUE,zlim=c(0,max(true.draw,predicted.combined.draw,na.rm=T)),col=grey(seq(0,1,,50)),legend.width = 2) 
}

dev.off()



save.image('EmulationCheckResults.Rdata')