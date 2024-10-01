rm(list=ls())


load('ground_ice_thickness_Ross_500.Rdata')

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

set.seed(5)
pca.result.b=sparse.logistic.pca(mat.Y.binary,start.mu=rnorm(n),start.A=matrix(rnorm(J.y*p),p,J.y),start.B=matrix(rnorm(n*J.y),n,J.y))

mat.Y.red1=pca.result.b$A
K.eta1=pca.result.b$B
mean.Y1=pca.result.b$mu

######################################################################
# Dimension Reduction for Ice Thickness                              #
######################################################################

library(pcaMethods)

mat.Y.log=mat.Y

mat.Y.log[mat.Y<1]=log(mat.Y[mat.Y<1])+1

mat.Y.log[mat.Y.log==-Inf]=NA

log.results=ppca(mat.Y.log,nPcs=10,threshold = 1e-03, maxIterations = 100000,seed=2018)

mat.Y.red2=scores(log.results)
K.eta2=loadings(log.results)

J.eta=ncol(K.eta2)


######################################################################
# Emulation for Ice Thickness                                        #
######################################################################
mat.Y.red=mat.Y.red1

source('PCAemulator_500.R')

######################################################################
# Emulation for Binary Patterns                                      #
######################################################################

mat.Y.red=mat.Y.red2

source('PCAemulatorConditional_ns_500_df=1.R')

simul=0

save.image('emulation_ns_500_5.Rdata')


