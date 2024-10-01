rm(list=ls())
###############################################################
#extracting gound ice thickness for RossIceShelf
###############################################################
load('thickness500.Rdata')

nx=50
ny=62

x.starting=87
y.starting=15

thickness=thickness[x.starting:(x.starting+(nx-1)),y.starting:(y.starting+ny-1),as.numeric(rownames(param.original))]*(1-mask[x.starting:(x.starting+(nx-1)),y.starting:(y.starting+ny-1),as.numeric(rownames(param.original))])

thickness.obs=thickness.obs[x.starting:(x.starting+(nx-1)),y.starting:(y.starting+ny-1)]*(1-mask.obs[x.starting:(x.starting+(nx-1)),y.starting:(y.starting+ny-1)])


n=nx*ny
p=dim(thickness)[3]

mat.Y=matrix(thickness,p,n,byrow=TRUE)

param=param.original

save(file='ground_ice_thickness_Ross_500.Rdata','thickness','thickness.obs','param.original','param','mat.Y')


#Finding the best run based on mse

mse=NULL
for(i in 1:500)
{
  mse[i]=sqrt(mean((mat.Y[i,]-thickness.obs)^2,na.rm=TRUE))
}
which.min(mse)

model.draw1=thickness[,,243]
model.draw1[model.draw1==0]=NA

model.draw2=thickness[,,313]
model.draw2[model.draw2==0]=NA

model.draw3=thickness[,,71]
model.draw3[model.draw3==0]=NA



###Generate Example Plots####################
par(mfrow=c(2,2),mar=c(1,1,3,1),oma=c(0,0,0,5))
thickness.obs.draw=thickness.obs
thickness.obs.draw[thickness.obs.draw==0]=NA


image(thickness.obs,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',main='Observation')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
image(thickness.obs.draw,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)

points.draw=matrix(NA,nrow(thickness.obs.draw),ncol(thickness.obs.draw))
library(plyr)
x.draw=as.matrix(expand.grid(seq(1,11,3),seq(1,11,3)))

points.draw[x.draw]=2000
image(points.draw,col='red',xaxt='n',yaxt='n',add=TRUE)

image(model.draw1,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',main='Parameter Setting #1')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
image(model.draw1,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)

image(model.draw2,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',main='Parameter Setting #2')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
image(model.draw2,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)

image(model.draw3,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',main='Parameter Setting #3')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
image(model.draw3,col=grey(seq(0,0.8,,50)),zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),xaxt='n',yaxt='n',add=TRUE)

library(fields)
par(oma=c(0,0,0,1.5))# reset margin to be much smaller.
image.plot( legend.only=TRUE,zlim=c(0,max(thickness.obs,model.draw1,model.draw2,model.draw3,na.rm=T)),col=grey(seq(0,0.8,,50)),legend.width = 2,main='height (m)') 


