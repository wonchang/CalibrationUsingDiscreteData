rm(list=ls())
library(deSolve)
library(invgamma)

T <- 300
iter <- 2000000
N <- 5000

# SIR model
SIR <- function(t, s, theta) {
  
  with(as.list(c(s, theta)), {
    dS <- birth - beta*I*S - vaccination*S - death*S
    dI <- beta*I*S - recovery*I - death*I
    dR <- recovery*I + vaccination*S - death*R
    return(list(c(dS,dI,dR)))
  })
}

### initial state
##  Susceptible 1, Infected 0.001, Recovered 0
init       <- c(S = 1, I = 0.001, R = 0)
## beta: infection parameter; gamma: recovery parameter
theta <- c(beta = 0.62, recovery = 0.018, death = 0.002, birth = 0.002, vaccination = 0.08)
## time steps
times <- seq(0,T-1,by = 1)

## Solve using ode function
out <- ode(y = init, times=times, func=SIR, parms=theta)

##Creating observational data
set.seed(2022)
#Z <- out[,3] + rnorm(T,0,0.002)
Z <- rpois(T,N*out[,3])

## Plotting underlying true process
library(tidyverse)
out <- as.tibble(out)
out %>% ggplot() + geom_line(aes(y=S,x=time,color="S"),size=2) +
  geom_line(aes(y=I,x=time,color="I"),color='blue',size=2) + 
  geom_line(aes(y=R,x=time,color="R"),color='red',size=2) +
  scale_color_manual(name='',values = c("S"="black","I"="blue","R"="red"))

##Forward model function
model <- function(param)
{
  theta <- c(beta = param[1], recovery = param[2], death = param[3], birth = param[4], vaccination = param[5])
  out <- ode(y = init, times=times, func=SIR, parms=theta)
  return(out)
}

##Log posterior density
posterior <- function(all_param)
{
  n <- length(Z)
  Y_theta <- model(all_param[1:5]) 
  
  within_01 <- prod(Y_theta[,2]>=0)*prod(Y_theta[,3]>=0)*prod(Y_theta[,4]>=0)*
    prod(Y_theta[,2]<=1)*prod(Y_theta[,3]<=1)*prod(Y_theta[,4]<=1)
  
  if(prod(all_param>0)==1&prod(all_param<1)==1&within_01==1)
  {
    log_p <- sum(Z*log(N*Y_theta[,3])-N*Y_theta[,3])
  } else
  {
    log_p <- -Inf
  }
  
  return(log_p)
}

#Running MCMC
init.time=proc.time()
library(mcmc)
mcmc.sample=metrop(posterior,c(0.6,0.019,0.0021,0.0019,0.075),scale=c(0.001,0.001,0.001,0.001,0.001),nbatch=iter)$batch
proc.time()-init.time

##Checking the generated chain
#plot(mcmc.sample[,3])
#plot(mcmc.sample[,1],mcmc.sample[,5])

mcmc.sample <- as.tibble(mcmc.sample)
colnames(mcmc.sample) <- c('beta','recovery','death','birth','vaccination')
plot(mcmc.sample[[3]]) 

saveRDS(mcmc.sample,file='mcmc_sample.rds')
mcmc.sample <- readRDS('mcmc_sample.rds')

# Density plot
mcmc.sample %>% ggplot(aes(x=beta))+
  geom_density(color="darkblue", fill="lightblue")

# Hexbin chart 
mcmc.sample %>% slice(which(row_number() %% 100 == 1)) %>%
  ggplot(aes(x=beta, y=vaccination) ) +
  geom_hex(bins = 25) +
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept=0.08) +
  geom_vline(xintercept=0.62) +
  theme_bw()

# Hexbin chart 
mcmc.sample %>% slice(which(row_number() %% 100 == 1)) %>%
  ggplot(aes(x=recovery, y=death) ) +
  geom_hex(bins = 25) +
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept=0.002) +
  geom_vline(xintercept=0.018) +
  theme_bw()

##Run the model for an example chain member
i=iter
theta <- c(mcmc.sample[i,1], mcmc.sample[i,2], mcmc.sample[i,3], mcmc.sample[i,4], mcmc.sample[i,5])
out <- ode(y = init, times=times, func=SIR, parms=theta)

plot(mcmc.sample[1:iter,c(1,4)],type='l')

library(tidyverse)
out <- as.tibble(out)
out %>% ggplot() + geom_line(aes(y=S,x=time,color="S"),size=2) +
  geom_line(aes(y=I,x=time,color="I"),color='blue',size=2) + 
  geom_line(aes(y=R,x=time,color="R"),color='red',size=2) +
  scale_color_manual(name='',values = c("S"="black","I"="blue","R"="red"))

#Calibrated runs based on MCMC sample
out.mat <- array(NA,dim=c(T,iter/100))
count <- 0
for(i in seq(1,iter,100))
{
  count <- count+1
  theta <- c(mcmc.sample[i,1], mcmc.sample[i,2], mcmc.sample[i,3], mcmc.sample[i,4], mcmc.sample[i,5])
  It <- ode(y = init, times=times, func=SIR, parms=theta)[,3]
  out.mat[,count] <- rpois(T,N*It)
}

saveRDS(out.mat,file='out_mat.rds')
out.mat <- readRDS('out_mat.rds')

##Plotting I
time.mat.I <- times
mean.mat.I <- apply(out.mat,1,mean)
lower.mat.I <- apply(out.mat,1,quantile,probs=0.025)
upper.mat.I <- apply(out.mat,1,quantile,probs=0.975)
mat.I <- tibble(time=time.mat.I,mean=mean.mat.I,lower=lower.mat.I,upper=upper.mat.I,Z=Z)

mat.I %>% ggplot() + geom_line(aes(y=Z,x=time),color='black',size=0.5) +
  geom_line(aes(y=lower,x=time),color='blue',alpha=0.5,size=0.5) + 
  geom_line(aes(y=upper,x=time),color='blue',alpha=0.5,size=0.5) +
  geom_line(aes(y=mean,x=time),color='blue',size=0.5)

 