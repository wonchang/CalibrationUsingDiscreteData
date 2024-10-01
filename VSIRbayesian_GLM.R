rm(list=ls())
library(deSolve)
library(invgamma)

T <- 80
iter <- 100000
N <- 50000

# SIR model
SIR <- function(t, s, theta) {
  
  with(as.list(c(s, theta)), {
    dS <- beta*I*S - vaccination*S 
    dI <- beta*I*S - recovery*I 
    dR <- recovery*I + vaccination*S 
    return(list(c(dS,dI,dR)))
  })
}

### initial state
##  Susceptible 1, Infected 0.001, Recovered 0
init       <- c(S = 1, I = 0.001, R = 0)
## beta: infection parameter; gamma: recovery parameter
theta <- c(beta = 0.12, recovery = 0.025, vaccination = 0.06)
## time steps
times <- seq(0,T-1,by = 1)

## Solve using ode function
out <- ode(y = init, times=times, func=SIR, parms=theta)
out

##Creating observational data
set.seed(2022)
#Z <- out[,3] + rnorm(T,0,0.002)
plot(out[,4])

Z <- rpois(T,N*out[,3])
plot(Z)

Z.I <- tibble(time=times,I=N*out$I,Z)

Z.I %>% ggplot() + geom_line(aes(y=Z,x=time),color='black',size=1) +
  geom_line(aes(y=I,x=time),color='blue',alpha=0.5,size=1) +
  ylab("# of Infected") +
  theme(text = element_text(size = 25), 
        panel.background = element_rect(fill=FALSE,color='black')) 

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
  theta <- c(beta = param[1], recovery = param[2], vaccination = param[3])
  out <- ode(y = init, times=times, func=SIR, parms=theta)
  return(out)
}

##Log posterior density
posterior <- function(all_param)
{
  n <- length(Z)
  Y_theta <- model(all_param[1:3]) 
  
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
mcmc.sample=metrop(posterior,c(0.13,0.023,0.059),scale=c(0.001,0.001,0.001),nbatch=iter)$batch
proc.time()-init.time

##Checking the generated chain

mcmc.sample <- as.tibble(mcmc.sample)
colnames(mcmc.sample) <- c('beta','recovery','vaccination')
plot(mcmc.sample[[3]]) 
plot(mcmc.sample[[1]],mcmc.sample[[2]]) 

saveRDS(mcmc.sample,file='mcmc_sample.rds')
mcmc.sample <- readRDS('mcmc_sample.rds')

# Density plot
mcmc.sample %>% ggplot(aes(x=beta))+
  geom_density(color="darkblue", fill="lightblue")

par(mfrow=c(2,2))
# Hexbin chart 
mcmc.sample %>% slice(which(row_number() %% 10 == 1)) %>%
  ggplot(aes(x=beta, y=vaccination) ) +
  geom_bin_2d(bins = 20, aes(fill = ..density..)) +
  scale_fill_continuous(type = "viridis") +
  theme(text = element_text(size = 25), 
        panel.background = element_rect(fill=FALSE,color='black')) +
  geom_hline(yintercept=0.06) +
  geom_vline(xintercept=0.12) 

mcmc.sample %>% slice(which(row_number() %% 10 == 1)) %>%
  ggplot(aes(x=beta, y=recovery) ) +
  geom_bin_2d(bins = 20, aes(fill = ..density..)) +
  scale_fill_continuous(type = "viridis") +
  theme(text = element_text(size = 25), 
        panel.background = element_rect(fill=FALSE,color='black')) +
  geom_hline(yintercept=0.025) +
  geom_vline(xintercept=0.12) 

mcmc.sample %>% slice(which(row_number() %% 10 == 1)) %>%
  ggplot(aes(x=vaccination, y=recovery) ) +
  geom_bin_2d(bins = 20, aes(fill = ..density..)) +
  scale_fill_continuous(type = "viridis") +
  theme(text = element_text(size = 25), 
        panel.background = element_rect(fill=FALSE,color='black')) +
  geom_hline(yintercept=0.025) +
  geom_vline(xintercept=0.06) 


##Run the model for an example chain member
i=iter
theta <- c(mcmc.sample[i,1], mcmc.sample[i,2], mcmc.sample[i,3], mcmc.sample[i,4], mcmc.sample[i,5])
out <- ode(y = init, times=times, func=SIR, parms=theta)

plot(mcmc.sample[1:iter,c(1,3)],type='l')

library(tidyverse)
out <- as.tibble(out)
out %>% ggplot() + geom_line(aes(y=S,x=time,color="S"),size=2) +
  geom_line(aes(y=I,x=time,color="I"),color='blue',size=2) + 
  geom_line(aes(y=R,x=time,color="R"),color='red',size=2) +
  scale_color_manual(name='',values = c("S"="black","I"="blue","R"="red"))

#Calibrated runs based on MCMC sample
out.mat <- array(NA,dim=c(T,iter/10))
count <- 0
for(i in seq(1,iter,10))
{
  count <- count+1
  theta <- c(mcmc.sample[i,1], mcmc.sample[i,2], mcmc.sample[i,3])
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

mat.I %>% ggplot() + geom_line(aes(y=Z,x=time),color='black') +
  theme(text = element_text(size = 25), 
        panel.background = element_rect(fill=FALSE,color='black')) +
  ylab("# of Infected") +
  geom_line(aes(y=lower,x=time),color='blue',alpha=0.5,lty=2) + 
  geom_line(aes(y=upper,x=time),color='blue',alpha=0.5,lty=2) +
  geom_line(aes(y=mean,x=time),color='blue')

 