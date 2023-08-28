knitr::opts_chunk$set(echo = TRUE)
require(deSolve)
require(rmarkdown)
require(ggplot2)
require(dplyr)
library(data.table)
library(magrittr)
library(fasstr)




seirmod=function(t, y, parms){
  #Pull state variables from y vector
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]
  x=y[5]
  #Pull the required parameter values from the parms vector
  beta=parms["beta"]
  kappa=parms["kappa"]
  gamma=parms["gamma"]
  mu=parms["mu"]
  #Define the equations 
  lambda = (beta*(1-mu*(0.5*tanh(x)+0.5)))
  dS = -lambda*S*I
  dE = lambda*S*I - kappa*E
  dI = kappa*E - gamma*I
  dR = gamma*I
  dx = (-0.01 + beta*I)
  res=c(dS, dE, dI, dR, dx)
  #Return list of gradients
  list(res)
}


times  = seq(0, 10000, by=1)
parms  = c(beta=1.1, kappa=1/6, gamma=1/6, mu=1)
start = c(S=0.95, E= 0, I=0.05, R = 0, x = 0.8)


out = ode(y = start, times = times, func = seirmod, 
          parms = parms, maxsteps = 1000, method = "ode45")
out = as.data.frame(out)

out$x <- 0.5 * tanh(out$x) + 0.5


plot(x = out$time, y = out$I, col = "black", ylab = "Fraction of individuals", 
     xlab = "Time (days)", type = "l", xlim = c(0, 1000), ylim = c(0,1))
    lines(x = out$time, y = out$S, col =2)
    lines(x = out$time, y = out$E, col =3)
    lines(x = out$time, y = out$R, col =4)
    lines(x = out$time, y = out$x)