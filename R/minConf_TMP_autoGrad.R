autoGrad <- function(x,type,numDiff,funObj,...){
  # [f,g] = autoGrad(x,useComplex,funObj,varargin)
  # 
  # Numerically compute gradient of objective function from function values
  # 
  # type =
  # 1 - forward-differencing (p+1 evaluations)
  # 3 - central-differencing (more accurate, but requires 2p evaluations)
  # 2 - complex-step derivative (most accurate and only requires p evaluations, but only works for certain objectives)
  
  p = length(x)
  mu = 1e-150
  
  if (numDiff == 2){ # Use Complex Differentials
    diff = rep(0,p)
    e_j = rep(0,p)
    for (j in 1:p){
      e_j[j] = 1
      diff[j] = funObj(x + mu*1i*e_j,...)
      e_j[j] = 0
    }
    ff = mean(Re(diff))
    g = Im(diff)/mu
  }else if(numDiff == 3){
    mu = 2*sqrt(1e-12)*(1+norm(x,"2"))/norm(p,"2")
    diff1 = rep(0,p)
    diff2 = rep(0,p)
    e_j = rep(0,p)
    for (j in 1:p){
      e_j[j] = 1
      diff1[j] = funObj(x + mu*e_j,...)
      diff2[j] = funObj(x - mu*e_j,...)
      e_j[j] = 0
    }
    ff = mean(c(diff1,diff2))
    g = (diff1 - diff2)/(2*mu)
  }else{ # Use Finite Differencing
    ff = as.double(funObj(x,...))
    mu = 2*sqrt(1e-12)*(1+norm(x,"2"))/norm(p,"2")
    diff = rep(0,p)
    e_j = rep(0,p)
    for (j in 1:p){
      e_j[j] = 1
      diff[j] = funObj(x + mu*e_j,...)
      e_j[j] = 0
    }
    g = (diff-ff)/mu
  }
  return(list(f=ff,g=g))
}