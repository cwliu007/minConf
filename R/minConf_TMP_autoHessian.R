autoHessian <- function(x,type,numDiff,funObj,...){
  # Numerically compute Hessian of objective function from function values
  # 
  # type =
  # 1 - forward-differencing
  # 3 - central-differencing

  p = length(x)
  
  eps = .Machine$double.eps
  
  H = matrix(0,p,p)
  
  if (numDiff == 2){ # Use Complex Differentials

  }else if(numDiff == 3){
    logp = funObj(x,...)
    ma = cbind(abs(x),rep(1,p))
    h2 = eps^(1/4)*sign(x)*apply(ma,1,max)
    h2[h2==0] = eps^(1/4)
    e_i = rep(0,p)
    e_j = rep(0,p)
    for (i in 1:p){
      e_i[i] = h2[i]
      for (j in i:p){
        e_j[j] = h2[j]
        if (i == j){
          Ui3 = funObj(x+e_i,...)
          Ui4 = funObj(x-e_i,...)
          FD = (Ui3 - 2*logp + Ui4)/(h2[i]^2)
        }else{
          Ui1 = funObj(x+e_i+e_j,...)
          Ui2 = funObj(x-e_i+e_j,...)
          Ui3 = funObj(x+e_i-e_j,...)
          Ui4 = funObj(x-e_i-e_j,...)
          FD = (Ui1 - Ui2 - Ui3 + Ui4)/(4*e_i[i]*e_j[j])
        }
        H[i,j] <- H[j,i] <- FD
        e_j[j] = 0
      }
      e_i[i] = 0
    }
    return(H)
    

  }else{ # Use Finite Differencing
    
    logp = funObj(x,...)
    ma = cbind(abs(x),rep(1,p))
    h2 = eps^(1/4)*sign(x)*apply(ma,1,max)
    h2[h2==0] = eps^(1/4)
    e_i = rep(0,p)
    e_j = rep(0,p)
    for (i in 1:p){
      e_i[i] = h2[i]
      logp_k = funObj(x+e_i,...)
      for (j in i:p){
        e_j[j] = h2[j]
        logp_h = funObj(x+e_j,...)
        
        logp_kh = funObj(x+e_i+e_j,...)
        
        FD = (logp_kh - logp_h - logp_k + logp)/(e_i[i]*e_j[j])
        
        H[i,j] <- H[j,i] <- FD
        
        e_j[j] = 0
      }
      e_i[i] = 0
    }

  }
  return(H)
}