lbfgs <- function(g,s,y,Hdiag){
  # BFGS Search Direction
  # 
  # This function returns the (L-BFGS) approximate inverse Hessian,
  # multiplied by the gradient
  # 
  # If you pass in all previous directions/sizes, it will be the same as full BFGS
  # If you truncate to the k most recent directions/sizes, it will be L-BFGS
  # 
  # s - previous search directions (p by k)
  # y - previous step sizes (p by k)
  # g - gradient (p by 1)
  # Hdiag - value of initial Hessian diagonal elements (scalar)
  
  p = dim(s)[1]
  k = dim(s)[2]
  
  ro = rep(0.0,k)
  if (k>=1){
    for (i in 1:k){
      ro[i] = 1 / (t(as.matrix(y[,i]))%*%as.matrix(s[,i]))
    }
  }

  
  q = matrix(0,p,k+1)
  r = matrix(0,p,k+1)
  al = rep(0,k)
  be = rep(0,k)
  
  q[,k+1] = g
  
  if (k>=1){
    for (i in k:1){
      al[i] = ro[i]*t(as.matrix(s[,i]))%*%as.matrix(q[,i+1])
      q[,i] = q[,i+1] - al[i] * y[,i]
    }
  }

  
  # Multiply by Initial Hessian
  r[,1] = Hdiag*q[,1]
  
  if (k>=1){
    for (i in 1:k){
      be[i] = ro[i]*t(as.matrix(y[,i]))%*%as.matrix(r[,i])
      r[,i+1] = r[,i] + s[,i]*(al[i]-be[i])
    }
  }

  d = r[,k+1]
}
