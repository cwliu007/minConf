dampedUpdate <- function(y,s,corrections,debug,old_dirs,old_stps,Hdiag){
  
  # B0 = eye(length(y))/Hdiag;

  row = dim(old_dirs)[1]
  row2 = dim(old_stps)[1]
  
  col = dim(old_dirs)[2]
  col2 = dim(old_stps)[2]
  
  if (col>=2){
    S = as.matrix(old_dirs[,2:col])
    Y = as.matrix(old_stps[,2:col])
  }else{
    col = 0
    col2 = 0
    S = matrix(0,row,col)
    Y = matrix(0,row2,col2)
  }

  k = dim(Y)[2]
  L = matrix(0,k,k)
  
  if (k>0){
    for (j in 1:k){
      if ((j+1)<=k){
        for (i in (j+1):k){
          L[i,j] = t(as.matrix(S[,i]))%*%as.matrix(Y[,j])
        }
      }
    }
  }

  temp = as.matrix(diag(t(S)%*%Y))
  temp2 = matrix(0,k,k)
  diag(temp2) <- temp
  D = as.matrix(temp2)
  if (dim(D)[1]==0) dim(D) <- c(0,0)
  N = cbind(S/Hdiag,Y)
  M = rbind(cbind(t(S)%*%S/Hdiag,L), cbind(t(L),-D))
  
  ys = t(y)%*%s
  if (col>=2){
    Bs = s/Hdiag - N%*%solve(M,t(N)%*%s)
  }else{
    Bs = s/Hdiag
  }
  sBs = t(s)%*%Bs
  
  eta = .02
  if (ys < (eta*sBs)){
    if (debug==1){
      print('Damped Update\n')
    }
    theta = min(max(0,((1-eta)*sBs)/(sBs - ys)),1)
    y = theta*y + (1-theta)*Bs
  }
  
  numCorrections = dim(old_dirs)[2]
  if (numCorrections < corrections){
    # Full Update
    old_dirs = cbind(old_dirs,s)
    old_stps = cbind(old_stps,y)
  }else{
    # Limited-Memory Update
    old_dirs = cbind(old_dirs[,2:corrections], s)
    old_stps = cbind(old_stps[,2:corrections], y)
  }
  
  # Update scale of initial Hessian approximation
  Hdiag = as.double( t(y)%*%s %*% solve(t(y)%*%y) )
  
  return(  list(old_dirs=old_dirs,old_stps=old_stps,Hdiag=Hdiag)  )
}