lbfgsUpdate <- function(y,s,corrections,debug,old_dirs,old_stps,Hdiag){
  ys = t(y)%*%s
  if (ys > 1e-10){
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
    Hdiag = as.double(ys / (t(y)%*%y))
  }else{
    if (debug==1) print('Skipping Update\n')
  }
  return(  list(old_dirs=old_dirs,old_stps=old_stps,Hdiag=Hdiag)  )
}
