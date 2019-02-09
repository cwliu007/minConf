minConf_TMP <- function(
  x, funObj, gr = NULL, he = NULL, LB = -Inf, UB = Inf, method = "lbfgs", ...
  ,verbose=0,numDiff=1,optTol=1e-6
  ,maxIter=500,suffDec=1e-4,interp=1,corrections=100,damped=0
  ){

  nVars = length(as.vector(x))

  if (length(LB)==1){
    LB = rep(LB,nVars)
    UB = rep(UB,nVars)
  }

  # Make objective function (if using numerical derivatives)
  funEvalMultiplier = 1;
  if (numDiff>0){
    if (numDiff == 2){
      useComplex = 1;
    }else{
      useComplex = 0;
    }
    # funObj = autoGrad(x,numDiff,funObj,...) 
    funEvalMultiplier = nVars+1-useComplex;
  }

  # Evaluate Initial Point
  x = projectBounds(x,LB,UB)

  f = funObj(x ,...)
  
  if (numDiff==0){
    g = gr(x ,...)
  }else{
    g = autoGrad(x,useComplex,numDiff,funObj,...) 
  }
  g = as.vector(g)

  if (method=="newton"){
    H = he(x ,...)
    secondOrder = TRUE
  }else{
    H = NULL
    secondOrder = FALSE
  }

  funEvals = 1

  # Compute Working Set
  working = rep(1,nVars);
  working[(x < LB+optTol*2) & g >= 0] = 0;
  working[(x > UB-optTol*2) & g <= 0] = 0;
  working = which(working != 0)

  # Check Optimality
  if (length(working)==0){
    message = "All variables are at their bound and no further progress is possible at initial point"
    if (verbose >= 1) print(message)

    return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )

  }else if( norm(g[working],"2") <= optTol ){
    message = "All working variables satisfy optimality condition at initial point"
    if (verbose >= 1) print(message)

    return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )

  }

  if (verbose >= 3){
    switch(method,
           sd = {
             print("Steepest Descent")
           },
           lbfgs = {
             print("L-BFGS")
           },
           bfgs = {
             print("BFGS")
           },
           newton = {
             print("Newton")
           },
           {
             message = paste("Unrecognized Method: ",method,sep="")
             return(message)
           }
          )
  }

  
  i = 1
  while (funEvals <= maxIter) {

    # Compute Step Direction
    d = rep(0,nVars)
    switch(method,
           sd = {
             d[working] = - g[working]
           },
           lbfgs = {
             if (i==1){
               d[working] = - g[working]
               old_dirs = matrix(0,nVars,0)
               old_stps = matrix(0,nVars,0)
               Hdiag = 1
             }else{
               if (damped!=0){
                 out = dampedUpdate(g-g_old,x-x_old,corrections,verbose==3,old_dirs,old_stps,Hdiag)
               }else{
                 out = lbfgsUpdate(g-g_old,x-x_old,corrections,verbose==3,old_dirs,old_stps,Hdiag)
               }
               old_dirs = out$old_dirs
               old_stps = out$old_stps
               Hdiag = out$Hdiag
               
               curvSat = which(colSums(as.matrix(old_dirs[working,])*as.matrix(old_stps[working,])) > 1e-10)
               d[working] = lbfgs(-g[working],as.matrix(old_dirs[working,curvSat]),as.matrix(old_stps[working,curvSat]),Hdiag)
             }
             g_old = g
             x_old = x
           },
           bfgs = {
             if (i == 1){
               d[working] = - g[working]
               B = diag(nVars)
             }else{
               y = g-g_old
               s = x-x_old

               ys = t(y)%*%s

               if (i == 2){
                 if (ys > 1e-10){
                   B = as.numeric(((t(y)%*%y)/(t(y)%*%s)))*diag(nVars)
                 }
               }
               if (ys > 1e-10){
                 B = B + (y%*%t(y))/as.numeric(t(y)%*%s) - (B%*%s%*%t(s)%*%B)/as.numeric(t(s)%*%B%*%s)
               }else{
                 if (verbose==2) print("Skipping Update")
               }
               d[working] = - solve(B[working,working],g[working])
             }
             g_old = g
             x_old = x
           },
           newton = {
             
             posDef = all( eigen(H[working,working],symmetric=TRUE,only.values=TRUE)$values > 0 )
             
             if (posDef){
               R = chol(H[working,working])
               d[working] = - solve(R,solve(t(R),g[working]))
             }else{
               if (verbose==3) print("Adjusting Hessian")
               H[working,working] = H[working,working] + diag(length(working)) * max(0,1e-12 - min(Re(eigen(H[working,working])$values)))
               d[working] = - solve(H[working,working],g[working])
             }
           }
           )

    # Check that Progress can be made along the direction
    f_old = f
    gtd = g%*%d
    if (gtd > -optTol){
      message = "Directional Derivative below optTol"
      if (verbose==3) print(message)
      return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )
    }


    # Select Initial Guess to step length
    if (i == 1 && !secondOrder){
      t = min(1,1/sum(abs(g[working])))
    }else{
      t = 1
    }

    # Evaluate the Objective and Projected Gradient at the Initial Step
    x_new = projectBounds(x+t*d,LB,UB)
    if (secondOrder){
      f_new = funObj(x_new ,...)
      if (numDiff==0){
        g_new = gr(x_new ,...)
        H = he(x_new ,...)
      }else{
        g_new = autoGrad(x_new,useComplex,numDiff,funObj,...) 
        H = autoHessian(x_new,useComplex,numDiff,funObj,...) 
      }
      g_new = as.vector(g_new)
      
    }else{
      f_new = funObj(x_new ,...)
      if (numDiff==0){
        g_new = gr(x_new ,...)
      }else{
        g_new = autoGrad(x_new,useComplex,numDiff,funObj,...) 
      }
      g_new = as.vector(g_new)
    }
    funEvals = funEvals+1

    ## Backtracking Line Search
    lineSearchIters = 1
    while ( (f_new > (f + suffDec*t(g)%*%(x_new-x)))  ||  (!isLegal(f_new)) ){
      temp = t
      if (interp == 0 || !isLegal(f_new) || !isLegal(g_new)){
        if (verbose==3) print("Halving Step Size")
        t = .5*t
      }else{
        if (verbose==3) print("Cubic Backtracking")
        t = polyinterp(t(matrix(c(0, f, gtd, t, f_new, t(g_new)%*%d),3,2)))$minPos
      }

      # Adjust if change is too small
      if (t < temp*1e-3){
        if (verbose==3) print("Interpolated value too small, Adjusting")
        t = temp*1e-3
      }else if(t > temp*0.6){
        if (verbose==3) print("Interpolated value too large, Adjusting")
        t = temp*0.6
      }

      # Check whether step has become too small
      if (sum(abs(t*d)) < optTol){
        message = "Line Search failed"
        if (verbose==3) print(message)
        t = 0
        f_new = f
        g_new = g

        return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )
      }

      # Evaluate New Point
      x_new = projectBounds(x+t*d,LB,UB)
      f_new = funObj(x_new ,...)
      if (numDiff==0){
        g_new = gr(x_new ,...)
      }else{
        g_new = autoGrad(x_new,useComplex,numDiff,funObj,...) 
      }
      g_new = as.vector(g_new)
      funEvals = funEvals+1
      lineSearchIters = lineSearchIters+1
    }

    # Take Step
    x = x_new;
    f = f_new;
    g = g_new;

    # Compute Working Set
    working = rep(1,nVars)
    working[(x < LB+optTol*2) & (g >= 0)] = 0
    working[(x > UB-optTol*2) & (g <= 0)] = 0
    working = which(working != 0)

    # Output Log
    if (verbose>=2) print(c(as.integer(i),as.integer(funEvals*funEvalMultiplier),t,f,sum(abs(g[working]))))

    # Check Optimality
    if (length(working)==0){
      message = "All variables are at their bound and no further progress is possible at initial point"
      if (verbose>=1) print(message)

      return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )

    }else if( norm(g[working],"2") <= optTol ){
      message = "All working variables satisfy optimality condition at initial point"
      if (verbose>=1) print(message)

      return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )
    }

    # Check for lack of progress
    if (sum(abs(t*d)) < optTol){
      message = "Step size below optTol"
      if (verbose>=1) print(message)

      return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )
    }

    if (abs(f-f_old) < optTol){
      message = "Function value changing by less than optTol"
      if (verbose>=1) print(message)

      return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )
    }

    if (funEvals*funEvalMultiplier > maxIter){
      message = "Function Evaluations exceeds maxIter"
      if (verbose>=1) print(message)

      return(  list(par=x,objective=f,funEvals=funEvals,g=g,H=H,message=message)  )
    }

    # If necessary, compute Hessian
    if (secondOrder && (lineSearchIters > 1) ){
      f_new = funObj(x ,...)
      if (numDiff==0){
        g_new = gr(x ,...)
        H = he(x ,...)
      }else{
        g_new = autoGrad(x,useComplex,numDiff,funObj,...) 
        H = autoHessian(x,useComplex,numDiff,funObj,...)
      }
      g_new = as.vector(g_new)
    }

    i = i + 1
  }

}


projectBounds <- function(x,LB,UB){
  ind = x < LB
  x[ind] = LB[ind]
  ind = x > UB
  x[ind] = UB[ind]
  return(x)
}

isLegal <- function(v){
  legal = sum(any(Im(v[])!=0))==0 && sum(is.na(v[]))==0 && sum(is.infinite(v[]))==0
}
