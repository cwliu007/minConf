polyinterp <- function(points,doPlot,xminBound,xmaxBound){
  # function [minPos] = polyinterp(points,doPlot,xminBound,xmaxBound)
  # 
  # Minimum of interpolating polynomial based on function and derivative
  # values
  # 
  # In can also be used for extrapolation if {xmin,xmax} are outside
  # the domain of the points.
  # 
  # Input:
  #   points(pointNum,[x f g])
  # doPlot: set to 1 to plot, default: 0
  # xmin: min value that brackets minimum (default: min of points)
  # xmax: max value that brackets maximum (default: max of points)
  # 
  # set f or g to sqrt(-1) if they are not known
  # the order of the polynomial is the number of known f and g values minus 1
  
  doPlot = 1
  
  nPoints = dim(points)[1]
  order = sum((Im(points[,2:3])==0))-1
  
  # Code for most common case:
  # - cubic interpolation of 2 points
  #      w/ function and derivative values for both
  # - no xminBound/xmaxBound
  
  if (nPoints == 2 && order ==3 && doPlot == 0){
    minVal = min(points[,1])
    minPos = which.min(points[,1])
    notMinPos = -minPos+3
    d1 = points[minPos,3] + points[notMinPos,3] - 3*(points[minPos,2]-points[notMinPos,2]) %*% solve((points[minPos,1]-points[notMinPos,1]))
    d2 = sqrt(d1^2 - points[minPos,3]*points[notMinPos,3])
    if (is.double(d2)){
      t = points[notMinPos,1] - (points[notMinPos,1] - points[minPos,1])*((points[notMinPos,3] + d2 - d1)/(points[notMinPos,3] - points[minPos,3] + 2*d2))
      minPos = min(max(t,points[minPos,1]),points[notMinPos,1])
    }else{
      minPos = mean(points[,1])
    }
    return(list(minPos=minPos))
  }
  
  
  xmin = min(points[,1])
  xmax = max(points[,1])
  
  # Compute Bounds of Interpolation Area
  if (exists("xminBound")){
    xminBound = xmin
  }
  if (exists("xmaxBound")){
    xmaxBound = xmax
  }
  
  # Constraints Based on available Function Values
  A = matrix(0,0,order+1)
  b = matrix(0,0,1)
  for (i in 1:nPoints){
    if (Im(points[i,2])==0){
      constraint = matrix(0,1,order+1)
      for (j in order:0){
        constraint[order-j+1] = points[i,1]^j
      }
      A = rbind(A,constraint)
      b = rbind(b,points[i,2])
    }
  }
  
  # Constraints based on available Derivatives
  for (i in 1:nPoints){
    if (is.double(points[i,3])){
      constraint = matrix(0,1,order+1)
      for (j in 1:order){
        constraint[j] = (order-j+1)*points[i,1]^(order-j)
      }
      A = rbind(A,constraint)
      b = rbind(b,points[i,3])
    }
  }
  
  # Find interpolating polynomial
  params = solve(A,b)
  
  # Compute Critical Points
  dParams = matrix(0,order,1)
  for (i in 1:(length(params)-1)){
    dParams[i] = params[i]*(order-i+1)
  }
  
  if (any(is.infinite(dParams))){
    cp = c(xminBound,xmaxBound,points[,1])
  }else{
    vv = sort(Re(polyroot(rev(dParams))),decreasing=TRUE) # equal to "roots(dParams)" in MATLAB
    cp = c(xminBound,xmaxBound,points[,1],  vv  )
  }
  
  # Test Critical Points
  fmin = Inf
  minPos = (xminBound+xmaxBound)/2 # Default to Bisection if no critical points valid
  for (xCP in cp){
    if (Im(xCP)==0 && xCP>=xminBound && xCP<=xmaxBound){
      fCP = polyval_from_pracma(params,xCP)
      if (Im(fCP)==0 && fCP<fmin){
        minPos = Re(xCP)
        fmin = Re(fCP)
      }
    }
  }
  # Plot Situation
  if (doPlot==1){
    
    # to be added
  }
  
  return(list(minPos=minPos,fmin=fmin))
}