# minConf_TMP
minConf_TMP in R for optimization

The code includes the R version of the minConf_TMP from minConf package in MATLAB. (See https://www.cs.ubc.ca/~schmidtm/Software/minConf.html)

Comments are welcome to improve the code.

--- EXAMPLE (also see https://www.cs.ubc.ca/~schmidtm/Software/TMP/examples.html#2)

1. Simulate data in MATLAB

rng(1)

% Set up Problem
nInst = 1000;
nVars = 100;
A = randn(nInst,nVars);
x = rand(nVars,1).*(rand(nVars,1) > .5);
b = A*x + randn;
funObj = @(x)SquaredError(x,A,b);

% Initial Value
x_init = zeros(nVars,1);

% Bounds
LB = zeros(nVars,1);
UB = inf(nVars,1);

save test.mat


2. Run minConf_TMP in R

funObj <- function(w,X,y){
  
  n = dim(X)[1]
  p = dim(X)[2]
  XX = t(X) %*% X
  
  if (n < p){
    Xw = X%*%w;
    res = Xw-y;
    f = sum(res^2);
    g = 2*(t(X)%*%res);
  }else{
    XXw = XX%*%w;
    Xy = t(X)%*%y;
    f = t(w)%*%XXw - 2*t(w)%*%Xy + t(y)%*%y;
    g = 2*XXw - 2*Xy;
  }

  H = 2*XX
  
  p = length(w)
  T = array(0,c(p,p,p))
  
  return(f)
}
gr <- function(w,X,y){
  
  n = dim(X)[1]
  p = dim(X)[2]
  XX = t(X) %*% X
  
  if (n < p){
    Xw = X%*%w;
    res = Xw-y;
    f = sum(res^2);
    g = 2*(t(X)%*%res);
  }else{
    XXw = XX%*%w;
    Xy = t(X)%*%y;
    f = t(w)%*%XXw - 2*t(w)%*%Xy + t(y)%*%y;
    g = 2*XXw - 2*Xy;
  }
  
  return(g)
}
he <- function(w,X,y){

  n = dim(X)[1]
  p = dim(X)[2]
  XX = t(X) %*% X

  
  H = 2*XX

  return(H)
}

library("R.matlab")

out = readMat("test.mat")

nInst = out$nInst
nVars = out$nVars
A = out$A
x = out$x
b = out$b

x_init = out$x.init

LB = out$LB
UB = out$UB


wTMP = minConf_TMP(x_init,funObj,gr=gr,he=he,LB=LB,UB=UB,method="lbfgs"
                   ,A,b,numDiff=0,verbose=3)

round(wTMP$par,4)
round(wTMP$g,4)
wTMP$objective
wTMP$funEvals
