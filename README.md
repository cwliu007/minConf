# minConf_TMP for minimization

R code of the minConf_TMP (converted from minConf package in MATLAB). (See usage at https://www.cs.ubc.ca/~schmidtm/Software/minConf.html)

---

How To Install in R:

library(devtools)

install_github("genwei007/minConf_TMP")

---

Example:

fw <- function(x) (1-x[1])^2 + 100*(x[2]-x[1]^2)^2 # Rosenbrock function where c(1,1) is minimum

x = c(0,0)

m1 <- optim(x, fw, method = "L-BFGS-B")

m1$par

#0.9998007, 0.9996013

m2 <- minConf_TMP(x, fw, method="lbfgs")

m2$par

#0.9999954, 0.9999901

m3 <- minConf_TMP(x, fw, method="newton")

m3$par

#0.9999458, 0.9998900
