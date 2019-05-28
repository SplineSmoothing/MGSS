#####--------------------------------------------
##### preamble

rm(list=ls())

library(gaussquad)                   # for curvature penalty
library(combinat)                    # for curvature penalty
library(Rcpp)                        # incorporate C++

source("./spline_functions.R")       # external R functions for P-splines
source("./multigrid_functions.R")    # external R functions for multigrid method
sourceCpp("./Rcpp_functions.cpp")    # external C++ functions

#####--------------------------------------------
##### simple test data

set.seed(123)
P <- 2                                          # number of covariates
n <- 10000                                      # number of observations
X <- sapply(1:P, function(p) runif(n,0,1))      # covariates
t <- sapply( 1:n, function(i) -16*( (sum(X[i,]^2) / length(X[i,])) -0.5) )
fx <- 1 / ( 1 + exp(t) )
y <- fx + rnorm(n, 0, 0.1)                      # observations

#####--------------------------------------------
##### spline system

G <- 5                                          # number of grids for multigrid
m <- lapply( 1:G, function(g) rep(2^g-1,P) )    # number of spline knots
q <- rep(3, P)                                  # spline degree
Omega <- lapply(1:P, function(p) c(0,1) )       # underlying space
J <- lapply(1:G, function(g) m[[g]]+q+1)        # number of directionla basis functions
K <- prod(J[[G]])                               # total number of basis functions
tPhi_list <- lapply(1:G, function(g) lapply(1:P, function(p) t( bspline_matrix(X[,p], m[[g]][p], q[p] ,Omega[[p]]) ) ) )    # spline matrices
Psi_list <- lapply(1:G, function(g)  curvature_penalty(m[[g]], q, Omega) )   # survature penalty
b <- MVP_khatrirao_rcpp(tPhi_list[[G]], y)      # right-hand side vector
norm_b <- sqrt( sum(b^2) )
lambda <- 10                                    # weight of the regularization (manually)

#####--------------------------------------------
##### multigrid setup

Prol <- lapply( 2:G, function(g) lapply( 1:P, function(p) prolongation_matrix(J[[g-1]][p],J[[g]][p],q[p]) ) )     # prolongation matrices
Rest <- lapply( 2:G, function(g) lapply( 1:P, function(p) restriction_matrix(J[[g]][p],J[[g-1]][p],q[p]) ) )      # restriction matrices
w <- 0.5         # damping parameter for (damped) Jacobi smoother
nu <- c(6,3)     # number of pre- and post-smoothing iterations

#####--------------------------------------------
##### matrix-free MGCG method

tol <- 10^(-6)          # tolerance for stopping criterion
alpha <- rep(0, K)      # starting value
if( max(alpha^2)!=0 ){
  Aalpha <- MVP_spline(tPhi_list[[G]], alpha) + lambda*MVP_penalty(Psi_list[[G]], alpha)
  r <- b-Aalpha
} else{
  r <- b
}
d <- r
z <- v_cycle(tPhi_list, Psi_list, Rest, Prol, lambda, b, nu, w, alpha)    # apply MG as preconditioner
d <- z
rz <- as.numeric( crossprod(r,z) )
cat(0, sqrt(sum(r^2))/norm_b, "\n")
for(i in 1:K){         # loop of the MGCG iteration
  
  Ad <- MVP_spline(tPhi_list[[G]], d) + lambda*MVP_penalty(Psi_list[[G]], d)
  t <- as.numeric( rz / crossprod(d, Ad) )
  alpha <- alpha+t*d
  rz_old <- rz
  r <- r-t*Ad
  cat(i, sqrt(sum(r^2))/norm_b, "\n")
  if(sqrt(sum(r^2))/norm_b <= tol){
    break
  }
  z <- v_cycle(tPhi_list, Psi_list, Rest, Prol, lambda, r, nu, w)     # apply MG as preconditioner
  rz <- crossprod(r,z)
  beta <-  as.numeric( rz / rz_old )
  d <- z + beta*d
  
}


