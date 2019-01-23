
###
library(Rcpp)
library(gaussquad)
library(combinat)
source("../functions/RLS_functions.R")
source("../functions/MG_functions.R")
sourceCpp("../functions/Rcpp_functions.cpp")

### data
P <- 3
data <- read.csv( paste("../data/data_",P,"D.csv",sep=""), dec="," )
X <- data[,1:P]
y <- data$y
n <- length(y)

### parameter
G <- 5
m <- lapply( 1:G, function(g) rep(2^g-1,P) )
q <- rep(3,P)
J <- lapply(1:G, function(g) m[[g]]+q+1)
K <- prod(J[[G]])
Omega <- lapply( 1:P, function(p) c(min(X[,p]),max(X[,p])) )

### basis matrices Phi_p^T for g=1,...,G
Phi_t_list <- lapply( 1:G, function(g) lapply( 1:P, function(p) t(my_bs_matrix(X[,p],m[[g]][p],q[p],Omega[[p]])) ) )

### thin plate penalty matrices Psi_rp^p for g=1,...,G
Psi_list <- lapply( 1:G, function(g) my_TP_regularization(m[[g]],q,Omega) )

### right-hand side vector
b <- MVP_kr_Rcpp(Phi_t_list[[G]],y)
norm_b <- sum(sqrt(b^2))

### smoothing parameter (manually)
lambda <- 0.2

### a priori weighting of penalty matrices
nreg <- length(Psi_list[[G]])
for(g in 1:G){
  for(j in 1:nreg){
    Psi_list[[g]][[j]][[1]] <- lambda*Psi_list[[g]][[j]][[1]]
  }
}

### setup for matrix-free Jacobi-method (digonal and weight)
diags <- lapply( 1:G, function(g) diag_kr_Rcpp_fct(Phi_t_list[[g]]) + rowSums(sapply( 1:nreg, function(j) diag_kron_Rcpp_fct(Psi_list[[g]][[j]]) ))  )
if(P==4){
  w <- 0.1
} else{
  w <- 0.2
}
nu <- c(6,3)

### grid transfer matrices for g=2,...,G
Prol <- lapply( 2:G, function(g) lapply( 1:P, function(p) my_prolongation(J[[g-1]][p],J[[g]][p],q[p]) ) )
Rest <- lapply( 2:G, function(g) lapply( 1:P, function(p) my_restriction(J[[g]][p],J[[g-1]][p],q[p]) ) )

### matrix-free MGCG
alpha <- rep(0,K)
max_iter <- 20
tol <- 10^(-5)
for(k in 1:max_iter){
  
  alpha <- my_vcycle( Phi_t_list, Psi_list, diags, Rest, Prol, b, nu, w, alpha)
  
  base <- MVP_kr_Rcpp( Phi_t_list[[G]], MVP_krtrans_Rcpp(Phi_t_list[[G]],alpha) )
  pen <- lapply( 1:nreg, function(j) MVP_kron_Rcpp(Psi_list[[G]][[j]],alpha) )
  res <- b- (base + Reduce("+",pen))
  cat(sum(sqrt(res^2))/norm_b ,"\n")
  if( sum(sqrt(res^2))/norm_b < tol ){
    break
  }
  
}



