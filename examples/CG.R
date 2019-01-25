
###
#install.packages(c("Rcpp","gaussquad","combinat","polynom","orthopolynom"))
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
m <- rep(36,P)
q <- rep(3,P)
Omega <- lapply( 1:P, function(p) c(min(X[,p]),max(X[,p])) )
K <- prod(m+q+1)

### basis matrices Phi_p^T
Phi_t_list <- lapply(1:P, function(p) t(my_bs_matrix(X[,p],m[p],q[p],Omega[[p]])) )

### thin plate penalty matrices Psi_rp^p
Psi_list <- my_TP_regularization(m,q,Omega)
n_reg <- length(Psi_list)

### right-hand side vector
b <- MVP_kr_Rcpp(Phi_t_list,y)
norm_b <- sqrt(sum(b^2))

### smoothing parameter (manually)
lambda <- 0.2

### CG method
max_iter <- K
tol <- 10^(-5)
alpha <- rep(0,K)
res <- b
dir <- res
norm_res <- sqrt(sum(res^2))
for(i in 1:max_iter){
  
  base <- MVP_kr_Rcpp( Phi_t_list, MVP_krtrans_Rcpp(Phi_t_list,dir) )
  pen <- Reduce("+" , lapply( 1:n_reg, function(p) MVP_kron_Rcpp(Psi_list[[p]],dir) ) )
  Adir <- base + lambda*pen
  a <- as.numeric(norm_res^2 / crossprod(dir,Adir))
  alpha <- alpha + a*dir
  res <- res - a*Adir
  norm_res_old <- norm_res
  norm_res <- sqrt(sum(res^2))
  norm_res_rel <- norm_res/norm_b
  cat(norm_res_rel, "\n")
  if( norm_res_rel < tol){
    break
  }
  beta <- norm_res^2 / norm_res_old^2
  dir <- res + beta*dir
}




