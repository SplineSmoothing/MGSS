### Truncated Power Function
tp_fct = function(x,t,q){
  (x-t)^q * (x>t)
}

### B-spline Basis Matrix (equidistant knots)
my_bs_matrix = function(x,m,q,Omega=c(min(x),max(x))){
  
  h <- (Omega[2]-Omega[1]) / (m+1)
  knots <- seq(Omega[1]-q*h, Omega[2]+q*h, by = h)
  H <- outer(x, knots, tp_fct, q)
  D <- diff(diag(dim(H)[2]), diff = q+1) / (gamma(q+1)*h^q)
  Phi <- (-1)^(q+1) * tcrossprod(H,D)
  Phi[Phi<1e-10] <- 0
  return(Phi)
  
}

### Gramian Matrix of varphi_{-q,q}^d,...,varphi_{m,q}^d (equidistant knots)
my_L2norm_matrix = function(m,q,d,Omega){

  # Parameter
  K <- m+q+1
  h <- (Omega[2]-Omega[1]) / (m+1)
  
  # Gauss-Quadratur
  n_gauss <- q+1
  gauss_legendre <- legendre.quadrature.rules(n_gauss)[[n_gauss]]
  w_ref <- gauss_legendre$w[order(gauss_legendre$x)]
  points_ref <- (h/2)*sort(gauss_legendre$x) + (Omega[1]+(h/2))
  points <- as.vector( sapply( 1:(m+1), function(j) points_ref+(j-1)*h ) )
  w <- rep(w_ref,(m+1))
  Phi <- my_bs_matrix(points, m, (q-d), Omega)
  G <- (h/2)*crossprod(Phi, Phi*w)
  
  # Gramian Matrix
  if(d==0){
    Psi <- G
  } else{
    Delta <- diff(diag(K),diff=d)
    Psi <- (1 / (h^(2*d))) * crossprod(Delta,G%*%Delta)
  }
  return(Psi)
  
}

### List of Grmaian Matrices of each spatial dimension
my_L2norm_matrix_list = function(m,q,d,Omega){
  
  P <- length(q)
  Reg_list <- lapply(1:P, function(p) my_L2norm_matrix(m[p],q[p],d[p],Omega[[p]]) )
  return(Reg_list) 
  
}

### Thin Plate Regularization (D=2)
my_TP_regularization = function(m,q,Omega){
  
  P <- length(m)
  d_mat <- 2*diag(P)
  if(P>1){
    vec <- c(1,1,rep(0,P-2))
    M <- matrix( unlist( unique(permn(vec)) ), nrow=P, byrow=F )
    d_mat <- cbind(d_mat,M)
  }
  Psi_list <- lapply(1:dim(d_mat)[2], function(p) my_L2norm_matrix_list(m,q,d_mat[,p],Omega) )
  nreg <- length(Psi_list)
  w_Psi <- sapply(1:nreg, function(j) 2/prod(factorial(d_mat[,j])) )
  for(j in 1:nreg){
    Psi_list[[j]][[1]] <- w_Psi[j]*Psi_list[[j]][[1]]
  }
  return(Psi_list)
  
}

### Difference Penalty
my_diff_regularization = function(K_p,l,Omega){
  
  P <- length(K_p)
  Delta <- lapply(1:P, function(p) diff(diag(K_p[p]),diff=l[p]))
  Identity <- lapply(1:P, function(p) diag(K_p[p]))
  Delta_list <- lapply(1:P, function(p) Identity)
  for(p in 1:P){
    Delta_list[[p]][[p]] <- t(Delta[[p]])%*%Delta[[p]]
  }
  return(Delta_list)
  
}

