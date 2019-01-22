

### 
set.seed(111)
P <- 2
n <- 100000

### Data
X <- sapply(1:P, function(p) runif(n,0,1))
e <- rnorm(n,0,0.1)
t <- sapply( 1:n, function(i) -16*( (sum(X[i,]^2) /length(X[i,])) -0.5) )
fx <- 1 / ( 1 + exp(t) )
y <- fx + e

###
data <- cbind(X,y,fx)
colnames(data) <- c( sapply(1:P, function(p) paste("x_",p,sep="")), "y", "fx" )
write.csv(data, file=paste("data_",P,"D.csv",sep="") )

