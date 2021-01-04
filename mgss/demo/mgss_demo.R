library(mgss)

#####--------------------------------------------
##### generate data (train + test)
n <- 100000
P <- 3
data <- generate_test_data(n, P)
X_train <- data$X_train
y_train <- data$y_train

#####--------------------------------------------
##### Setup of spline parameters
#P <- ncol(X_train)
G <- 5
m <- rep(2^G-1, P)
q <- rep(3, P)
#penalty_type <- "curve"
penalty_type <- "diff"

#####--------------------------------------------
##### regularization parameter
lambda <- 0.1


#####--------------------------------------------
##### Provided algorithms
#model <- CG_smooth(m, q, lambda, X_train, y_train, pen_type = penalty_type, tolerance = 0.1)
#model <- PCG_smooth(m, q, lambda, X_train, y_train, pen_type = penalty_type)
model <- MGCG_smooth(G, q, lambda, X_train, y_train, tolerance = 0.1)


#####--------------------------------------------
##### Model performance on training data
model$rmse
model$R_squared
n <- length(y_train)
res <- model$residuals
RSS <- sum(res^2)
df <- estimate_trace(m, q, lambda, X_train, pen_type = penalty_type)
AIC <- log(RSS) + 2*df/n
AIC_C <- log(RSS) + 2*(df+1) / (n-df-2)


#####--------------------------------------------
##### New predictions on test data and model validation
X_test <- data$X_test
y_test <- data$y_test

y_pred <- predict_smooth(model, X_test)
res <- y_pred - y_test
RSS <- sum(res^2)
mean_res <- mean(res)
sd_res <- sd(res)
rmse <- sqrt( mean( res^2 ) )
R_squared <- 1 - ( sum(res^2) / sum( (y_test-mean(y_test))^2 ) )

