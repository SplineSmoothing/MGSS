# MGSS: A Matrix-Free Multigrid Preconditioner for Spline Smoothing
Description: Data smoothing with penalized splines is a popular method and is well established for one- or two-dimensional covariates. The extension to multiple covariates is straightforward but suffers from exponentially increasing memory requirements and computational complexity. This toolbox provides a matrix-free implementation of a conjugate gradient (CG) method for the regularized least squares problem resulting from tensor product B-spline smoothing with multivariate and scattered data. It further provides matrix-free preconditioned versions of the CG-algorithm where the user can choose between a simpler diagonal preconditioner and an advanced geometric multigrid preconditioner. The main advantage is that all algorithms are performed matrix-free and therefore require only a small amount of memory. For further detail see Siebenborn & Wagner (2019) <arXiv:1901.00654>.

## Demo
The demo code allows to choose between a pure CG and a multigrid preconditioned CG solver.
A train and test data set with `n` samples in `P` dimensions is generated via the function `data <- generate_test_data(n, P)`.

### Generate data (train + test)

```R
n <- 100000
P <- 3
data <- generate_test_data(n, P)
X_train <- data$X_train
y_train <- data$y_train
```
### Setup of spline and regularization parameters

```R
G <- 5
m <- rep(2^G-1, P)
q <- rep(3, P)
#penalty_type <- "curve"
penalty_type <- "diff"

lambda <- 0.1
```
### Solver selection

```R
#model <- CG_smooth(m, q, lambda, X_train, y_train, pen_type = penalty_type, tolerance = 0.1)
#model <- PCG_smooth(m, q, lambda, X_train, y_train, pen_type = penalty_type)
model <- MGCG_smooth(G, q, lambda, X_train, y_train, tolerance = 0.1)
```

#### Model performance on training data

```R
model$rmse
model$R_squared
n <- length(y_train)
res <- model$residuals
RSS <- sum(res^2)
df <- estimate_trace(m, q, lambda, X_train, pen_type = penalty_type)
AIC <- log(RSS) + 2*df/n
AIC_C <- log(RSS) + 2*(df+1) / (n-df-2)
```

#### New predictions on test data and model validation
```R
X_test <- data$X_test
y_test <- data$y_test

y_pred <- predict_smooth(model, X_test)
res <- y_pred - y_test
RSS <- sum(res^2)
mean_res <- mean(res)
sd_res <- sd(res)
rmse <- sqrt( mean( res^2 ) )
R_squared <- 1 - ( sum(res^2) / sum( (y_test-mean(y_test))^2 ) )
```
