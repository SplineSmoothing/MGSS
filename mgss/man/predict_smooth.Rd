% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mgss_main.R
\name{predict_smooth}
\alias{predict_smooth}
\title{Predictions from model}
\usage{
predict_smooth(model_smooth, X)
}
\arguments{
\item{model_smooth}{A spline model resulting from \code{CG_smooth}, \code{PCG_smooth}, or \code{MGCG_smooth}.}

\item{X}{Matrix containing the new observations.}
}
\value{
Vector of length \code{nrow(X)} of predictions.
}
\description{
Makes predictions of new observations from a fitted spline model.
}
\examples{
data <- generate_test_data(100, 2)
X <- data$X_train
y <- data$y_train
result <- PCG_smooth(m = c(7,7), q = c(3,3), lambda = 0.1, X = X, y = y, print_error = FALSE)
X_test <- data$X_test
predict_smooth(model_smooth = result, X = X_test)

}
