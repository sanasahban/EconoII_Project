
#Test Data
ytest <- data.frame(y = c(1,2,4,6,4,12,11,13,14,16,17))
ytest.mat <- as.matrix(ytest)

#OLS
#create matrix of regressors where i is column of 1s and t is trend variable
X <- data.frame(i = matrix(1,10,1),t=data$t) 
x.mat <- as.matrix(X) #convert to matrix
beta_ols <- (solve(t(x.mat)%*%x.mat))%*%(t(x.mat)%*%y) #compute ols estimates
rownames(beta_ols) <- c("intercept", "t")

l <- function(theta, data){ 
  -T/2*log(2*pi) - T*log(theta[4]) + 1/2*log(1 - theta[1]^2) -
    (1/(2*theta[4]^2))*((1 - theta[1]^2)*
                          (data$y[1] - cbind(data$i, data$t)[1,]%*%
                             theta[2:3])^2 + 
                          sum((data$y[-1] - cbind(data$i, data$t)[-1,]
                               %*%theta[2:3] - theta[1]*
                                 (data$ylag1[-1] - cbind(data$ilag1[-1], 
                                                         data$tlag1[-1])%*%theta[2:3]))^2))
}


# FOC wrt PHI from book
phi1_mle * (y[1] - x[1,] %*% beta_mle) ^ 2 - 
  (phi1_mle * sigma_e2_mle)/(1 - phi1_mle ^ 2) +
  sum(
    (ylag1 - xlag %*% beta_mle) * ( y[-1]-x[-1,] %*% beta_mle - 
                                      phi1_mle * (ylag1 - xlag %*% beta_mle))
  )

# FOC wrt beta from book
(1 - phi1_mle ^ 2) * x[1,] %*% (y[1] - x[1,] %*% beta_mle) + 
  sum(
    t( x[-1,] - phi1_mle * xlag ) %*% 
      ( y[-1]-x[-1,] %*% beta_mle - 
          phi1_mle * (ylag1 - xlag %*% beta_mle) 
      )
  )
