library(dplyr)

###-----Data Generating Process (DGP)---------------------------------

## Data Generating Process (DGP)
aplha = 0
beta = 0
phi2 = 0
mean_e = 0
sigma_e =  3

y0 = 0 #only if phi1=1
phi <- c(0.5, 0.9, 0.99, 1) #use phi1=phi for each cases

T = 10

## Replication of data
data <- data.frame(epsilon = rnorm(T, mean = mean_e, sd = sigma_e), 
                   i = rep(1,T), t = 1:T)

y <- NULL

for (t in 1:T) {
  ifelse(t == 1,
         ifelse(abs(phi[1]) < 1, 
                y[t] <- rnorm(1, mean = 0, 
                              sd = sqrt((sigma_e^2)/(1-phi[1]^2))),
                y[t] <- c(phi[1] * y0 + data$epsilon[t])),
         y[t] <- c(phi[1] * y[as.numeric(t - 1)] + data$epsilon[t]))
}

data <- data.frame(data, y)

## calcuating first lags of dependant and independant variables

data <- data %>% 
  dplyr::mutate(ilag1 = lag(i), tlag1 = lag(t), ylag1 = lag(y))

###---- OLS - ordinary least square estimaion ------------------------

x <- cbind(data$i, data$t)
beta_ols <- (solve(t(x) %*% x)) %*% (t(x) %*% y) #compute ols estimates
rownames(beta_ols) <- c("intercept", "t")

res_ols  <- data.frame(res_ols = y - x %*% beta_ols) 
library(dplyr)
res_ols <- res_ols %>% 
  dplyr::mutate(res_ols_lag1 = lag(res_ols))
phi1_ols <- (solve(t(res_ols_lag1[-1]) %*% res_ols_lag1[-1])) %*% 
  (t(res_ols_lag1[-1]) %*% res_ols$res_ols[-1]) #compute ols estimates

sigma_e2_ols <- sum(res_ols$res_ols ^ 2) / (T - 2) #st. error of OLS 

###---- MLE - maximum likelihood estimation --------------------------

## Define Log-likelihood function

# Arguments: theta - vector of parameters as specied below
#               theta = (phi1_mle, alpha_mle, beta_mle, sigma_e_mle)

#            data - data.frame with dependant and independant vars 
#                 and their first lags (see DGP above)

loglikelihood <- function(theta, data){ 
  
  y <- data$y # vector of regressand/dependant variable
  ylag1 <- data$ylag1[-1] # first lag of y
  x <- cbind(data$i, data$t) # regressors i = column of 1s, t = trend variable
  xlag <- cbind(data$ilag1[-1], data$tlag1[-1]) # first lag of regressors
  
  # parameters
  phi1_mle <- theta[1]
  beta_mle <- theta[2:3]
  sigma_e2_mle <- theta[4] ^ 2 #variance of epsilon
  
  # loglikelihood function of y with AR(1) errors and normal iid 
  # innovations epsilon (negative of loglikelihood for maximisation)
  
  - (- T/2*log(2 * pi) - T/2 * log(sigma_e2_mle) + 
    1/2 * log(1 - phi1_mle ^ 2) - 1/(2 * sigma_e2_mle) * (
      (1 - phi1_mle ^ 2) * (y[1] - x[1,] %*% beta_mle) ^ 2 + sum(
        (y[-1] - x[-1,] %*% beta_mle - phi1_mle * 
           (ylag1 - xlag %*% beta_mle)) ^ 2
        )
      )
    )
}

## Gradient (first derivative of the loglikelihood function)

score <- function(theta, data){
  y <- data$y
  ylag1 <- data$ylag1[-1]
  x <- cbind(data$i, data$t)
  xlag <- cbind(data$ilag1[-1] ,data$tlag1[-1])
  phi1_mle <- theta[1]
  beta_mle <- theta[2:3]
  sigma_e2_mle <- theta[4] ^ 2 #variance of epsilon
  
  c( 
    # FOC wrt phi1
    - phi1_mle/(1 - phi1_mle ^ 2) + 1/sigma_e2_mle * (
      phi1_mle * (y[1] - x[1,] %*% beta_mle) ^ 2 + sum(
        (ylag1 - xlag %*% beta_mle) * 
          (y[-1] - x[-1,] %*% beta_mle - phi1_mle * 
             (ylag1 - xlag %*% beta_mle)
          )
      )
    )
    ,
    
    # FOC wrt beta
    1/sigma_e2_mle * (
      (1 - phi1_mle ^ 2) * x[1,] %*% (y[1] - x[1,] %*% beta_mle) +
        sum(
          t(x[-1,] - phi1_mle * xlag) %*% 
            (y[-1] - x[-1,] %*% beta_mle - 
               phi1_mle * (ylag1 - xlag %*% beta_mle))
        )
    )
    ,
    
    # FOC wrt sigma_e2
    - T/(2 * sigma_e2_mle) + 1/(2 * sigma_e2_mle ^ 2) * (
      (1 - phi1_mle ^ 2) * (y[1] - x[1,] %*% beta_mle) ^ 2 + sum(
        (y[-1] - x[-1,] %*% beta_mle - 
           phi1_mle * (ylag1 - xlag %*% beta_mle)) ^ 2
      )
    )
  )
}

## Use OLS for initial conditions

## Optimisation [using ols estimates as initial conditions]
mle <- optim(par = c(phi1_ols, beta_ols, sigma_e2_ols), 
             fn = loglikelihood, gr = score, data = data, 
             method = "BFGS", hessian = TRUE)

as.matrix(mle$par)
var_mle <- (- mle$hessian)

# Different initial values lead to different estimates and variance is 
# either too underestimated or hessian has *negative variance

###---- NLS - nonlinear least square estimation ----------------------

nls_func <- function(theta, data){
  
  y <- data$y
  ylag1 <- data$ylag1[-1]
  x <- cbind(data$i, data$t)
  xlag <- cbind(data$ilag1[-1] ,data$tlag1[-1])
  
  phi1_nls <- theta[1]
  beta_nls <- theta[2:3]
#  sigma_e2_nls <- theta[4] ^ 2 #variance of epsilon
  
  # nls objective function [epsilon*epsilon]
  sum(
    ((y[-1] - x[-1,] %*% beta_nls) - phi1_nls * (ylag1 - xlag %*% beta_nls)) ^ 2
  )
}

foc_nls <- function(theta, data){
  c(
    # FOC wrt phi1
    - 2 * t((ylag1 - xlag %*% beta_nls)) %*% 
      ((y[-1] - x[-1,] %*% beta_nls) - phi1_nls * 
         (ylag1 - xlag %*% beta_nls))
    ,
    # FOC wrt beta (intercept and slope)
    2 * t((-x[-1,] + phi1_nls * xlag)) %*% 
      ((y[-1] - x[-1,] %*% beta_nls) - phi1_nls * 
         (ylag1 - xlag %*% beta_nls))
  )
}

nls_est <- optim(par = c(0.4, 0, 0), fn = nls_func, gr = foc_nls, 
             data = data, method = "BFGS")

res_nls <- y - x %*% nls$est[2:3] 
sigma_e2_nls <- sum((res_nls ^ 2) / (T - 4))
nls_est$par
sigma_e2_nls