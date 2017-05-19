###---- Test Data ----------------------------------------------------

T <- 12

data <- data.frame(i = rep(1,T), t = c(1:T), 
                   y = c(1, 2, 4, 6, 4, 12, 11, 13,11, 14, 16, 17))


## calcuating first lags of dependant and independant variables

data <- data %>% 
  dplyr::mutate(ilag1 = lag(i), tlag1 = lag(t), ylag1 = lag(y))

###---- OLS - ordinary least square estimaion ------------------------

ols_func <- function(depvar, indvar, res = NULL){
  y <- as.matrix(depvar)
  x <- cbind(1, indvar)
  
  beta <- (solve(t(x) %*% x)) %*% (t(x) %*% y) #compute ols estimates
  
  # st. error of beta_ols
  res_ols  <- data.frame(res = y - x %*% beta) 
  library(dplyr)
  res_ols <- res_ols %>% 
    dplyr::mutate(res_lag1 = lag(res))
  
  sigma_e2_ols <- sum(res_ols$res ^ 2)/(nrow(y) - ncol(x)) # df = n - k

  var_cov <- sigma_sq_ols * (solve(t(x) %*% x))
  beta <- data.frame(parameter = c("intercept", "slope"), 
                     coefficient = beta, 
                     st.error = c(sqrt(var_cov[1,1]), 
                                  sqrt(var_cov[2,2]))) 
  
  # T-statistic for null hypothesis beta = 0
  beta <- beta %>%
    dplyr::mutate(t_stat = (coefficient/st.error))

  # Compile results for output

  ifelse(res == FALSE, 
         result <- list(estimates = beta, 
                        sigma_e2_ols = as.data.frame(sigma_e2_ols)),
         result <- list(estimates = beta, 
                        sigma_e2_ols = as.data.frame(sigma_e2_ols), 
                        residuals_ols = res_ols))
  return(result)
  
}

beta_ols <- ols_func(depvar = data$y, indvar = data$t)
phi1_ols <- ols_func(depvar = result$residuals_ols[-1,1], 
                     indvar = result$residuals_ols[-1,2], res = FALSE)


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
