###---- Test Data ----------------------------------------------------

T <- 12

data <- data.frame(i = rep(1,T), t = c(1:T), 
                   y = c(1, 2, 4, 6, 4, 12, 11, 13,11, 14, 16, 17))

## calcuating first lags of dependant and independant variables
library(dplyr)
data <- data %>% 
  dplyr::mutate(ilag1 = lag(i), tlag1 = lag(t), ylag1 = lag(y))

###---- OLS - ordinary least square estimaion ------------------------

ols <- function(depvar, indvar, res = NULL){
  y <- as.matrix(depvar)
  x <- cbind(1, indvar)

  beta <- (solve(t(x) %*% x)) %*% (t(x) %*% y) #compute ols estimates
  
  # st. error of beta_ols
  res_ols  <- data.frame(res = y - x %*% beta) 
  library(dplyr)
  res_ols <- res_ols %>% 
    dplyr::mutate(res_lag1 = lag(res))
  
  sigma_e2_ols <- sum(res_ols$res ^ 2)/(nrow(y) - ncol(x)) # df = n - k

  var_cov_ols <- sigma_e2_ols * (solve(t(x) %*% x))
  beta <- data.frame(parameter = c("Intercept", "Slope"), 
                     coefficient = beta, 
                     st.error = c(sqrt(var_cov_ols[1,1]), 
                                  sqrt(var_cov_ols[2,2]))) 
  
  # T-statistic for null hypothesis beta = 0
  beta <- beta %>%
    dplyr::mutate(t_stat = (coefficient/st.error))

  # Compile results for output

  ifelse(res == FALSE, 
         result <- list(estimates = beta, 
                        sigma_e2_ols = sigma_e2_ols),
         result <- list(estimates = beta, 
                        sigma_e2_ols = sigma_e2_ols, 
                        residuals_ols = res_ols))
  return(result)
}


ols_est <- ols(depvar = data$y, indvar = data$t, res = TRUE)

# estimating the coefficient of AR(1)
phi1_ols <- ols(depvar = ols_est$residuals_ols[-1,1], 
                     indvar = ols_est$residuals_ols[-1,2], res = FALSE)

# updating the list to include required output only
ols_est <- list(beta_ols = ols_est$estimates, 
                phi1_ols = phi1_ols$estimates[2,], 
                sigma_e2_ols = ols_est$sigma_e2_ols)
rm(phi1_ols)


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
  y <- data$y
  ylag1 <- data$ylag1[-1]
  x <- cbind(data$i, data$t)
  xlag <- cbind(data$ilag1[-1] ,data$tlag1[-1])
  
  phi1_nls <- theta[1]
  beta_nls <- theta[2:3]

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


nls <- function(nls_func, foc_nls, data, ols_est){
  
  y <- as.matrix(data$y)
  ylag1 <- data$ylag1[-1]
  x <- cbind(data$i, data$t)
  xlag <- cbind(data$ilag1[-1] ,data$tlag1[-1])
  
  phi1_nls <- ols_est$phi1_ols$coefficient
  beta_nls <- ols_est$beta_ols$coefficient
  
  nls_par <- optim(par = c(phi1_nls, beta_nls), 
                   fn = nls_func, gr = foc_nls, data = data, 
                   method = "BFGS")
  
  res_nls <- y[-1] - x[-1,] %*% beta_nls 
  sigma_e2_nls <- sum((res_nls ^ 2) / (T - 4))
  
  var_cov_nls <- sigma_e2_nls * 
    solve(
      rbind(cbind(t(x[-1,] - nls_par$par[1] * xlag) %*% 
                    (x[-1,] - nls_par$par[1] * xlag), 
                  t(x[-1,] - nls_par$par[1] * xlag) %*% 
                    (ylag1 - xlag %*% nls_par$par[2:3])), 
            cbind(t(t(x[-1,] - nls_par$par[1] * xlag) %*% 
                      (ylag1 - xlag %*% nls_par$par[2:3])), 
                  t((ylag1 - xlag %*% nls_par$par[2:3])) %*% 
                    (ylag1 - xlag %*% nls_par$par[2:3])))
    )
  
  phi1_nls <- data.frame(parameter = "Slope", estimate = nls_par$par[1], 
                         st.error = sqrt(var_cov_nls[3,3]), 
                         t_stat = nls_par$par[1] / 
                           sqrt(var_cov_nls[3,3]))
  
  beta_nls <- data.frame(parameter = c("Intercept", "Slope"),
                           estimate = nls_par$par[2:3], 
                         st.error = c(sqrt(var_cov_nls[1,1]), 
                                      sqrt(var_cov_nls[2,2])), 
                         t_stat = nls_par$par[2:3] / 
                           c(sqrt(var_cov_nls[1,1]), 
                             sqrt(var_cov_nls[2,2])))
  
  return(list(beta_nls = beta_nls, phi1_nls = phi1_nls, 
                  sigma_e2_nls = sigma_e2_nls))
}  

nls_est <- nls(nls_func, foc_nls, data, ols_est)


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
  phi1_mle <- (exp(theta[1]) - 1) / (exp(theta[1]) + 1)
  beta_mle <- theta[2:3]
  sigma_e2_mle <- (exp(theta[4]/2)) ^ 2 #transformed variance of epsilon
  
  # loglikelihood function of y with AR(1) errors and normal iid 
  # innovations epsilon (negative of loglikelihood for maximisation)
  u <- y - x %*% beta_mle
  ulag1 <- ylag1 - xlag %*% beta_mle
    
  - (- T/2*log(2 * pi) - T/2 * log(sigma_e2_mle) + 
       1/2 * log(1 - phi1_mle ^ 2) - 1/(2 * sigma_e2_mle) * (
         (1 - phi1_mle ^ 2) * (u[1]) ^ 2 + sum(
           (u[-1] - phi1_mle * (ulag1)) ^ 2)
       )
  )
}

mle <- function(loglikelihood, score, data, ols_est){
    
  ## Use OLS for initial conditions
  theta <- c(ols_est$phi1_ols$coefficient, ols_est$beta_ols$coefficient, 
           ols_est$sigma_e2_ols)

  ## Optimisation [using ols estimates as initial conditions]
  mle <- optim(par = theta, 
               fn = loglikelihood, data = data, 
               method = "BFGS", hessian = TRUE)
  
  r_theta <- matrix(
    c(2 * exp(mle$par[1])/(exp(mle$par[1]) + 1) ^ 2, rep(0,4), 1, 
      rep(0,4), 1, rep(0,4), exp(mle$par[4])), 4, 4
  )
  
  # covariance matrix using delta method
  var_mle <- r_theta %*% solve(mle$hessian)  %*% t(r_theta)
  
  #Storing the results and transforming phi1 and sigma_e2 from theta
  
  result = list(beta_mle = data.frame(
    parameter = c("Intercept","Slope"), coefficient = mle$par[2:3], 
    st.error = c(sqrt(var_mle[2,2]),sqrt(var_mle[3,3])), 
    t_stat = c(mle$par[2]/sqrt(var_mle[2,2]), 
               mle$par[3]/sqrt(var_mle[3,3]))), 
    phi1_mle = data.frame(parameter = "Slope", 
                          coefficient = (exp(mle$par[1]) - 1)/
                            (exp(mle$par[1]) + 1), 
                          st.error = sqrt(var_mle[1,1]), 
                          t_stat = mle$par[1]/sqrt(var_mle[1,1])), 
    sigma_e2_mle = (exp(mle$par[4]/2)) ^ 2)
  
  return(result)
  
}

mle_est <- mle(loglikelihood, score, data, ols_est)

