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


###---- MLE - maximum likelihood estimation --------------------------

## Define Log-likelihood function

# Arguments: theta - vector of parameters as specied below
#               theta = (phi1_mle, alpha_mle, beta_mle, sigma_e_mle)

#            data - data.frame with dependant and independant vars 
#                 and their first lags (see DGP above)

loglikelihood <- function(theta, data){ 
  
  y <- data$y # vector of regresand/dependant variable
  ylag1 <- data$ylag1[-1] # first lag of y
  x <- cbind(data$i, data$t) # regressors i = column of 1s, t = trend variable
  xlag <- cbind(data$ilag1[-1], data$tlag1[-1]) # first lag of regressors
  
  # parameters
  phi1_mle <- theta[1]
  beta_mle <- theta[2:3]
  sigma_e2_mle <- theta[4] ^ 2 #variance of epsilon
  
  # loglikelihood function of y with normal iid innovations epsilon
  
  (- T/2*log(2 * pi) - T/2 * log(sigma_e2_mle) + 
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

## initial values for optimisation algorithm
theta = as.matrix(data.frame(phi1_mle = 0.5, alpha_mle = 0, 
                             beta_mle = 0, sigma_e_mle = 3))

## Optimisation [not working]
mle <- optim(par = c(0.9, 0, 0, 8), fn = loglikelihood, gr = score, 
             data = data, method = "BFGS", hessian = TRUE, control = 
               list(fnscale = -1))


