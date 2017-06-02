#install.packages("numDeriv")

library(dplyr)

###-----Data Generating Process (DGP)---------------------------------

aplha = 0
beta = 0
phi2 = 0
mean_e = 0
sigma_e =  3

y0 = 0 #only if phi1=1
phi <- c(0.5, 0.9, 0.99, 1) #use phi1=phi for each cases

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


mle <- function(loglikelihood, data, ols_est){
  
  ## Use OLS for initial conditions
  theta <- c(ols_est$phi1_ols$coefficient, ols_est$beta_ols$coefficient, 
             ols_est$sigma_e2_ols)
  
  ## Optimisation [using ols estimates as initial conditions]
  mle <- optim(par = theta, 
               fn = loglikelihood, data = data, 
               method = "BFGS")
  
  r_theta <- matrix(
    c(2 * exp(mle$par[1])/(exp(mle$par[1]) + 1) ^ 2, rep(0,4), 1, 
    rep(0,4), 1, rep(0,4), exp(mle$par[4])), 4, 4
  )
  
  library(numDeriv)
  hessian <- hessian(func = loglikelihood, x = mle$par, data = data)
  
  # covariance matrix using delta method
  var_mle <- r_theta %*% solve(hessian)  %*% t(r_theta)
  
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

###----SIMULATION-----------------------------------------------------

# R = No. of replications and T = sample size

simulation <- function(R, T){

  phi <- c(0.5, 0.9, 0.99, 1) #use phi1=phi for each cases
  
  # declaring lists for output
  null <- data.frame(phi1 = rep(0,R), st.error = rep(0,R), 
                     t_stat = rep(0,R))
  
  phi1_ols_sim <- list(phi1=null, phi2=null, phi3=null, phi4=null)
  names(phi1_ols_sim) <- paste0("phi_", phi[1:4])
  
  alpha_mle_sim <- list(phi1=null, phi2=null, phi3=null, phi4=null)
  names(alpha_mle_sim) <- paste0("phi_", phi[1:4])
  
  beta_mle_sim <- list(phi1=null, phi2=null, phi3=null, phi4=null)
  names(beta_mle_sim) <- paste0("phi_", phi[1:4])
  
  phi1_mle_sim <- list(phi1=null, phi2=null, phi3=null, phi4=null)
  names(phi1_mle_sim) <- paste0("phi_", phi[1:4])
  
  # Simulation
  
  for (i in 1:4){
    
    phi1 <- phi[i] # defining value of phi1
  
    for (r in 1:R){
      
      ## Replications of data
      data <- data.frame(epsilon = rnorm(T, mean = mean_e, 
                                         sd = sigma_e), 
                         i = rep(1,T), t = 1:T)
      
      y <- NULL
      
      for (t in 1:T) {
        ifelse(t == 1,
               ifelse(abs(phi1) < 1, 
                      y[t] <- rnorm(1, mean = 0, 
                                    sd = sqrt((sigma_e^2)/(1-phi1^2))),
                      y[t] <- c(phi1 * y0 + data$epsilon[t])),
               y[t] <- c(phi1 * y[as.numeric(t - 1)] + data$epsilon[t]))
      }
      
      data <- data.frame(data, y)
      
      # calcuating first lags of dependant and independant variables
      
      data <- data %>% 
        dplyr::mutate(ilag1 = lag(i), tlag1 = lag(t), ylag1 = lag(y))
      
      ## OLS
      ols_est <- ols(data$y, data$t, res = TRUE)
      phi1_ols <- ols(ols_est$residuals_ols$res[-1], 
                      ols_est$residuals_ols$res_lag1[-1], res = FALSE)
      ols_est <- list(beta_ols = ols_est$estimates, 
                      phi1_ols = phi1_ols$estimates[2,], 
                      sigma_e2_ols = ols_est$sigma_e2_ols)
      
      # Saving OLS estimate of phi1 for this replication
      
      phi1_ols_sim[[i]]$phi1[r] <- phi1_ols$estimates$coefficient[2] 
      phi1_ols_sim[[i]]$st.error[r] <- phi1_ols$estimates$st.error[2]
      phi1_ols_sim[[i]]$t_stat[r] <- phi1_ols$estimates$t_stat[2]
      
      
      ## MLE
      mle_est <- mle(loglikelihood = loglikelihood, data = data, 
                     ols_est = ols_est)
      
      # Saving MLE estimate of alpha
      alpha_mle_sim[[i]]$phi1[r] <- mle_est$beta_mle$coefficient[1] 
      alpha_mle_sim[[i]]$st.error[r] <- mle_est$beta_mle$st.error[1]
      alpha_mle_sim[[i]]$t_stat[r] <- mle_est$beta_mle$t_stat[1]
      
      # Saving MLE estimate of beta
      beta_mle_sim[[i]]$phi1[r] <- mle_est$beta_mle$coefficient[2] 
      beta_mle_sim[[i]]$st.error[r] <- mle_est$beta_mle$st.error[2]
      beta_mle_sim[[i]]$t_stat[r] <- mle_est$beta_mle$t_stat[2]
      
      # Saving MLE estimate of phi
      phi1_mle_sim[[i]]$phi1[r] <- mle_est$phi1_mle$coefficient[1] 
      phi1_mle_sim[[i]]$st.error[r] <- mle_est$phi1_mle$st.error[1]
      phi1_mle_sim[[i]]$t_stat[r] <- mle_est$phi1_mle$t_stat[1]
      
    }
  }
  
  return(list(phi1_ols_sim = phi1_ols_sim, 
              alpha_mle_sim = alpha_mle_sim, 
              beta_mle_sim = beta_mle_sim, 
              phi1_mle_sim = phi1_mle_sim))
}

###---- Sample  100: Tables and Graphs ---------------------------
# 3.2.1

T = 100

phi1_ols_3.2 <- data.frame(Phi_1 = rep(0,4), coefficient = rep(0,4), 
                           st.dev = rep(0,4), t_stat = rep(0,4))
phi1_mle_3.2 <- data.frame(Phi_1 = rep(0,4), coefficient = rep(0,4), 
                           st.dev = rep(0,4), t_stat = rep(0,4))
alpha_mle_3.2 <- data.frame(Phi_1 = rep(0,4), coefficient = rep(0,4), 
                           st.dev = rep(0,4), t_stat = rep(0,4))
beta_mle_3.2 <- data.frame(Phi_1 = rep(0,4), coefficient = rep(0,4), 
                           st.dev = rep(0,4), t_stat = rep(0,4))

for (i in 1:4){
  
  phi1 <- phi[i] # defining value of phi1
  
  for (r in 1:2){
    
    ## Replications of data
    data <- data.frame(epsilon = rnorm(T, mean = mean_e, 
                                       sd = sigma_e), 
                       i = rep(1,T), t = 1:T)
    
    y <- NULL
    
    for (t in 1:T) {
      ifelse(t == 1,
             ifelse(abs(phi1) < 1, 
                    y[t] <- rnorm(1, mean = 0, 
                                  sd = sqrt((sigma_e^2)/(1-phi1^2))),
                    y[t] <- c(phi1 * y0 + data$epsilon[t])),
             y[t] <- c(phi1 * y[as.numeric(t - 1)] + data$epsilon[t]))
    }
    
    data <- data.frame(data, y)
    
    # calcuating first lags of dependant and independant variables
    
    data <- data %>% 
      dplyr::mutate(ilag1 = lag(i), tlag1 = lag(t), ylag1 = lag(y))
    
    mypath <- file.path(paste("yt_r_", r, "_phi1_", phi[i],".jpeg", sep = ""))
    jpeg(file = mypath)
    plot(data$t, data$y, xlab = "t", ylab = "yt",
         main = paste("Scatterplot: yt for phi1 = ", phi[i]))
    dev.off()
    
    # 3.2.2
    
    # OLS estimation using or function
    ols_est <- ols(depvar = data$y, indvar = data$t, res = TRUE)
    phi1_ols <- ols(depvar = ols_est$residuals_ols$res[-1], 
                    indvar = ols_est$residuals_ols$res_lag1[-1], 
                    res = FALSE)
    ols_est <- list(beta_ols = ols_est$estimates, 
                    phi1_ols = phi1_ols$estimates[2,], 
                    sigma_e2_ols = ols_est$sigma_e2_ols)
    
    # MLE estimation using our function
    mle_est <- mle(loglikelihood = loglikelihood, data = data, 
                   ols_est = ols_est)
  }
  
  # Saving the results for second replication
  phi1_ols_3.2[i,] <- c(paste("Phi1 = ", phi[i]),
                        round(phi1_ols$estimates$coefficient[2], 4), 
                        round(phi1_ols$estimates$st.error[2], 4), 
                        round(phi1_ols$estimates$t_stat[2], 4))
  phi1_mle_3.2[i,] <- c(paste("Phi1 = ", phi[i]),
                        round(mle_est$phi1_mle$coefficient[1], 4), 
                        round(mle_est$phi1_mle$st.error[1], 4), 
                        round(mle_est$phi1_mle$t_stat[1], 4))
  alpha_mle_3.2[i,] <- c(paste("Phi1 = ", phi[i]),
                         round(mle_est$beta_mle$coefficient[1], 4), 
                         round(mle_est$beta_mle$st.error[1], 4), 
                         round(mle_est$beta_mle$t_stat[1], 4))
  beta_mle_3.2[i,] <- c(paste("Phi1 = ", phi[i]),
                        round(mle_est$beta_mle$coefficient[2], 4), 
                        round(mle_est$beta_mle$st.error[2], 4), 
                        round(mle_est$beta_mle$t_stat[2], 4))
  
}

Table_2 <- list(alpha_mle = alpha_mle_3.2, beta_mle = beta_mle_3.2, 
                phi1_mle = phi1_mle_3.2, phi1_ols = phi1_ols_3.2)

###----3.2.3 Simulation Tables---------------------------------------------------

## Saving results of simulations with all sample size
#sim_results_T_10 <- simulation(R = 500, T = 10)
#sim_results_T_80 <- simulation(R = 500, T = 80)
#sim_results_T_200 <- simulation(R = 500, T = 200)
#sim_results_T_320 <- simulation(R = 500, T = 320)

# save(sim_results_T_10, sim_results_T_80, sim_results_T_200, 
#     sim_results_T_320, file="simulation_results.RData")

# Results have been saved use load("simulation_results.RData") to load
# Please note that the simulation function often gives error for large 
# number of replications (due to singular hessian encoutered in MLE).
# The results have been saved to avoid running into errors each time. 

load("simulation_results.RData")
sim_results_all <- list(sim_results_T_10, sim_results_T_80, 
                        sim_results_T_200, sim_results_T_320)

# Tables

sample_size <- c(10, 80, 200, 320)

# SUMMARY STATISTICS PHI 1 OLS ESTIMATE
summary_phi1_ols <- data.frame(Phi_1 = rep(0, 16), 
                               Sample_size = rep(0, 16),
                               Mean = rep(0, 16),
                               St.dev = rep(0, 16), 
                               t_stat = rep(1, 16),
                               Mean_Bias = rep(0, 16),
                               RMSE = rep(0, 16),
                               Quantile_5 = rep(0, 16), 
                               Median = rep(0, 16), 
                               Quantile_95 = rep(0, 16),
                               asymptotic_var = rep(0, 16))

j <- NULL

for (t in 1:4){
  ifelse(t == 1, j <- 1, ifelse(t == 2, j<-5, 
                                ifelse(t == 3, j <- 9, j <- 13)))
  
  for (i in 1:4){
    
    sim_res <- sim_results_all[[t]]$phi1_ols_sim[[i]]$phi1
    
    summary_phi1_ols[j,] <- c(paste("Phi1 = ", phi[i]), sample_size[t],
                              round(mean(sim_res), 4), 
                              round(sd(sim_res), 4),
                              round((mean(sim_res) - phi[i])/ 
                                sd(sim_res), 4),
                              round(mean(sim_res - phi[i]), 4),
                              round(sqrt(mean((sim_res) ^ 2)), 4),
                              round(quantile(sim_res, 0.05), 4),
                              round(median(sim_res), 4),
                              round(quantile(sim_res, 0.95), 4),
                              round(var(sqrt(sample_size[t]) * 
                                      (sim_res - phi[i])), 4))
    j=j+1
  }
}

# SUMMARY STATISTICS PHI 1 MLE ESTIMATE

summary_phi1_mle <- data.frame(Phi_1 = rep(0,16), 
                               Sample_size = rep(0,16),
                               Mean = rep(0,16), 
                               St.dev = rep(0,16), 
                               t_stat = rep(1, 16),
                               Mean_Bias = rep(0,16),
                               RMSE = rep(0,16),
                               Quantile_5 = rep(0,16), 
                               Median = rep(0,16), 
                               Quantile_95 = rep(0,16),
                               asymptotic_var = rep(0,16))

j <- NULL

for (t in 1:4){
  ifelse(t == 1, j <- 1, ifelse(t == 2, j<-5, 
                                ifelse(t == 3, j <- 9, j <- 13)))
  
  for (i in 1:4){
    
    sim_res <- sim_results_all[[t]]$phi1_mle_sim[[i]]$phi1
    
    summary_phi1_mle[j,] <- c(paste("Phi1 = ", phi[i]), sample_size[t],
                              round(mean(sim_res), 4), 
                              round(sd(sim_res), 4),
                              round((mean(sim_res) - phi[i])/
                                      sd(sim_res), 4),
                              round(mean(sim_res - phi[i]), 4),
                              round(sqrt(mean((sim_res) ^ 2)), 4),
                              round(quantile(sim_res, 0.05), 4),
                              round(median(sim_res), 4),
                              round(quantile(sim_res, 0.95), 4),
                              round(var(sqrt(sample_size[t]) * 
                                          (sim_res - phi[i])), 4))
    j=j+1
  }
}


# SUMMARY STATISTICS ALPHA MLE ESTIMATE

summary_alpha_mle <- data.frame(Phi_1 = rep(0,16), 
                               Sample_size = rep(0,16),
                               Mean = rep(0,16), 
                               St.dev = rep(0,16), 
                               t_stat = rep(1, 16),
                               Mean_Bias = rep(0,16),
                               RMSE = rep(0,16),
                               Quantile_5 = rep(0,16), 
                               Median = rep(0,16), 
                               Quantile_95 = rep(0,16),
                               asymptotic_var = rep(0,16))

j <- NULL

for (t in 1:4){
  ifelse(t == 1, j <- 1, ifelse(t == 2, j<-5, 
                                ifelse(t == 3, j <- 9, j <- 13)))
  
  for (i in 1:4){
    
    sim_res <- sim_results_all[[t]]$alpha_mle_sim[[i]]$phi1
    
    summary_alpha_mle[j,] <- c(paste("Phi1 = ", phi[i]), sample_size[t],
                              round(mean(sim_res), 4), 
                              round(sd(sim_res), 4),
                              round(mean(sim_res) / sd(sim_res), 4),
                              round(mean(sim_res - phi[i]), 4),
                              round(sqrt(mean((sim_res) ^ 2)), 4),
                              round(quantile(sim_res, 0.05), 4),
                              round(median(sim_res), 4),
                              round(quantile(sim_res, 0.95), 4),
                              round(var(sqrt(sample_size[t]) * 
                                          (sim_res - phi[i])), 4))
    j=j+1
  }
}


# SUMMARY STATISTICS BETA MLE ESTIMATE

summary_beta_mle <- data.frame(Phi_1 = rep(0,16), 
                               Sample_size = rep(0,16),
                               Mean = rep(0,16), 
                               St.dev = rep(0,16), 
                               t_stat = rep(1, 16),
                               Mean_Bias = rep(0,16),
                               RMSE = rep(0,16),
                               Quantile_5 = rep(0,16), 
                               Median = rep(0,16), 
                               Quantile_95 = rep(0,16),
                               asymptotic_var = rep(0,16))

j <- NULL

for (t in 1:4){
  ifelse(t == 1, j <- 1, ifelse(t == 2, j<-5, 
                                ifelse(t == 3, j <- 9, j <- 13)))
  
  for (i in 1:4){
    
    sim_res <- sim_results_all[[t]]$beta_mle_sim[[i]]$phi1
    
    summary_beta_mle[j,] <- c(paste("Phi1 = ", phi[i]), sample_size[t],
                              round(mean(sim_res), 4), 
                              round(sd(sim_res), 4),
                              round(mean(sim_res) / sd(sim_res), 4),
                              round(mean(sim_res - phi[i]), 4),
                              round(sqrt(mean((sim_res) ^ 2)), 4),
                              round(quantile(sim_res, 0.05), 4),
                              round(median(sim_res), 4),
                              round(quantile(sim_res, 0.95), 4),
                              round(var(sqrt(sample_size[t]) * 
                                          (sim_res - phi[i])), 4))
    j=j+1
  }
}

Table_3 <- list(alpha_mle = summary_alpha_mle, 
                beta_mle = summary_beta_mle, 
                phi1_mle = summary_phi1_mle, 
                phi1_ols = summary_phi1_ols)

###---- 3.2.3 Graphs ---------------------------

## Generating and saving the plots for T = 200

# alpha ML

y1 <- sim_results_T_200$alpha_mle_sim$phi_0.5$phi1
y2 <- sim_results_T_200$alpha_mle_sim$phi_0.9$phi1
y3 <- sim_results_T_200$alpha_mle_sim$phi_0.99$phi1
y4 <- sim_results_T_200$alpha_mle_sim$phi_1$phi1

xmax = 40
ymax = 0.47

jpeg('alpha_mle.jpg')

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Sample Size = 200",
     col='red', 
     main="Maximum Likelihood Estimate of alpha")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', col='violetred', main = "")

legend("topright", legend = c(paste("Phi_1 = ", phi[1:4])), 
       col = c("red", "forestgreen", "darkblue", "violetred"), 
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1, 
       text.width = 20)

dev.off()

# beta ML

y1 <- sim_results_T_200$beta_mle_sim$phi_0.5$phi1
y2 <- sim_results_T_200$beta_mle_sim$phi_0.9$phi1
y3 <- sim_results_T_200$beta_mle_sim$phi_0.99$phi1
y4 <- sim_results_T_200$beta_mle_sim$phi_1$phi1

xmax = 0.7
ymax = 20

jpeg('beta_mle.jpg')

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Sample Size = 200",
     col='red', 
     main="Maximum Likelihood Estimate of beta")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', col='violetred', main = "")

legend("topright", legend = c(paste("Phi_1 = ", phi[1:4])), 
       col = c("red", "forestgreen", "darkblue", "violetred"), 
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1, 
       text.width = 0.3)

dev.off()

# phi1 ML

y1 <- sim_results_T_200$phi1_mle_sim$phi_0.5$phi1
y2 <- sim_results_T_200$phi1_mle_sim$phi_0.9$phi1
y3 <- sim_results_T_200$phi1_mle_sim$phi_0.99$phi1
y4 <- sim_results_T_200$phi1_mle_sim$phi_1$phi1

xmax = 1
ymax = 15

jpeg('phi1_mle.jpg')

plot(density(y1), ylim = c(0,ymax), xlim = c(0, xmax),
     type='l', xlab = "No. of Replications = 500, Sample Size = 200",
     col='red', 
     main="Maximum Likelihood Estimate of Phi 1")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(0, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(0, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(0, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', col='violetred', main = "")

legend("topleft", legend = c(paste("Phi_1 = ", phi[1:4])), 
       col = c("red", "forestgreen", "darkblue", "violetred"), 
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1, 
       text.width = 0.2)

dev.off()

# phi1 OLS

y1 <- sim_results_T_200$phi1_ols_sim$phi_0.5$phi1
y2 <- sim_results_T_200$phi1_ols_sim$phi_0.9$phi1
y3 <- sim_results_T_200$phi1_ols_sim$phi_0.99$phi1
y4 <- sim_results_T_200$phi1_ols_sim$phi_1$phi1

xmax = 1
ymax = 15

jpeg('phi1_ols.jpg')

plot(density(y1), ylim = c(0,ymax), xlim = c(0, xmax),
     type='l', xlab = "No. of Replications = 500, Sample Size = 200",
     col='red', 
     main="Ordinary Least Square Estimate of Phi 1")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(0, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(0, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', 
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(0, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', col='violetred', main = "")

legend("topleft", legend = c(paste("Phi_1 = ", phi[1:4])), 
       col = c("red", "forestgreen", "darkblue", "violetred"), 
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1, 
       text.width = 0.2)

dev.off()

###---- Graphs for all sample sizes ---------------------------

#Phi1 = 0.5
# alpha ML

y1 <- sqrt(10) * (sim_results_T_10$alpha_mle_sim$phi_0.5$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$alpha_mle_sim$phi_0.5$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$alpha_mle_sim$phi_0.5$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$alpha_mle_sim$phi_0.5$phi1 - 0)

xmax = 40
ymax = 0.05

png('alpha_mle_finite_sample_phi_0.5.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.5",
     col='red', lty = 1,
     main="ML estimate of alpha for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 18)

dev.off()

# Beta ML

y1 <- sqrt(10) * (sim_results_T_10$beta_mle_sim$phi_0.5$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$beta_mle_sim$phi_0.5$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$beta_mle_sim$phi_0.5$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$beta_mle_sim$phi_0.5$phi1 - 0)

xmax = 5
ymax = 6.2

png('beta_mle_finite_sample_phi_0.5.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.5",
     col='red', lty = 1,
     main="ML estimate of beta for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 3)

dev.off()

# phi1 ML

y1 <- sqrt(10) * (sim_results_T_10$phi1_mle_sim$phi_0.5$phi1 - 0.5)
y2 <- sqrt(80) * (sim_results_T_80$phi1_mle_sim$phi_0.5$phi1 - 0.5)
y3 <- sqrt(200) * (sim_results_T_200$phi1_mle_sim$phi_0.5$phi1 - 0.5)
y4 <- sqrt(320) * (sim_results_T_320$phi1_mle_sim$phi_0.5$phi1 - 0.5)

xmin = -5
xmax = 2
ymax = 0.55

png('phi1_mle_finite_sample_phi_0.5.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.5",
     col='red', lty = 1,
     main="ML estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

# phi1 OLS

y1 <- sqrt(10) * (sim_results_T_10$phi1_ols_sim$phi_0.5$phi1 - 0.5)
y2 <- sqrt(80) * (sim_results_T_80$phi1_ols_sim$phi_0.5$phi1 - 0.5)
y3 <- sqrt(200) * (sim_results_T_200$phi1_ols_sim$phi_0.5$phi1 - 0.5)
y4 <- sqrt(320) * (sim_results_T_320$phi1_ols_sim$phi_0.5$phi1 - 0.5)

xmin = -5
xmax = 2
ymax = 0.52

png('phi1_ols_finite_sample_phi_0.5.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.5",
     col='red', lty = 1,
     main="OLS estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

#Phi1 = 0.9
# alpha ML

y1 <- sqrt(10) * (sim_results_T_10$alpha_mle_sim$phi_0.9$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$alpha_mle_sim$phi_0.9$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$alpha_mle_sim$phi_0.9$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$alpha_mle_sim$phi_0.9$phi1 - 0)

xmax = 160
ymax = 0.015

png('alpha_mle_finite_sample_phi_0.9.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.9",
     col='red', lty = 1,
     main="ML estimate of alpha for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 70)

dev.off()

# Beta ML

y1 <- sqrt(10) * (sim_results_T_10$beta_mle_sim$phi_0.9$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$beta_mle_sim$phi_0.9$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$beta_mle_sim$phi_0.9$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$beta_mle_sim$phi_0.9$phi1 - 0)

xmax = 7
ymax = 1.5

png('beta_mle_finite_sample_phi_0.9.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.9",
     col='red', lty = 1,
     main="ML estimate of beta for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 4)

dev.off()

# phi1 ML

y1 <- sqrt(10) * (sim_results_T_10$phi1_mle_sim$phi_0.9$phi1 - 0.9)
y2 <- sqrt(80) * (sim_results_T_80$phi1_mle_sim$phi_0.9$phi1 - 0.9)
y3 <- sqrt(200) * (sim_results_T_200$phi1_mle_sim$phi_0.9$phi1 - 0.9)
y4 <- sqrt(320) * (sim_results_T_320$phi1_mle_sim$phi_0.9$phi1 - 0.9)

xmin = -6.5
xmax = 1
ymax = 0.8

png('phi1_mle_finite_sample_phi_0.9.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.9",
     col='red', lty = 1,
     main="ML estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

# phi1 OLS

y1 <- sqrt(10) * (sim_results_T_10$phi1_ols_sim$phi_0.9$phi1 - 0.9)
y2 <- sqrt(80) * (sim_results_T_80$phi1_ols_sim$phi_0.9$phi1 - 0.9)
y3 <- sqrt(200) * (sim_results_T_200$phi1_ols_sim$phi_0.9$phi1 - 0.9)
y4 <- sqrt(320) * (sim_results_T_320$phi1_ols_sim$phi_0.9$phi1 - 0.9)

xmin = -7
xmax = 1.2
ymax = 0.9

png('phi1_ols_finite_sample_phi_0.9.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.9",
     col='red', lty = 1,
     main="OLS estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

#Phi1 = 0.99
# alpha ML

y1 <- sqrt(10) * (sim_results_T_10$alpha_mle_sim$phi_0.99$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$alpha_mle_sim$phi_0.99$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$alpha_mle_sim$phi_0.99$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$alpha_mle_sim$phi_0.99$phi1 - 0)

xmax = 1000
ymax = 0.007

png('alpha_mle_finite_sample_phi_0.99.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.99",
     col='red', lty = 1,
     main="ML estimate of alpha for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 450)

dev.off()

# Beta ML

y1 <- sqrt(10) * (sim_results_T_10$beta_mle_sim$phi_0.99$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$beta_mle_sim$phi_0.99$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$beta_mle_sim$phi_0.99$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$beta_mle_sim$phi_0.99$phi1 - 0)

xmax = 9
ymax = 0.25

png('beta_mle_finite_sample_phi_0.99.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.99",
     col='red', lty = 1,
     main="ML estimate of beta for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 4)

dev.off()

# phi1 ML

y1 <- sqrt(10) * (sim_results_T_10$phi1_mle_sim$phi_0.99$phi1 - 0.99)
y2 <- sqrt(80) * (sim_results_T_80$phi1_mle_sim$phi_0.99$phi1 - 0.99)
y3 <- sqrt(200) * (sim_results_T_200$phi1_mle_sim$phi_0.99$phi1 - 0.99)
y4 <- sqrt(320) * (sim_results_T_320$phi1_mle_sim$phi_0.99$phi1 - 0.99)

xmin = -6.5
xmax = 0.5
ymax = 1

png('phi1_mle_finite_sample_phi_0.99.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.99",
     col='red', lty = 1,
     main="ML estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

# phi1 OLS

y1 <- sqrt(10) * (sim_results_T_10$phi1_ols_sim$phi_0.99$phi1 - 0.99)
y2 <- sqrt(80) * (sim_results_T_80$phi1_ols_sim$phi_0.99$phi1 - 0.99)
y3 <- sqrt(200) * (sim_results_T_200$phi1_ols_sim$phi_0.99$phi1 - 0.99)
y4 <- sqrt(320) * (sim_results_T_320$phi1_ols_sim$phi_0.99$phi1 - 0.99)

xmin = -6.5
xmax = 0.5
ymax = 1.5

png('phi1_ols_finite_sample_phi_0.99.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 0.99",
     col='red', lty = 1,
     main="OLS estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

#Phi1 = 1
# alpha ML

y1 <- sqrt(10) * (sim_results_T_10$alpha_mle_sim$phi_1$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$alpha_mle_sim$phi_1$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$alpha_mle_sim$phi_1$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$alpha_mle_sim$phi_1$phi1 - 0)

xmax = 900
ymax = 0.035

png('alpha_mle_finite_sample_phi_1.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 1",
     col='red', lty = 1,
     main="ML estimate of alpha for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 450)

dev.off()

# Beta ML

y1 <- sqrt(10) * (sim_results_T_10$beta_mle_sim$phi_1$phi1 - 0)
y2 <- sqrt(80) * (sim_results_T_80$beta_mle_sim$phi_1$phi1 - 0)
y3 <- sqrt(200) * (sim_results_T_200$beta_mle_sim$phi_1$phi1 - 0)
y4 <- sqrt(320) * (sim_results_T_320$beta_mle_sim$phi_1$phi1 - 0)

xmax = 10
ymax = 0.15

png('beta_mle_finite_sample_phi_1.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 1",
     col='red', lty = 1,
     main="ML estimate of beta for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(-xmax, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topright", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 5)

dev.off()

# phi1 ML

y1 <- sqrt(10) * (sim_results_T_10$phi1_mle_sim$phi_1$phi1 - 1)
y2 <- sqrt(80) * (sim_results_T_80$phi1_mle_sim$phi_1$phi1 - 1)
y3 <- sqrt(200) * (sim_results_T_200$phi1_mle_sim$phi_1$phi1 - 1)
y4 <- sqrt(320) * (sim_results_T_320$phi1_mle_sim$phi_1$phi1 - 1)

xmin = -6.5
xmax = 0.5
ymax = 1

png('phi1_mle_finite_sample_phi_1.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 1",
     col='red', lty = 1,
     main="ML estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

# phi1 OLS

y1 <- sqrt(10) * (sim_results_T_10$phi1_ols_sim$phi_1$phi1 - 1)
y2 <- sqrt(80) * (sim_results_T_80$phi1_ols_sim$phi_1$phi1 - 1)
y3 <- sqrt(200) * (sim_results_T_200$phi1_ols_sim$phi_1$phi1 - 1)
y4 <- sqrt(320) * (sim_results_T_320$phi1_ols_sim$phi_1$phi1 - 1)

xmin = -6.5
xmax = 0.5
ymax = 1.5

png('phi1_ols_finite_sample_phi_1.jpg', width = 600, height = 600)

plot(density(y1), ylim = c(0,ymax), xlim = c(xmin, xmax),
     type='l', xlab = "No. of Replications = 500, Phi 1 = 1",
     col='red', lty = 1,
     main="OLS estimate of phi 1 for different sample sizes")
par(new = TRUE)
plot(density(y2), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1, 
     col='forestgreen', main = "")
par(new = TRUE)
plot(density(y3), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "",
     type='l', lty = 1,
     col='darkblue', main = "")
par(new = TRUE)
plot(density(y4), ylim = c(0,ymax), xlim = c(xmin, xmax),
     axes = FALSE, xlab = "", ylab = "", lty = 1, 
     type='l', col='sienna4', main = "")

legend("topleft", legend = c(paste("Sample size = ", sample_size[1:4])),
       col = c("red", "forestgreen", "darkblue", "sienna4"),
       cex = 0.8, text.font = 4, lty = 1, y.intersp = 1,
       text.width = 2)

dev.off()

###---- QQ plots ---------------------------

# for phi1 = 0.5

jpeg('qqplot_phi1_ols_for_phi0.5.jpg')
qqnorm(sim_results_T_200$phi1_ols_sim$phi_0.5$phi1, 
       main = "OLS estimate of Phi 1 for phi1 = 0.5")
qqline(sim_results_T_200$phi1_ols_sim$phi_0.5$phi1)
dev.off()

jpeg('qqplot_phi1_mle_for_phi0.5.jpg')
qqnorm(sim_results_T_200$phi1_mle_sim$phi_0.5$phi1,
       main = "ML estimate of Phi 1 for phi1 = 0.5")
qqline(sim_results_T_200$phi1_mle_sim$phi_0.5$phi1)
dev.off()

jpeg('qqplot_alpha_mle_for_phi0.5.jpg')
qqnorm(sim_results_T_200$alpha_mle_sim$phi_0.5$phi1, 
       main = "ML estimate of alpha for phi1 = 0.5")
qqline(sim_results_T_200$alpha_mle_sim$phi_0.5$phi1)
dev.off()

jpeg('qqplot_beta_mle_for_phi0.5.jpg')
qqnorm(sim_results_T_200$beta_mle_sim$phi_0.5$phi1, 
       main = "ML estimate of beta for phi1 = 0.5")
qqline(sim_results_T_200$beta_mle_sim$phi_0.5$phi1)
dev.off()

# for phi1 = 0.9

jpeg('qqplot_phi1_ols_for_phi0.9.jpg')
qqnorm(sim_results_T_200$phi1_ols_sim$phi_0.9$phi1, 
       main = "OLS estimate of Phi 1 for phi1 = 0.9")
qqline(sim_results_T_200$phi1_ols_sim$phi_0.9$phi1)
dev.off()

jpeg('qqplot_phi1_mle_for_phi0.9.jpg')
qqnorm(sim_results_T_200$phi1_mle_sim$phi_0.9$phi1, 
       main = "ML estimate of phi 1 for phi1 = 0.9")
qqline(sim_results_T_200$phi1_mle_sim$phi_0.9$phi1)
dev.off()

jpeg('qqplot_alpha_mle_for_phi0.9.jpg')
qqnorm(sim_results_T_200$alpha_mle_sim$phi_0.9$phi1, 
       main = "ML estimate of alpha for phi1 = 0.9")
qqline(sim_results_T_200$alpha_mle_sim$phi_0.9$phi1)
dev.off()

jpeg('qqplot_beta_mle_for_phi0.9.jpg')
qqnorm(sim_results_T_200$beta_mle_sim$phi_0.9$phi1, 
       main = "ML estimate of beta for phi1 = 0.9")
qqline(sim_results_T_200$beta_mle_sim$phi_0.9$phi1)
dev.off()

# for phi1 = 0.99

jpeg('qqplot_phi1_ols_for_phi0.99.jpg')
qqnorm(sim_results_T_200$phi1_ols_sim$phi_0.99$phi1, 
       main = "OLS estimate of Phi 1 for phi1 = 0.99")
qqline(sim_results_T_200$phi1_ols_sim$phi_0.99$phi1)
dev.off()

jpeg('qqplot_phi1_mle_for_phi0.99.jpg')
qqnorm(sim_results_T_200$phi1_mle_sim$phi_0.99$phi1, 
       main = "ML estimate of phi 1 for phi1 = 0.99")
qqline(sim_results_T_200$phi1_mle_sim$phi_0.99$phi1)
dev.off()

jpeg('qqplot_alpha_mle_for_phi0.99.jpg')
qqnorm(sim_results_T_200$alpha_mle_sim$phi_0.99$phi1, 
       main = "ML estimate of alpha for phi1 = 0.99")
qqline(sim_results_T_200$alpha_mle_sim$phi_0.99$phi1)
dev.off()

jpeg('qqplot_beta_mle_for_phi0.99.jpg')
qqnorm(sim_results_T_200$beta_mle_sim$phi_0.99$phi1, 
       main = "ML estimate of beta for phi1 = 0.99")
qqline(sim_results_T_200$beta_mle_sim$phi_0.99$phi1)
dev.off()

# for phi1 = 1

jpeg('qqplot_phi1_ols_for_phi1.jpg')
qqnorm(sim_results_T_200$phi1_ols_sim$phi_1$phi1, 
       main = "OLS estimate of Phi 1 for phi1 = 1")
qqline(sim_results_T_200$phi1_ols_sim$phi_1$phi1)
dev.off()

jpeg('qqplot_phi1_mle_for_phi1.jpg')
qqnorm(sim_results_T_200$phi1_mle_sim$phi_1$phi1, 
       main = "ML estimate of Phi 1 for phi1 = 1")
qqline(sim_results_T_200$phi1_mle_sim$phi_1$phi1)
dev.off()

jpeg('qqplot_alpha_mle_for_phi1.jpg')
qqnorm(sim_results_T_200$alpha_mle_sim$phi_1$phi1, 
       main = "ML estimate of alpha for phi1 = 1")
qqline(sim_results_T_200$alpha_mle_sim$phi_1$phi1)
dev.off()

jpeg('qqplot_beta_mle_for_phi1.jpg')
qqnorm(sim_results_T_200$beta_mle_sim$phi_1$phi1, 
       main = "ML estimate of beta for phi1 = 1")
qqline(sim_results_T_200$beta_mle_sim$phi_1$phi1)
dev.off()

# The end