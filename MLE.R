#---------------------------------------------------------------------

#Data Generating Process (DGP)

aplha = 0
beta = 0
phi2 = 0
mean_e = 0
sigma_e =  3

y0 = 0 #only if phi1=1
phi <- c(0.5, 0.9, 0.99, 1) #use phi1=phi for each cases

T = 10

data <- data.frame(epsilon = rnorm(T, mean = mean_e, sd = sigma_e), 
                   i = rep(1,T), t = 1:T)

y <- NA 

for (t in 1:T) {
  ifelse (t == 1,
          ifelse(abs(phi[1]) < 1, 
                 y[t] <- rnorm(1, mean = 0, sd = 
                                 sqrt((sigma_e^2)/(1-phi[1]^2))),
                 y[t] <- c(phi[1] * y0 + data$epsilon[t])
                 ),
          y[t] <- c(phi[1] * y[as.numeric(t - 1)] + data$epsilon[t])
  )
}

data <- data.frame(data, y)
library(dplyr)
data <- data %>% mutate(ilag1 = lag(i), tlag1 = lag(t), ylag1 = lag(y))

#---------------------------------------------------------------------

#MLE - maximum likelihood estimation

#create matrix of regressors where i is column of 1s and t is trend variable
x <- cbind(data$i, data$t)
xlag <- cbind(data$ilag1[-1] ,data$tlag1[-1])

#loglikelihood FOCS
phi1_mle = 0.5 #test
beta_mle <- matrix( c(0,0), ncol=1)  #test
sigma_e_mle = 3

#FOC wrt beta

(1-phi1_mle^2)*x[1,]%*%(y[1,]-x[1,]%*%beta_mle) + 
  sum(t((x[-1,]-phi1_mle*xlag[-1,]))%*%
        (y[-1]-x[-1,]%*%beta_mle-phi1_mle*(data$ylag1[-1]-xlag[-1,]%*%beta_mle)))==0

#Log-likelihood function

l <- function(data, x, xlag, theta){ #theta is a vector of parameters as specied above
  -T/2*log(2*pi) - T*log(theta[4]) + 1/2*log(1 - theta[1]^2) -
    (1/(2*theta[4]))*((1 - theta[1]^2)*
                           (data$y[1] - x[1,]%*%theta[2:3])^2) + 
                              sum((data$y[-1] - x[-1,]%*%theta[2:3] - 
                                     theta[1]*(data$ylag1[-1] - 
                                        xlag[-1,]%*%theta[2:3]))^2)
}

#initial value
theta = as.matrix(data.frame(phi1_mle = 0.5, alpha_mle = 0, 
                             beta_mle = 0, sigma_e_mle = 3))
#Optimisation
mle <- optim(theta = c(0.5, 0, 0, 3), l, data = data, method = "BFGS")
