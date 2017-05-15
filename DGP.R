
#Data Generating Process (DGP)
aplha = 0
beta = 0
phi2 = 0
mean_e = 0
sigma_e =  3

y0 = 0
phi <- c(0.5, 0.9, 0.99, 1)

T = 10

data <- data.frame(t = 1:T, epsilon = rnorm(T, mean = mean_e, sd = sigma_e))

y <- matrix(0, T, 1)

for (t in 1:T) {
  ifelse (t == 1,
          y[t] <- c(phi[1] * y0 + data$epsilon[t]),
          y[t] <- c(phi[1] * y[as.numeric(t - 1)] + data$epsilon[t])
  )
}

data <- data.frame(data, y)
library(dplyr)
data <- data %>% mutate(ylag1=lag(y,1),ylag2=lag(y,2))
as.matrix(data)


