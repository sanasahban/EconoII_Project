
#Test Data
ytest <- data.frame(y = c(1,2,4,6,4,12,11,13,14,16,17))
ytest.mat <- as.matrix(ytest)

#OLS
X <- data.frame(i = matrix(1,10,1),t=data$t)
x.mat <- as.matrix(X)
beta_ols <- (solve(t(x.mat)%*%x.mat))%*%(t(x.mat)%*%y)
rownames(beta_ols) <- c("intercept", "t")
