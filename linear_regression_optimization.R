# Import data 
abalone <- read.csv(file = 'abalone.csv') 

# Transform dummy variables
library(fastDummies)
abalone <- dummy_cols(abalone, select_columns = "Sex")

# Create response variable and set of predictors 
y <- abalone$Rings
x <- abalone[, 2:11]
x <- subset(x, select = -Rings)

## Here starts the package
n_obs <- length(y)
x <- cbind('1' = 1, x) # create a column with ones
x <- data.matrix(x) # transform the previous column to a matrix
n_predictors <- ncol(x) # obtained the number of predictors plus the constant

betas <- rep(NA, n_predictors) # vector to create the initial solution

betas[1] <- mean(y) # fill the first element with the mean of the response

# This for loop fill the initial solution
for (i in 2:n_predictors) {
  betas[i] <- cov(y, x[,i])/ var(x[,i])
}

# Computing the error or loss function
error <- t(y - (x %*% betas))%*%(y - (x %*% betas))

# Start gradient descent 

learning_rate <- 0.3

for (it in 1:30000) { # This for loop repeats the algorithm
betas_update <- rep(NA, n_predictors) # initialize the betas for update

for(i in 1:n_predictors){ # This for loop gets the gradient for each beta
  gradient_beta <- (-2/n_obs) * sum(x[, i] * data.matrix(y - x %*% betas)) # This formula comes from gradient descent 
  betas_update[i] <- betas[i] - learning_rate * gradient_beta # Update the ith beta
}

error_gradient <-  t(y - (x %*% betas_update))%*%(y - x %*% betas_update) # Compute the loss function

if (error_gradient < error) { # Evaluate if there is improvement
  betas <- betas_update
  error <- error_gradient
}
} # Finish gradient descent

y_hat <- x %*% betas
y_mean <- mean(y)

SSE <- sum(data.matrix(y - x %*% betas)^2)
SST <- sum((y-y_mean)^2)

sigma_hat_sq = as.numeric(t(y-x%*%betas)%*%(y-x%*%betas)/(n_obs-n_predictors))

## R^2
Rsq = 1-(SSE/SST)

# C_p
Cp = SSE + 2*n_predictors*sigma_hat_sq

# confidence intervals
alpha <- 0.05
x_t_x = solve(t(x)%*%x)
D = rep(NA,n_predictors) # Initializing the diagonal 
for (i in 1 : n_predictors){
  D[i] = x_t_x[i,i] # Filling the diagonal up
}
ss = qnorm((1-alpha/2),mean=0,sd=1)*sqrt(sigma_hat_sq*D) # the right side of the interval
conf_int = cbind(betas-ss,betas+ss) # confidence intervals

DFM = n_predictors-1
DFE = n_obs - n_predictors
SSM = 0
for (i in 1 : n_obs){
  SSM = SSM + (y_hat[i] - y_mean)^2
}
MSM = SSM/DFM
MSE = SSE/DFE

# F_star
F_star = MSM/MSE

# P-value
p_value = 1-pf(F_star, DFM, DFE)

# Plots


