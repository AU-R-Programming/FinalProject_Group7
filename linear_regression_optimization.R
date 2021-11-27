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

initial_betas <- rep(NA, n_predictors) # vector to create the initial solution

initial_betas[1] <- mean(y) # fill the first element with the mean of the response

# This for loop fill the initial solution
for (i in 2:n_predictors) {
  initial_betas[i] <- cov(y, x[,i])/ var(x[,i])
}

# Computing the error or loss function
error <- t(y - (x %*% initial_betas))%*%(y - x %*% initial_betas)

# Start gradient descent 

learning_rate <- 0.1

for (it in 1:50) { # This for loop repeats the algorithm
betas_update <- rep(NA, n_predictors) # initialize the betas for update

for(i in 1:n_predictors){ # This for loop gets the gradient for each beta
  gradient_beta <- -2 * sum(x[,i] * (y - sum(initial_betas * x)))/length(x)/n_obs # This formula comes from gradient descent 
  betas_update[i] <- initial_betas[i] - learning_rate * gradient_beta # Update the ith beta
}

error_gradient <-  t(y - (x %*% betas_update))%*%(y - x %*% betas_update) # Compute the loss function

if (error_gradient < error) { # Evaluate if there is improvement
  initial_betas <- betas_update
  error <- error_gradient
}
print(error) 
}







