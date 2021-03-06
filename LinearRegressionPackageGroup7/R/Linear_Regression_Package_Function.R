#' @title Linear regression using gradient descent method
#'
#' @description Perform a linear regression using the gradient descent
#' method to find the coefficients of the regression. The user can
#' indicates the learning rate and the number of iterations.
#' @param response A \code{vector} of dimension 1 that represent the response
#' variable. It is also known as the "y" variable.
#' @param predictor A \code{dataframe} that contains all the predictors to
#' perform the regression. It is also known as the X variable.
#' @param alpha A \code{numeric} (continuous) used to denote the level of
#' significance when performing the confidence intervals. The default value is
#' 0.05. The alpha parameter must be between 0 and 1.
#' @param learning_rate A \code{numeric} (continuous) used to control the
#' algorithm rate of learning. The default value is 0.1.
#' @param iterations A \code{numeric} (integer) that represents the number of
#' iterations that the gradient descent will be running. As this value
#' increases the error will reduce.
#'
#' @return A \code{list} containing the following attributes:
#' \describe{
#'      \item{Betas}{Regression coefficients}
#'      \item{Confidence interval}{Confidence interval for all the predictors}
#'      \item{Cp}{Cp-mallows}
#'      \item{F statistic}{Statistic F to compute the p-value}
#'      \item{P-value}{p-value}
#'      \item{Resid_plot}{Plot of residuals against fitted values}
#'      \item{qq-plot}{qq-plot of residulas}
#'      \item{hist_resid}{Histogram of residuals}}
#' @export
#' @importFrom stats qnorm
#' @author  Group 7
#' @examples lin_regr(response = x, predictor = y, alpha = 0.05, learning_rate = 0.1,
#' iterations = 1000)
lin_regr = function(response, predictor, alpha = 0.05, learning_rate = 0.1, iterations = 1000){
  if (alpha >= 1 | alpha <= 0) {
    print("Alpha must be between 0 and 1")
  } else if (learning_rate >= 1 | learning_rate <= 0) {
    print("The learning rate must be between 0 and 1")
  } else if (iterations%%1 != 0) {
    print("The number of iterations must be a positive integer")
  } else if (iterations <= 0){
    print("The number of iterations must be a positive integer")
  } else {
    predictor = as.data.frame(predictor)
    col_predictor = colnames(predictor)
    response = as.vector(response)
    n_obs <- length(response)
    predictor <- cbind('1' = 1, predictor) # create a column with ones
    predictor <- data.matrix(predictor) # transform the previous column to a matrix
    n_predictors <- ncol(predictor) # obtained the number of predictors plus the constant

    betas <- rep(NA, n_predictors) # vector to create the initial solution

    betas[1] <- mean(response) # fill the first element with the mean of the response

    # This for loop fill the initial solution
    for (i in 2:n_predictors) {
      betas[i] <- cov(response, predictor[,i])/ var(predictor[,i])
    }

    # Computing the error or loss function
    error <- t(response - (predictor %*% betas))%*%(response - (predictor %*% betas))

    # Start gradient descent

    for (it in 1:iterations) { # This for loop repeats the algorithm
      betas_update <- rep(NA, n_predictors) # initialize the betas for update

      for(i in 1:n_predictors){ # This for loop gets the gradient for each beta
        gradient_beta <- (-2/n_obs) * sum(predictor[, i] * data.matrix(response - predictor %*% betas)) # This formula comes from gradient descent
        betas_update[i] <- betas[i] - learning_rate * gradient_beta # Update the ith beta
      }

      error_gradient <-  t(response - (predictor %*% betas_update))%*%(response - predictor %*% betas_update) # Compute the loss function

      if (error_gradient < error) { # Evaluate if there is improvement
        betas <- betas_update
        error <- error_gradient
      }
    } # Finish gradient descent

    y_hat <- predictor %*% betas
    y_mean <- mean(response)

    SSE <- sum(data.matrix(response - predictor %*% betas)^2)
    SST <- sum((response-y_mean)^2)

    sigma_hat_sq = as.numeric(t(response-predictor%*%betas)%*%(response-predictor%*%betas)/(n_obs-n_predictors))

    ## R^2
    Rsq = 1-(SSE/SST)

    # C_p
    Cp = SSE + 2*n_predictors*sigma_hat_sq

    # confidence intervals
    x_t_x = solve(t(predictor)%*%predictor)
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

    residual = response - y_hat

    # Plots
    par(mfrow=c(1,1))
    plot1 = plot(y_hat,residual,xlab = "Fitted values",ylab = "Residual",main = "Residuals vs Fitted-values")
    plot2 = abline(h=0,col="red")
    plot3 = qqnorm(residual,main = "QQ Plot of Residual")
    plot4 = qqline(residual, col = "red")
    plot5 = hist(residual,main = "Histogram of Residual")

    # Nice way showing of results
    # Coefficients
    table_betas = matrix(betas,1,n_predictors)
    colnames(table_betas) = c("intercept",col_predictor)
    rownames(table_betas) = "coefficient"
    table_betas = as.table(table_betas)
    # Confidence intervals
    tab_ci = matrix(conf_int,n_predictors,2)
    colnames(tab_ci) = c("lower bound","Upper bound")
    rownames(tab_ci) = c(1:n_predictors)
    tab_ci = as.table(tab_ci)
    # R^2
    tab_R_sq = matrix(Rsq,1,1)
    colnames(tab_R_sq) = "Value"
    rownames(tab_R_sq) = c("R-squared")
    tab_R_sq = as.table(tab_R_sq)
    # C_p
    tab_C_p = matrix(Cp,1,1)
    colnames(tab_C_p) = "Value"
    rownames(tab_C_p) = c("C_p")
    tab_C_p = as.table(tab_C_p)
    # F_star
    tab_F_star = matrix(F_star,1,1)
    colnames(tab_F_star) = "Value"
    rownames(tab_F_star) = c("F-statistic")
    tab_F_star = as.table(tab_F_star)
    # p-value
    tab_p = matrix(p_value,1,1)
    colnames(tab_p) = "Value"
    rownames(tab_p) = c("p-value")
    tab_p = as.table(tab_p)

    result = list(table_betas, tab_ci, tab_R_sq, tab_C_p, tab_F_star, tab_p)

    return(result)
  } # validation bracket
}



