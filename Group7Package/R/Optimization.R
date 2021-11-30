



#' @title Linear Regression Via Gradient Descent Optimization
#'
#' @description Performs linear regression via gradient descent optimization, and outputs the estimated coefficients,
#' confidence intervals for the coefficients, R square, cp, F statistics, and its p-value, plots of residual
#' against fitted values, qq-plot of residuals, and histogram of residuals.
#'
#' @param X An n x p \code{matrix} or a \code{column vector} of length \code{n} whose columns represent the values of
#' p explanatory variables, and \code{n} rows are the observations.
#'
#' @param y A \code{vector} of length n,whose entries are the values of the
#' response variable.
#'
#' @param alpha A \code{number} that lies between \code{0} and \code{1}, which represent the confidence level. If not
#' specified, \code{alpha} is assigned a default value of 0.05.
#' @return Returns a \code{list} containing the following attributes:
#' \describe{
#'      \item{betas}{A vector of estimated value of the coefficients}
#'
#'      \item{CI}{Confidence intervals for the coefficients}
#'
#'      \item{Rsq}{R-square value, which is a measure of goodness of fit of the model}
#'
#'      \item{cp}{cp value, which is a measure of goodness of fit of the model}
#'      \item{F_stat}{F statistic for the model}
#'
#'      \item{p-value}{p-value of the F_Stat}
#'
#'      \item{Resid_plot}{Plot of residuals against fitted values}
#'
#'      \item{qq-plot}{qq-plot of residulas}
#'
#'      \item{hist_resid}{Histogram of residuals}
#' }
#' @author Carlos Matus, Ridvan Ozdemir, and Ogonnaya Romanus
#' @importFrom stats rnorm
#' @export
#' @examples
#' M<-matrix(rnorm(1000),ncol=10) #explanatory variable values
#'  y<-rnorm(100)               #response
#' optimization(y,M,0.01)
#'
#' # Import data
#' abalone <- read.csv(file = 'abalone.csv')

#' # Transform dummy variables
#' library(fastDummies)
#' abalone <- dummy_cols(abalone, select_columns = "Sex")

#' # Create response variable and set of predictors
#' y <- abalone$Rings
#' x <- abalone[, 2:11]
#' x <- subset(x, select = -Rings)
#' optimization(y,x)

optimization<-function(y,x,alpha=0.05){
  y<-as.vector(y) #making the response a vector
  x<-as.matrix(x) #making the exploratory a matrix
  X<-cbind(rep(1,dim(x)[2]),x)  #Adding column vector of 1s to x
  #beta initial
  beta0<-rep(NA,dim(X)[2])
  beta0[1]<-mean(y)
  for (i in 2:length(beta0)) { #This loop will give other components of beta0
    beta0[i]<-cov(y,X[,i])/var(X[,i])}
  #loss
  loss<-function(y,x,beta0){
    beta0<-as.vector(beta0)
    f<- t(y-X%*%beta0)%*%(y-X%*%beta0)

  }
  #minimization
  mini<-optim(beta0,loss,y=y,x=X) #applying optim function to obtain betas that minimize the loss.
  #Note: by default optim function minimizes
  beta_hat<-mini$par  # return the betas.
  #return(beta_hat)

  n<-dim(X)[1]  #number of observations
  p<-dim(X)[2]  #number of predictors
  DFE<- n-p       #degree of freedom
  DFM<- p-1
  y_hat <- X %*%as.matrix(beta_hat)
  y_mean <- mean(y)

  #Residuals
  residu <- y - y_hat

  SSE <- sum(data.matrix(y - y_hat)^2)
  SST <- sum((y-y_mean)^2)

  sigma_hat <- (1/DFE)*as.numeric(t(residu)%*%(residu))

  ## R^2
  Rsq <- 1-(SSE/SST)

  # C_p
  cp <- SSE + 2*p*sigma_hat


  X_t_X = solve(t(X)%*%X)
  D = rep(NA,p) # Initializing the diagonal
  for (i in 1 : p){
    D[i] = X_t_X[i,i] # Filling the diagonal up
  }
  ss <- qnorm((1-alpha/2),mean=0,sd=1)*sqrt(sigma_hat*D) # the right side of the interval
  CI_beta <- cbind(beta_hat-ss,beta_hat+ss) # confidence intervals

#SSM
  SSM <- 0
  for (i in 1 : n){
    SSM <- SSM + (y_hat[i] - y_mean)^2
  }
  MSM <- SSM/DFM
  MSE <- SSE/DFE

  # F_star
  F_star <- MSM/MSE

  # P-value
  p_value <- 1-pf(F_star, DFM, DFE)

  # Plots

  #Return all estimated values
return(list(beta = beta_hat, sigma = sigma_hat, CI = CI_beta, Rsq = Rsq, cp_value = cp,  F_stat=F_star, p_value=p_value ))

}













