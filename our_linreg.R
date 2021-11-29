# # Import data 
# abalone <- read.csv(file = 'abalone.csv') 
# 
# # Transform dummy variables
# library(fastDummies)
# abalone <- dummy_cols(abalone, select_columns = "Sex")
# 
# # Create response variable and set of predictors 
# response <- abalone$Rings
# explanatory <- abalone[, 2:11]
# explanatory <- subset(explanatory, select = -Rings)



lin_reg = function(explanatory,response,alpha = 0.05){
  y = as.vector(response)
  n = length(y)
  X = cbind(rep(1,n),explanatory)
  X = as.matrix(X)
  ## calculation of beta_hat
  p = dim(X)[2]
  y_mean = mean(y)
  ## initial beta
  beta0 = rep(NA,p)
  beta0[1] = y_mean
  if (p>1){
    for (i in 2 : p){
      beta0[i] = cov(y,X[,i])/var(X[,i])
    }
  }
  
  fnc = function(aa) t(y-X%*%aa)%*%(y-X%*%aa)
  beta_hat = optim(beta0,fnc,method = "L-BFGS-B")$par
  
  y_hat = X%*%beta_hat
  
  sigma_hat_sq = as.numeric(t(y-X%*%beta_hat)%*%(y-X%*%beta_hat)/(n-p))
  SSE = 0
  SST = 0
  for (i in 1 : n){
    SSE = SSE + (y[i] - y_hat[i])^2
    SST = SST + (y[i]-y_mean)^2
  }
  
  ## R^2
  Rsq = 1-(SSE/SST)
  
  # C_p
  Cp = SSE + 2*p*sigma_hat_sq
  
  # confidence intervals
  X_t_X = solve(t(X)%*%X)
  D = rep(NA,p)
  for (i in 1 : p){
    D[i] = X_t_X[i,i]
  }
  ss = qnorm((1-alpha/2),mean=0,sd=1)*sqrt(sigma_hat_sq*D)
  conf_int = cbind(beta_hat-ss,beta_hat+ss)
  
  DFM = p-1
  DFE = n-p
  SSM = 0
  for (i in 1 : n){
    SSM = SSM + (y_hat[i] - y_mean)^2
  }
  MSM = SSM/DFM
  MSE = SSE/DFE
  
  # F_star
  F_star = MSM/MSE
  
  # P-value
  P = 1-pf(F_star,DFM,DFE)
  
  residual = y - y_hat
  
  # Plots
  
  
  plot1 = plot(y_hat,residual)
  plot2 = abline(h=0,col="red")
  plot3 = qqnorm(residual)
  plot4 = qqline(residual, col = 'red',datax = 0)
  
  return(list(beta_hat,conf_int,Rsq,Cp,F_star,P,plot1,plot2,plot3,plot4))
}


