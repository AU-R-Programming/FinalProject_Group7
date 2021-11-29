
optimization<-function(y,x){
  y<-as.vector(y) #making the response a vector
  x<-as.matrix(x) #making the exploratory a matrix 
  X<-cbind(rep(1,dim(x)[2]),x)  #Adding column vector of 1s to x
  #beta initial
  b0<-rep(NA,dim(X)[2])
  b0[1]<-mean(y)
  for (i in 2:length(b0)) { #This loop will give other components of beta0 
    b0[i]<-cov(y,X[,i])/var(X[,i])}
  #loss
  loss<-function(y,x,b0){
    b0<-as.vector(b0)
    f<- t(y-X%*%b0)%*%(y-X%*%b0)
    
  }
  #minimization
  mini<-optim(b0,loss,y=y,x=X) #applying optim function to obtain betas that minimize the loss.
  #Note: by default optim function minimizes
  return(mini$par)    # return the betas. 
  
}
# Example
M<-matrix(rnorm(400),ncol=4) #explanatory variable values
y<-rnorm(100)               #response
optimization(y,M)


#R based linear regression
lm(y ~ M)

