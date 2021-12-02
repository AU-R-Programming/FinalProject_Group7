# # Import data
abalone <- read.csv(file = 'abalone.csv')
# 
# Transform dummy variables
library(fastDummies)
abalone <- dummy_cols(abalone, select_columns = "Sex")
# 
# Create response variable and set of predictors
response <- abalone$Rings
predictor <- abalone[, 2:11]
predictor <- subset(predictor, select = -Rings)

library(group7project)
lin_reg(response, predictor, iterations = 2000)