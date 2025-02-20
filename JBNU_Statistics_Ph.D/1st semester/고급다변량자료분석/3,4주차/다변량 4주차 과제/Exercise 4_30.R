# Exercise 4.30

## Consider the used-car data. 
## This data has x1(age, measured in years) and 
## x2(as well as the selling price, measured in thousands of dollars) variables.
x1 <- c(1,2,3,3,4,5,6,8,9,11)
x2 <- c(18.95, 19.00, 17.95, 15.54, 14.00, 12.95, 8.94, 7.49, 6.00, 3.99)
used_car <- data.frame(x1=x1, x2=x2)
print(used_car)

### (a) Determine the power transformation λ1 that makes the x1 values approximately normal. Construct the Q-Q plot for the transformed data.
library(MASS)
library(rcompanion)
Box <- boxcox(used_car$x1~1, lambda=seq(-6,6,0.001))
Cox <- data.frame(Box$x, Box$y) # x : lambda의 값, y : log-likelihood의 값         
Cox2 <- Cox[with(Cox, order(-Cox$Box.y)),] # minus log-likelihood의 크기순으로 정렬
Cox2[1,] # minus log-likelihood의 값이 최대가 되는 경우                                
(lambda <- Cox2[1, "Box.x"]) # minus log-likelihood의 값이 최대가 되는 경우의 lambda 값                 
x1_trans <- (used_car$x1^lambda - 1)/lambda # power transformation  
plotNormalHistogram(x1_trans)
qqnorm(x1_trans) ; qqline(x1_trans)

### (b) Determine the power transformation λ2 that makes the x2 values approximately normal. Construct the Q-Q plot for the transformed data.
Box <- boxcox(used_car$x2~1, lambda=seq(-6,6,0.001))
Cox <- data.frame(Box$x, Box$y) # x : lambda의 값, y : log-likelihood의 값         
Cox2 <- Cox[with(Cox, order(-Cox$Box.y)),] # minus log-likelihood의 크기순으로 정렬
Cox2[1,] # minus log-likelihood의 값이 최대가 되는 경우                                
(lambda <- Cox2[1, "Box.x"]) # minus log-likelihood의 값이 최대가 되는 경우의 lambda 값                 
x2_trans <- (used_car$x2^lambda - 1)/lambda # power transformation   
plotNormalHistogram(x2_trans)
qqnorm(x2_trans) ; qqline(x2_trans)

### (c) Determine the power transformations [λ1,λ2] that make the [x1,x2] values jointly using (4-40). Compare the results with those obtained in Parts a and b.
library(car)
attach(used_car)
summary(powerTransform(cbind(x1,x2)~1))
