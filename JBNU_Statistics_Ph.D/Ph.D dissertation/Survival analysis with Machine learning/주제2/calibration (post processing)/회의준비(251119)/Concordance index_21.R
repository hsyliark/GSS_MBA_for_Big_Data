# reference : https://www.erikdrysdale.com/survConcordance/

## Function for the concordance index (C-Index)

# Function for the i'th concordance
cindex_i <- function(So,eta,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  idx.k <- which(So[,1] > tt_i)
  conc <- sum(eta[i] > eta[idx.k]  )
  disc <- sum(eta[i] < eta[idx.k]  )
  return(c(conc,disc))
}

# Wrapper for total concordance
cindex <- function(So,eta) {
  conc.disc <- c(0,0)
  for (i in which(So[,2] == 1)) {
    conc.disc <- conc.disc + cindex_i(So,eta,i)
  }
  names(conc.disc) <- c('concordant','discordant')
  return(conc.disc)
}

# Example of the 'veteran' data
library(survival)
So <- with(veteran,Surv(time,status==1))
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior),data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
Som <- as.matrix(So)
Som[,1] <- Som[,1] + (1-Som[,2])*min(Som[,1])/2
survConcordance.fit(So, eta)[1:2]
cindex(Som, eta)
concordance(coxph(So ~ X))


## Equation (2)

sigmoid <- function(x) { 1/(1+exp(-x)) }

l_eta_i <- function(So,eta,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  idx.k.i <- which(So[,1] > tt_i)
  loss.i <- sum(1 + log( sigmoid(eta[i] - eta[idx.k.i]) )/log(2) )
  return(loss.i)
}

l_eta <- function(So, eta) {
  loss <- 0
  for (i in which(So[,2] == 1)) {
    loss <- loss + l_eta_i(So,eta,i)
  }
  return(-loss / nrow(So)^2)
}


## Equation (3)

sigmoid2 <- function(x) { 1/(1+exp(x)) }

dl_eta_i <- function(eta,So,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  idx.k <- which(So[,1] > tt_i)
  idx.j <- which(tt_i > So[,1] & So[,2]==1) 
  res.i <- dd_i*sum(sigmoid2(eta[i] - eta[idx.k])) - sum(sigmoid2(eta[idx.j] - eta[i]))
  return(res.i)
}

dl_eta <- function(X,eta,So) {
  grad <- rep(0, ncol(X))
  for (i in seq_along(eta)) {
    grad <- grad + X[i,] * dl_eta_i(eta, So, i)
  }
  grad <- -1 * grad / nrow(X)^2
  return(grad)
}













#--------------------------------------------------------------------------#


# real data (veteran)
library(survival) 
So <- with(veteran, Surv(time,status==1))
veteran <- veteran[order(So[,1]),] # sorting observed time=min(Y,C)
So <- So[order(So[,1]),] # sorting observed time=min(Y,C)
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior), data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5, ncol(X)))
time.sort <- sort(So[,1]) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation

# real data (lung)
library(survival)
lung_edit <- lung
lung_edit$status <- lung_edit$status - 1
lung_edit <- na.omit(lung_edit)
lung_edit <- lung_edit[,-1]
So <- with(lung_edit, Surv(time, status==1))
lung_edit <- lung_edit[order(So[,1]),] # sorting observed time=min(Y,C)
So <- So[order(So[,1]),] # sorting observed time=min(Y,C)
X <- model.matrix(~age+sex+ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss, 
                  data=lung_edit)[,-1]
X <- scale(X)
lung_edit_1 <- as.data.frame(cbind(lung_edit$time, lung_edit$status, X))
library(dplyr)
lung_edit_1 <- rename(lung_edit_1, time=V1, status=V2)
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5, ncol(X)))
time.sort <- sort(So[,1]) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation

1-mean(lung_edit$status)
time.sort0 <- quantile(So[,1])
time.sort1 <- quantile(So[,1], probs=c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                                       0.6, 0.7, 0.8, 0.9))
time.sort2 <- quantile(So[,1], probs=c(0.1, 0.2, 0.5, 0.8, 0.9))
time.sort3 <- quantile(So[,1], probs=c(0.25, 0.5, 0.75))



# Simulation data generation (Assuming Cox PH model)
# Assume parametric basis distribution is Exponential distribution... 
lambda1 <- 0.005 # parameter for generate Y   
lambda2 <- 0.5/lambda1 # parameter for generate C 
beta0 <- 1:5
beta0 <- beta0/sqrt(sum(beta0^2)) # true coefficients
Y <- c() # real survival time
C <- c() # censoring time
t <- c() # observed survival time
censor <- c() # censoring indicator
X <- matrix(, nrow=100, ncol=length(beta0)) # explanatory variables
U1 <- c() # random variable from U(0,1) for generate Y
U2 <- c() # random variable from U(0,1) for generate C 

set.seed(1004)

for (i in 1:100) {
  U1[i] <- runif(1)
  U2[i] <- runif(1)
  X[i,] <- rnorm(5,mean=0,sd=1)
  Y[i] <- -(1/lambda1)*log(U1[i])*exp(-X[i,]%*%beta0)
  C[i] <- lambda2-U2[i]
  t[i] <- min(Y[i],C[i])
  censor[i] <- ifelse(Y[i]<=C[i],1,0)
}

mean(1-censor) # censoring rate

dt5 <- data.frame(X=X, Y=Y, C=C, time=t, censor=censor)
head(dt5,10)

library(survival)
So <- with(dt5, Surv(time,censor==1))
dt5 <- dt5[order(So[,1]),] # sorting observed time=min(Y,C)
So <- So[order(So[,1]),] # sorting observed time=min(Y,C)
X <- model.matrix(~X.1+X.2+X.3+X.4+X.5, data=dt5)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
beta_hat <- model1[["coefficients"]]
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation
censor_time <- dt5$time[dt5$censor==0]
dt_censor <- data.frame(censor_time=censor_time, y=1)

# case1
# censoring rate = 0.13
# lambda1 = 0.005, lambda2 = 4 / lambda1
# set.seed(1004)
time.sort0 <- seq(from=20, to=770, length.out=100) # sorting survival time
time.sort1 <- seq(from=20, to=770, length.out=50) # sorting survival time
time.sort2 <- seq(from=20, to=770, length.out=20) # sorting survival time
time.sort3 <- seq(from=20, to=770, length.out=5) # sorting survival time
time.sort4 <- seq(from=30, to=50, by=2) # sorting survival time
# case2
# censoring rate = 0.29
# lambda1 = 0.005, lambda2 = 1.5 / lambda1
# set.seed(1004)
time.sort0 <- seq(from=20, to=270, length.out=100) # sorting survival time
time.sort1 <- seq(from=20, to=270, length.out=50) # sorting survival time
time.sort2 <- seq(from=20, to=270, length.out=20) # sorting survival time
time.sort3 <- seq(from=20, to=270, length.out=5) # sorting survival time
time.sort4 <- seq(from=30, to=50, by=2) # sorting survival time
# case3
# censoring rate = 0.56
# lambda1 = 0.005, lambda2 = 0.5 / lambda1
# set.seed(1004)
time.sort0 <- seq(from=20, to=80, length.out=100) # sorting survival time
time.sort1 <- seq(from=20, to=80, length.out=50) # sorting survival time
time.sort2 <- seq(from=20, to=80, length.out=20) # sorting survival time
time.sort3 <- seq(from=20, to=80, length.out=5) # sorting survival time
time.sort4 <- seq(from=30, to=50, by=2) # sorting survival time





## C-Index with beta

C_Index_beta <- function(So, beta, X) {
  loss <- 0
  ip <- 0
  for (i in which(So[,2] == 1)) {
    for (j in i:nrow(So)) {
      loss <- loss + ifelse(t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) > 0, 1, 0) 
      ip <- ip + 1
    }
  }
  return(loss / (ip-length(which(So[,2] == 1))))
}

C_Index_t_beta <- function(So, beta, X, t) {
  loss <- 0
  ip <- 0
  for (i in which((So[,2] == 1) & (So[,1] <= t))) {
    for (j in which(So[,1] > t)) {
      loss <- loss + ifelse(t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) > 0, 1, 0) 
      ip <- ip + 1
    }
  }
  return(loss / ip)
}


## Generalization

# U(\beta)
# Need sorting observed time=min(Y,C)

sigmoid <- function(x) { 1/(1+exp(-x)) }

C_tilde_beta <- function(So, beta, X) {
  loss <- 0
  ip <- 0
  for (i in which(So[,2] == 1)) {
    for (j in i:nrow(So)) {
    loss <- loss + 
      (1 + log( sigmoid( t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) ) )/log(2))
    ip <- ip + 1
    }
  }
  return(loss / (ip-length(which(So[,2] == 1))))
}

C_tilde_t_beta <- function(So, beta, X, t) {
  loss <- 0
  ip <- 0
  for (i in which((So[,2] == 1) & (So[,1] <= t))) {
    for (j in which(So[,1] > t)) {
    loss <- loss + 
      (1 + log( sigmoid( t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) ) )/log(2))
    ip <- ip + 1
    }
  }
  return(loss / ip)
}

sum_C_tilde_t_beta <- function(So, beta, X, time.sort) {
  loss <- 0
  for ( k in 2:length(time.sort) ) {
    loss <- loss + 
      (C_tilde_t_beta(So, beta, X, time.sort[k]) - C_tilde_t_beta(So, beta, X, time.sort[k-1]))^2
  }
  return(loss)
}

U_beta <- function(So, beta, X, lambda, time.sort) {
  loss_U <- -C_tilde_beta(So, beta, X) + 
    lambda*sum_C_tilde_t_beta(So, beta, X, time.sort)
  return(loss_U)
}


# Take the derivative U(\beta) by \beta
# Need sorting observed time=min(Y,C)

sigmoid <- function(x) { 1/(1+exp(-x)) }

der_C_tilde_beta <- function(So, beta, X) {
  loss <- 0
  ip <- 0
  for (i in which(So[,2] == 1)) {
    for (j in i:nrow(So)) {
    loss <- loss - 
      (as.matrix(X[j,]) - as.matrix(X[i,]))%*%sigmoid( t(beta)%*%(as.matrix(X[j,]) - as.matrix(X[i,])) )
    ip <- ip + 1
    }
  }
  return(loss / (ip-length(which(So[,2] == 1))) / log(2))
}

der_C_tilde_t_beta <- function(So, beta, X, t) {
  loss <- 0
  ip <- 0
  for ( i in which((So[,2] == 1) & (So[,1] <= t)) ) {
    for (j in which(So[,1] > t)) {
    loss <- loss - 
      (as.matrix(X[j,]) - as.matrix(X[i,]))%*%sigmoid( t(beta)%*%(as.matrix(X[j,]) - as.matrix(X[i,])) )
    ip <- ip + 1
    }
  }
  return(loss / ip / log(2))
}

sum_der_C_tilde_t_beta <- function(So, beta, X, time.sort) {
  loss <- 0
  for ( k in 2:length(time.sort) ) {
    loss <- loss + 
      as.numeric(C_tilde_t_beta(So, beta, X, time.sort[k]) - C_tilde_t_beta(So, beta, X, time.sort[k-1]))*
      (der_C_tilde_t_beta(So, beta, X, time.sort[k]) - der_C_tilde_t_beta(So, beta, X, time.sort[k-1]))
  }
  return(loss)
}

der_U_beta <- function(So, beta, X, lambda, time.sort) {
  der_loss_Q <- -der_C_tilde_beta(So, beta, X) + 
    2*lambda*sum_der_C_tilde_t_beta(So, beta, X, time.sort)
  return(der_loss_Q)
}


# Gradient descent algorithm

my.gradient.descent.U <- function(So, X, alpha, lambda, time.sort) {
  
  # alpha : learning rate
  # lambda : penalty parameter
  
  # if(lambda < 0) 
    # stop("Lambda is negative value. Please insert non-negative value of lambda.")
  
  X <-as.matrix(X)
  
  beta.old <- rep(0.5,ncol(X)) # The initial value of coefficient vector
  beta.old <- beta.old/sqrt(sum(beta.old^2))
  iteration <- c()
  difference <- c()
  c_index <- c()
  c_tilde <- c()
  sum_C_tilde_t <- c()
  
  for (i in 1:10000) {
    
    beta.new <- as.vector(beta.old - alpha*der_U_beta(So, beta.old, X, lambda, time.sort))
    diff <- sqrt(sum((beta.new - beta.old)^2))
    Cindex <- C_Index_beta(So, beta.new, X)
    Ctilde <- C_tilde_beta(So, beta.new, X)
    Ctil <- sum_C_tilde_t_beta(So, beta.new, X, time.sort)
    
    iteration <- c(iteration, i)
    difference <- c(difference, diff)
    c_index <- c(c_index, Cindex)
    c_tilde <- c(c_tilde, Ctilde)
    sum_C_tilde_t <- c(sum_C_tilde_t, Ctil)
    
    cat("( iteration , difference, C_Index, C_tilde, sum_C_tilde_t ) 
        = (", i, ",", diff, ",", Cindex, ",", Ctilde, ",", Ctil, ")\n")
    
    if (diff < 1E-4) break
    
    beta.old <- beta.new
    
  }
  
  c_index_t <- c()
  c_tilde_t <- c()
  
  for (j in 1:length(So[,1])) {
    c_index_t[j] <- C_Index_t_beta(So, beta.new, X, So[,1][j])
    c_tilde_t[j] <- C_tilde_t_beta(So, beta.new, X, So[,1][j])
  }
  
  cat("Algorithm converged...","\n\n")
  
  return(list(beta.new=beta.new, iteration=iteration, difference=difference,
              c_index=c_index, c_tilde=c_tilde, sum_C_tilde_t=sum_C_tilde_t,
              time=So[,1], c_index_t=c_index_t, c_tilde_t=c_tilde_t))
  
}  


# Mini-batch(Stochastic) gradient descent algorithm
# reference : https://jjeongil.tistory.com/577 (stochastic)
# reference : https://data-science-hi.tistory.com/164 (mini-batch)
# reference : https://brunch.co.kr/@linecard/561 (mini-batch)
# reference : https://sonny-daily-story.tistory.com/5 (mini-batch)
# reference : https://www.youtube.com/watch?v=87Q2LIlMWoY (mini-batch)
# reference : https://light-tree.tistory.com/133 (mini-batch)
# reference : https://89douner.tistory.com/43 (mini-batch)

my.mini.gradient.U <- function(So, X, alpha, lambda, k, time.sort) {
  
  # alpha : learning rate
  # lambda : penalty parameter
  # k : number of batch
  
  # if(lambda < 0) 
    # stop("Lambda is negative value. Please insert non-negative value of lambda.")
  
  X <-as.matrix(X)
  
  beta.old <- rep(0.5,ncol(X)) # The initial value of coefficient vector
  beta.old <- beta.old/sqrt(sum(beta.old^2))
  iteration <- c()
  difference <- c()
  c_index <- c()
  c_tilde <- c()
  sum_C_tilde_t <- c()
  
  for (i in 1:10000) {
    
    idx <- sample(1:nrow(X), nrow(X), replace=F)
    grad <- rep(0,ncol(X))
    
    for (j in 0:(k-1)) {
      batch.idx <- idx[(1:nrow(X))%%k==j]
      batch.So <- So[batch.idx,]
      batch.X <- X[batch.idx,]
      
      grad <- grad + der_U_beta(batch.So, beta.old, batch.X, lambda, time.sort)
      
      cat("( batch ) = ", j+1, "\n")
    }
    
    mean_grad <- grad/k
      
    beta.new <- beta.old - alpha*mean_grad
    diff <- sqrt(sum((beta.new - beta.old)^2))  
    Cindex <- C_Index_beta(So, beta.new, X)
    Ctilde <- C_tilde_beta(So, beta.new, X)
    Ctil <- sum_C_tilde_t_beta(So, beta.new, X, time.sort)
    
    iteration <- c(iteration, i)
    difference <- c(difference, diff)
    c_index <- c(c_index, Cindex)
    c_tilde <- c(c_tilde, Ctilde)
    sum_C_tilde_t <- c(sum_C_tilde_t, Ctil)
    
    cat("( iteration , difference, C_Index, C_tilde, sum_C_tilde_t ) 
        = (", i, ",", diff, ",", Cindex, ",", Ctilde, ",", Ctil, ")\n")
      
    if (diff < 1E-4) break
      
    beta.old <- beta.new
    }
    

  cat("Algorithm converged...","\n\n")
  
  return(list(beta.new=beta.new, iteration=iteration, difference=difference,
              c_index=c_index, c_tilde=c_tilde, sum_C_tilde_t=sum_C_tilde_t))
  
}  


res1 <- my.mini.gradient.U(So, X, alpha=0.2, lambda=3, k=5, time.sort1)
res2 <- my.mini.gradient.U(So, X, alpha=0.2, lambda=3, k=5, time.sort2)

beta_hat
C_Index_beta(So, beta_hat, X)
C_tilde_beta(So, beta_hat, X)

# real data (veteran)
time.sort1 <- quantile(So[,1], probs=c(0.05, 0.2, 0.3, 0.5, 0.7))
time.sort2 <- quantile(So[,1], probs=c(0.1, 0.2, 0.5, 0.8))
time.sort3 <- quantile(So[,1], probs=c(0.25, 0.5, 0.75))

risk.score_cox <- predict(object=coxph(So ~ X), newdata=veteran, type="risk")

res1 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=0, time.sort1)
res1$beta.new
risk.score_01 <- exp(X%*%res1$beta.new)
res2 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=5, time.sort1)
res2$beta.new
risk.score_02 <- exp(X%*%res2$beta.new)
res3 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=10, time.sort1)
res3$beta.new
risk.score_03 <- exp(X%*%res3$beta.new)
res4 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=20, time.sort1)
res4$beta.new
risk.score_04 <- exp(X%*%res4$beta.new)

res5 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=0, time.sort2)
res5$beta.new
risk.score_05 <- exp(X%*%res5$beta.new)
res6 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=5, time.sort2)
res6$beta.new
risk.score_06 <- exp(X%*%res6$beta.new) 
res7 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=10, time.sort2)
res7$beta.new
risk.score_07 <- exp(X%*%res7$beta.new)
res8 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=20, time.sort2)
res8$beta.new
risk.score_08 <- exp(X%*%res8$beta.new)

res9 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=0, time.sort3)
res9$beta.new
risk.score_09 <- exp(X%*%res9$beta.new)
res10 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=5, time.sort3)
res10$beta.new
risk.score_10 <- exp(X%*%res10$beta.new)
res11 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=10, time.sort3)
res11$beta.new
risk.score_11 <- exp(X%*%res11$beta.new)
res12 <- my.gradient.descent.U(So, X, alpha=0.001, lambda=20, time.sort3)
res12$beta.new
risk.score_12 <- exp(X%*%res12$beta.new)


# Select time.sort
library(timeROC)

risk.score_cox <- predict(object=coxph(So ~ X), newdata=lung_edit_1, type="risk")
coxph(So ~ X)$coefficients
C_Index_beta(So, coxph(So ~ X)$coefficients, X)
C_tilde_beta(So, coxph(So ~ X)$coefficients, X)
ROC.00 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_cox,
                  cause=1,weighting="marginal",times=So[,1])
case00 <- ROC.00[["AUC"]]
mean(case00, na.rm=T)
sd(case00, na.rm=T)

res1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0, time.sort2)
res1$beta.new
risk.score_01 <- exp(X%*%res1$beta.new)
C_Index_beta(So, res1$beta.new, X)
C_tilde_beta(So, res1$beta.new, X)
ROC.01 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_01,
                  cause=1,weighting="marginal",times=So[,1])
case01 <- ROC.01[["AUC"]]
mean(case01, na.rm=T)
sd(case01, na.rm=T)

res2 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=50, time.sort2)
res2$beta.new
risk.score_02 <- exp(X%*%res2$beta.new)
C_Index_beta(So, res2$beta.new, X)
C_tilde_beta(So, res2$beta.new, X)
ROC.02 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_02,
                  cause=1,weighting="marginal",times=So[,1])
case02 <- ROC.02[["AUC"]]
mean(case02, na.rm=T)
sd(case02, na.rm=T)

res3 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=100, time.sort2)
res3$beta.new
risk.score_03 <- exp(X%*%res3$beta.new)
C_Index_beta(So, res3$beta.new, X)
C_tilde_beta(So, res3$beta.new, X)
ROC.03 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_03,
                  cause=1,weighting="marginal",times=So[,1])
case03 <- ROC.03[["AUC"]]
mean(case03, na.rm=T)
sd(case03, na.rm=T)

res4 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=0, time.sort2)
res4$beta.new
risk.score_04 <- exp(X%*%res4$beta.new)
C_Index_beta(So, res4$beta.new, X)
C_tilde_beta(So, res4$beta.new, X)
ROC.04 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_04,
                  cause=1,weighting="marginal",times=So[,1])
case04 <- ROC.04[["AUC"]]
mean(case04, na.rm=T)
sd(case04, na.rm=T)

res5 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=5, time.sort2)
res5$beta.new
risk.score_05 <- exp(X%*%res5$beta.new)
C_Index_beta(So, res5$beta.new, X)
C_tilde_beta(So, res5$beta.new, X)
ROC.05 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_05,
                  cause=1,weighting="marginal",times=So[,1])
case05 <- ROC.05[["AUC"]]
mean(case05, na.rm=T)
sd(case05, na.rm=T)

res6 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=10, time.sort2)
res6$beta.new
risk.score_06 <- exp(X%*%res6$beta.new)
C_Index_beta(So, res6$beta.new, X)
C_tilde_beta(So, res6$beta.new, X)
ROC.06 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_06,
                  cause=1,weighting="marginal",times=So[,1])
case06 <- ROC.06[["AUC"]]
mean(case06, na.rm=T)
sd(case06, na.rm=T)

library(reshape)
dt <- data.frame(time=ROC.00[["times"]], censor=lung_edit_1$status, case00=case00,
                 case01=case01, case02=case02, case03=case03, case04=case04,
                 case05=case05, case06=case06)
dt <- na.omit(dt)
dt1 <- melt(data = dt,
            id.vars = c("time", "censor"),
            measure.vars = c("case00", "case01", "case02",
                             "case03", "case04", "case05",
                             "case06"))
colnames(dt1) <- c("time", "censor", "case", "AUC")
dt2 <- na.omit(dt1)
aggregate(dt2$AUC, list(dt2$case), FUN=mean, na.action = na.omit)

library(ggplot2)
ggplot(dt2, aes(x=time, y=AUC, group=case, color=case)) +
  geom_point(aes(x=time, y=1, color=as.factor(censor))) +
  geom_line(aes(group=case, color=case, linetype=case)) +
  geom_vline(xintercept=time.sort2, col="darkgreen", linetype="dashed") +
  ggtitle("Time dependent AUC with risk score (gradient descent) : 'lung' data") +
  xlab("time") + ylab("AUC") +
  #scale_x_continuous(limits = c(0, 200)) +
  #scale_y_continuous(limits = c(0.5, 0.8)) +
  theme(plot.title = element_text(hjust = 0.5))






















# change lambda
risk.score_cox <- predict(object=coxph(So ~ X), newdata=dt5, type="risk")
res1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0, time.sort3)
res1$beta.new
risk.score_01 <- exp(X%*%res1$beta.new)
res2 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=50, time.sort3)
res2$beta.new
risk.score_02 <- exp(X%*%res2$beta.new)
res3 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=100, time.sort3)
res3$beta.new
risk.score_03 <- exp(X%*%res3$beta.new)

# change number of sorting time point
risk.score_cox <- predict(object=coxph(So ~ X), newdata=dt5, type="risk")
res1 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=3, time.sort0)
res1$beta.new
risk.score_01 <- exp(X%*%res1$beta.new)
res2 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=3, time.sort1)
res2$beta.new
risk.score_02 <- exp(X%*%res2$beta.new)
res3 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=3, time.sort2)
res3$beta.new
risk.score_03 <- exp(X%*%res3$beta.new)
res4 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=3, time.sort3)
res4$beta.new
risk.score_04 <- exp(X%*%res4$beta.new)

# change lambda with restricted time point area
risk.score_cox <- predict(object=coxph(So ~ X), newdata=dt5, type="risk")
res1 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=0, time.sort4)
res1$beta.new
risk.score_01 <- exp(X%*%res1$beta.new)
res2 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=3, time.sort4)
res2$beta.new
risk.score_02 <- exp(X%*%res2$beta.new)
res3 <- my.gradient.descent.U(So, X, alpha=0.2, lambda=5, time.sort4)
res3$beta.new
risk.score_03 <- exp(X%*%res3$beta.new)

library(ggplot2)
dat <- data.frame(iteration=res3$iteration, difference=res3$difference,
                  c_index=res3$c_index, c_tilde=res3$c_tilde, 
                  sum_C_tilde_t=res3$sum_C_tilde_t)
ggplot(data=dat, aes(x=iteration, y=sum_C_tilde_t)) +
  geom_line() +
  geom_point() +
  ggtitle('gradient descent algorithm') +
  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold'))

diff_c_index_t <- c()
diff_c_tilde_t <- c()
for (i in 2:length(res3$time)) {
  diff_c_index_t[i] <- res3$c_index_t[i] - res3$c_index_t[i-1]
  diff_c_tilde_t[i] <- res3$c_tilde_t[i] - res3$c_tilde_t[i-1]
}
diff_c_index_t_2 <- diff_c_index_t^2
diff_c_tilde_t_2 <- diff_c_tilde_t^2

dat1 <- data.frame(time=res3$time, c_index_t=res3$c_index_t,
                  c_tilde_t=res3$c_tilde_t)
ggplot(data=dat1, aes(x=time, y=c_tilde_t)) +
  geom_line() +
  geom_point() +
  ggtitle('gradient descent algorithm') +
  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold'))

dat2 <- data.frame(diff_c_index_t=diff_c_index_t,
                   diff_c_tilde_t=diff_c_tilde_t,
                   diff_c_index_t_2=diff_c_index_t_2,
                   diff_c_tilde_t_2=diff_c_tilde_t_2,
                   time=res3$time)

ggplot(data=dat2, aes(x=diff_c_index_t, y=diff_c_tilde_t)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept=0, slope=1, colour="red") +
  ggtitle('gradient descent algorithm') +
  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold'))

ggplot(data=dat2) +
  geom_line(aes(x=time, y=diff_c_index_t, color="diff_c_index_t")) +
  geom_line(aes(x=time, y=diff_c_tilde_t, color="diff_c_tilde_t")) +
  ggtitle('gradient descent algorithm') +
  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold'))

ggplot(data=dat2, aes(x=diff_c_index_t_2, y=diff_c_tilde_t_2)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept=0, slope=1, colour="red") +
  ggtitle('gradient descent algorithm') +
  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold'))

ggplot(data=dat2) +
  geom_line(aes(x=time, y=diff_c_index_t_2, color="diff_c_index_t_2")) +
  geom_line(aes(x=time, y=diff_c_tilde_t_2, color="diff_c_tilde_t_2")) +
  ggtitle('gradient descent algorithm') +
  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold'))




#---------------------------------------------------------------------------------------#

library(survival)
So <- with(veteran, Surv(time,status==1))
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior), data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
risk.score_cox <- predict(object=coxph(So ~ X), newdata=veteran, type="risk") # risk score
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5, ncol(X)))
time.sort <- sort(So[,1]) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation


# alpha=0.01, lambda=0, k=5, number of sorting time=2 (case01)
time.sort_01 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res01 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0, k=5, time.sort_01)
risk.score_01 <- exp(X%*%res01$beta.new)
res01_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0, time.sort_01)
risk.score_01_1 <- exp(X%*%res01_1$beta.new)

# alpha=0.01, lambda=0, k=5, number of sorting time=10 (case02)
time.sort_02 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res02 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0, k=5, time.sort_02)
risk.score_02 <- exp(X%*%res02$beta.new)
res02_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0, time.sort_02)
risk.score_02_1 <- exp(X%*%res02_1$beta.new)

# alpha=0.01, lambda=0, k=5, number of sorting time=100 (case03)
time.sort_03 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res03 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0, k=5, time.sort_03)
risk.score_03 <- exp(X%*%res03$beta.new)
res03_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0, time.sort_03)
risk.score_03_1 <- exp(X%*%res03_1$beta.new)


# alpha=0.01, lambda=0.1, k=5, number of sorting time=2 (case04)
time.sort_04 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res04 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0.1, k=5, time.sort_04)
risk.score_04 <- exp(X%*%res04$beta.new)
res04_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0.1, time.sort_04)
risk.score_04_1 <- exp(X%*%res04_1$beta.new)

# alpha=0.01, lambda=0.1, k=5, number of sorting time=10 (case05)
time.sort_05 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res05 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0.1, k=5, time.sort_05)
risk.score_05 <- exp(X%*%res05$beta.new)
res05_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0.1, time.sort_05)
risk.score_05_1 <- exp(X%*%res05_1$beta.new)

# alpha=0.01, lambda=0.1, k=5, number of sorting time=100 (case06)
time.sort_06 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res06 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0.1, k=5, time.sort_06)
risk.score_06 <- exp(X%*%res06$beta.new)
res06_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0.1, time.sort_06)
risk.score_06_1 <- exp(X%*%res06_1$beta.new)


# alpha=0.01, lambda=1, k=5, number of sorting time=2 (case07)
time.sort_07 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res07 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=1, k=5, time.sort_07)
risk.score_07 <- exp(X%*%res07$beta.new)
res07_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=1, time.sort_07)
risk.score_07_1 <- exp(X%*%res07_1$beta.new)

# alpha=0.01, lambda=1, k=5, number of sorting time=10 (case08)
time.sort_08 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res08 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=1, k=5, time.sort_08)
risk.score_08 <- exp(X%*%res08$beta.new)
res08_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=1, time.sort_08)
risk.score_08_1 <- exp(X%*%res08_1$beta.new)

# alpha=0.01, lambda=1, k=5, number of sorting time=100 (case09)
time.sort_09 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res09 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=1, k=5, time.sort_09)
risk.score_09 <- exp(X%*%res09$beta.new)
res09_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=1, time.sort_09)
risk.score_09_1 <- exp(X%*%res09_1$beta.new)


# alpha=0.01, lambda=10, k=5, number of sorting time=2 (case10)
time.sort_10 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res10 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=10, k=5, time.sort_10)
risk.score_10 <- exp(X%*%res10$beta.new)
res10_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=10, time.sort_10)
risk.score_10_1 <- exp(X%*%res10_1$beta.new)

# alpha=0.01, lambda=10, k=5, number of sorting time=10 (case11)
time.sort_11 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res11 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=10, k=5, time.sort_11)
risk.score_11 <- exp(X%*%res11$beta.new)
res11_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=10, time.sort_11)
risk.score_11_1 <- exp(X%*%res11_1$beta.new)

# alpha=0.01, lambda=10, k=5, number of sorting time=100 (case12)
time.sort_12 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res12 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=10, k=5, time.sort_12)
risk.score_12 <- exp(X%*%res12$beta.new)
res12_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=10, time.sort_12)
risk.score_12_1 <- exp(X%*%res12_1$beta.new)


## Time-dependent AUC


# Loading packages

library(survival) 
library(timeROC)
library(nsROC)
library(timereg)
library(ggplot2)
library(reshape)
library(plyr)
library(tidyverse)
library(coxed)
library(gridExtra)
library(mgcv)


# Time dependent AUC (Cumulative/Dynamic time-dependent ROC curve)

# gradient descent

case00 <- c()
case01 <- c()
case02 <- c()
case03 <- c()
case04 <- c()
case05 <- c()
case06 <- c()
case07 <- c()
case08 <- c()
case09 <- c()
case10 <- c()
case11 <- c()
case12 <- c()
time.set <- seq(1,999,length.out=20)
for (i in time.set) {
  case00[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_cox,
                     predict.time=i)[["auc"]]
  case01[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_01,
                     predict.time=i)[["auc"]]
  case02[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_02,
                     predict.time=i)[["auc"]]
  case03[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_03,
                     predict.time=i)[["auc"]]
  case04[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_04,
                     predict.time=i)[["auc"]]
  case05[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_05,
                     predict.time=i)[["auc"]]
  case06[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_06,
                     predict.time=i)[["auc"]]
  case07[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_07,
                     predict.time=i)[["auc"]]
  case08[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_08,
                     predict.time=i)[["auc"]]
  case09[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_09,
                     predict.time=i)[["auc"]]
  case10[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_10,
                     predict.time=i)[["auc"]]
  case11[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_11,
                     predict.time=i)[["auc"]]
  case12[i] <- cdROC(stime=So[,1],status=So[,2],marker=risk.score_12,
                     predict.time=i)[["auc"]]
}
case00 <- case00[!is.na(case00)]
case01 <- case01[!is.na(case01)]
case02 <- case02[!is.na(case02)]
case03 <- case03[!is.na(case03)]
case04 <- case04[!is.na(case04)]
case05 <- case05[!is.na(case05)]
case06 <- case06[!is.na(case06)]
case07 <- case07[!is.na(case07)]
case08 <- case08[!is.na(case08)]
case09 <- case09[!is.na(case09)]
case10 <- case10[!is.na(case10)]
case11 <- case11[!is.na(case11)]
case12 <- case12[!is.na(case12)]



ROC.00 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_cox,
                  cause=1,weighting="marginal",times=So[,1])
case00 <- ROC.00[["AUC"]]

ROC.td01 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_01,
                    cause=1,weighting="marginal",times=So[,1])
case01 <- ROC.td01[["AUC"]]

ROC.td02 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_02,
                    cause=1,weighting="marginal",times=So[,1])
case02 <- ROC.td02[["AUC"]]

ROC.td03 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_03,
                    cause=1,weighting="marginal",times=So[,1])
case03 <- ROC.td03[["AUC"]]

ROC.td04 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_04,
                    cause=1,weighting="marginal",times=So[,1])
case04 <- ROC.td04[["AUC"]]

ROC.td05 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_05,
                    cause=1,weighting="marginal",times=So[,1])
case05 <- ROC.td05[["AUC"]]

ROC.td06 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_06,
                    cause=1,weighting="marginal",times=So[,1])
case06 <- ROC.td06[["AUC"]]

ROC.td07 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_07,
                    cause=1,weighting="marginal",times=So[,1])
case07 <- ROC.td07[["AUC"]]

ROC.td08 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_08,
                    cause=1,weighting="marginal",times=So[,1])
case08 <- ROC.td08[["AUC"]]

ROC.td09 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_09,
                    cause=1,weighting="marginal",times=So[,1])
case09 <- ROC.td09[["AUC"]]

ROC.td10 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_10,
                    cause=1,weighting="marginal",times=So[,1])
case10 <- ROC.td10[["AUC"]]

ROC.td11 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_11,
                    cause=1,weighting="marginal",times=So[,1])
case11 <- ROC.td11[["AUC"]]

ROC.td12 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_12,
                    cause=1,weighting="marginal",times=So[,1])
case12 <- ROC.td12[["AUC"]]


dt <- data.frame(time=ROC.00[["times"]], censor=dt5$censor, case00=case00,
                 case01=case01, case02=case02, case03=case03, case04=case04,
                 case05=case05, case06=case06, case07=case07, case08=case08,
                 case09=case09, case10=case10, case11=case11, case12=case12)

dt <- data.frame(time=veteran$time, censor=veteran$status, case00=case00,
                 case01=case01, case02=case02, case03=case03, case04=case04,
                 case05=case05, case06=case06, case07=case07, case08=case08,
                 case09=case09, case10=case10, case11=case11, case12=case12)

dt <- data.frame(time=time.set, case00=case00,
                 case01=case01, case02=case02, case03=case03, case04=case04,
                 case05=case05, case06=case06, case07=case07, case08=case08,
                 case09=case09, case10=case10, case11=case11, case12=case12)

dt <- data.frame(time=1:999, case00=case00,
                 case01=case01, case02=case02, case03=case03, case04=case04,
                 case05=case05, case06=case06, case07=case07, case08=case08,
                 case09=case09, case10=case10, case11=case11, case12=case12)
dt$censor <- ifelse(dt$time == 25 | dt$time == 83
                    | dt$time == 87 | dt$time == 97
                    | dt$time == 100 | dt$time == 103
                    | dt$time == 123 | dt$ time == 182
                    | dt$time == 231, 0, 1)
dt <- na.omit(dt)
dt1 <- melt(data = dt,
            id.vars = c("time", "censor"),
            measure.vars = c("case00", "case01", "case02",
                             "case03", "case04", "case05",
                             "case06", "case07", "case08",
                             "case09", "case10", "case11",
                             "case12"))
colnames(dt1) <- c("time", "censor", "case", "AUC")
dt2 <- na.omit(dt1)
aggregate(dt2$AUC, list(dt2$case), FUN=mean, na.action = na.omit)

dt1 <- data.frame(time=rep(ROC.00[["times"]],13),
                  censor=rep(dt5$censor,13),
                  case=c(rep("case00",137), rep("case01",137), rep("case02",137),
                         rep("case03",137), rep("case04",137),
                         rep("case05",137), rep("case06",137),
                         rep("case07",137), rep("case08",137),
                         rep("case09",137), rep("case10",137),
                         rep("case11",137), rep("case12",137)),
                  AUC=c(case00, case01, case02, case03, case04, case05, case06,
                        case07, case08, case09, case10, case11, case12))

dt1 <- data.frame(time=rep(veteran$time,13),
                  censor=rep(veteran$status,13),
                  case=c(rep("case00",137), rep("case01",137), rep("case02",137),
                         rep("case03",137), rep("case04",137),
                         rep("case05",137), rep("case06",137),
                         rep("case07",137), rep("case08",137),
                         rep("case09",137), rep("case10",137),
                         rep("case11",137), rep("case12",137)),
                  AUC=c(case00, case01, case02, case03, case04, case05, case06,
                        case07, case08, case09, case10, case11, case12))

dt1 <- data.frame(time=rep(time.set,13),
                  case=c(rep("case00",20), rep("case01",20), rep("case02",20),
                         rep("case03",20), rep("case04",20),
                         rep("case05",20), rep("case06",20),
                         rep("case07",20), rep("case08",20),
                         rep("case09",20), rep("case10",20),
                         rep("case11",20), rep("case12",20)),
                  AUC=c(case00, case01, case02, case03, case04, case05, case06,
                        case07, case08, case09, case10, case11, case12))

dt2 <- na.omit(dt1)

aggregate(dt2$AUC, list(dt2$case), FUN=mean, na.action = na.omit)

ggplot(dt2_12, aes(x=time, y=AUC, group=case, color=case)) +
  geom_point(aes(x=time, y=1, color=as.factor(censor))) +
  geom_line(aes(group=case, color=case, linetype=case)) +
  geom_vline(xintercept=time.sort3, col="darkgreen", linetype="dashed") +
  ggtitle("Time dependent AUC with risk score (gradient descent)") +
  xlab("time") + ylab("AUC") +
  #scale_x_continuous(limits = c(0, 200)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(dt2_12, aes(x=time, y=AUC, group=case, color=case)) +
  geom_point(aes(x=time, y=1, color="green")) +
  geom_line(aes(group=case, color=case, linetype=case)) +
  geom_vline(xintercept=time.sort3, col="darkgreen", linetype="dashed") +
  ggtitle("Time dependent AUC with risk score (gradient descent)") +
  xlab("time") + ylab("AUC") +
  scale_x_continuous(limits = c(0, 200)) +
  theme(plot.title = element_text(hjust = 0.5))

dt2_1 <- dt2[dt2$case=="case00" | dt2$case=="case01" | 
               dt2$case=="case02" | dt2$case=="case03" |
               dt2$case=="case04",] 
dt2_01 <- dt2[dt2$case=="case00" | dt2$case=="case01",]
dt2_02 <- dt2[dt2$case=="case01" | dt2$case=="case02",]
dt2_03 <- dt2[dt2$case=="case01" | dt2$case=="case03",]
dt2_04 <- dt2[dt2$case=="case01" | dt2$case=="case04",]

dt2_2 <- dt2[dt2$case=="case00" | dt2$case=="case05" | 
               dt2$case=="case06" | dt2$case=="case07" |
               dt2$case=="case08",]
dt2_05 <- dt2[dt2$case=="case00" | dt2$case=="case05",]
dt2_06 <- dt2[dt2$case=="case05" | dt2$case=="case06",]
dt2_07 <- dt2[dt2$case=="case05" | dt2$case=="case07",]
dt2_08 <- dt2[dt2$case=="case05" | dt2$case=="case08",]

dt2_3 <- dt2[dt2$case=="case00" | dt2$case=="case09" | 
               dt2$case=="case10" | dt2$case=="case11" |
               dt2$case=="case12",]
dt2_09 <- dt2[dt2$case=="case00" | dt2$case=="case09",]
dt2_10 <- dt2[dt2$case=="case09" | dt2$case=="case10",]
dt2_11 <- dt2[dt2$case=="case09" | dt2$case=="case11",]
dt2_12 <- dt2[dt2$case=="case09" | dt2$case=="case12",]



















# gradient descent

ROC.td00_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_cox,
                      cause=1,weighting="marginal",times=So[,1])
case00_1 <- ROC.td00_1[["AUC"]]

ROC.td01_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_01_1,
                      cause=1,weighting="marginal",times=So[,1])
case01_1 <- ROC.td01_1[["AUC"]]

ROC.td02_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_02_1,
                      cause=1,weighting="marginal",times=So[,1])
case02_1 <- ROC.td02_1[["AUC"]]

ROC.td03_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_03_1,
                      cause=1,weighting="marginal",times=So[,1])
case03_1 <- ROC.td03_1[["AUC"]]

ROC.td04_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_04_1,
                      cause=1,weighting="marginal",times=So[,1])
case04_1 <- ROC.td04_1[["AUC"]]

ROC.td05_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_05_1,
                      cause=1,weighting="marginal",times=So[,1])
case05_1 <- ROC.td05_1[["AUC"]]

ROC.td06_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_06_1,
                      cause=1,weighting="marginal",times=So[,1])
case06_1 <- ROC.td06_1[["AUC"]]

ROC.td07_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_07_1,
                      cause=1,weighting="marginal",times=So[,1])
case07_1 <- ROC.td07_1[["AUC"]]

ROC.td08_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_08_1,
                      cause=1,weighting="marginal",times=So[,1])
case08_1 <- ROC.td08_1[["AUC"]]

ROC.td09_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_09_1,
                      cause=1,weighting="marginal",times=So[,1])
case09_1 <- ROC.td09_1[["AUC"]]

ROC.td10_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_10_1,
                      cause=1,weighting="marginal",times=So[,1])
case10_1 <- ROC.td10_1[["AUC"]]

ROC.td11_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_11_1,
                      cause=1,weighting="marginal",times=So[,1])
case11_1 <- ROC.td11_1[["AUC"]]

ROC.td12_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_12_1,
                      cause=1,weighting="marginal",times=So[,1])
case12_1 <- ROC.td12_1[["AUC"]]

dt_1 <- data.frame(time=ROC.td00_1[["times"]], case00_1=case00_1,
                   case01=case01_1, case02=case02_1, case03=case03_1, case04=case04_1,
                   case05=case05_1, case06=case06_1, case07=case07_1, case08=case08_1,
                   case09=case09_1, case10=case10_1, case11=case11_1, case12=case12_1)
dt1_1 <- data.frame(time=rep(ROC.td00_1[["times"]],13),
                    case=c(rep("case00_1",137), rep("case01_1",137), rep("case02_1",137),
                           rep("case03_1",137), rep("case04_1",137),
                           rep("case05_1",137), rep("case06_1",137),
                           rep("case07_1",137), rep("case08_1",137),
                           rep("case09_1",137), rep("case10_1",137),
                           rep("case11_1",137), rep("case12_1",137)),
                    AUC=c(case00_1, case01_1, case02_1, case03_1, 
                          case04_1, case05_1, case06_1, case07_1, 
                          case08_1, case09_1, case10_1, case11_1, case12_1))

dt2_1 <- na.omit(dt1_1)

aggregate(dt2_1$AUC, list(dt2_1$case), FUN=mean, na.action = na.omit)

ggplot(dt1_1, aes(x=time, y=AUC, group=case, color=case)) +
  geom_line() +
  ggtitle("Time dependent AUC with risk score (gradient descent)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))


# ROC curve

TP1 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=11"]
FP1 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=11"]
TP2 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=25"]
FP2 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=25"]
TP3 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=36"]
FP3 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=36"]
TP4 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=411"]
FP4 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=411"]
TP5 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=587"]
FP5 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=587"]
TP <- c(TP1,TP2,TP3,TP4,TP5)
FP <- c(FP1,FP2,FP3,FP4,FP5)

time.survival <- c(rep("t=11",138),rep("t=25",552),rep("t=36",138),
                   rep("t=411",138),rep("t=587",138))
dt.ROC <- data.frame(times=time.survival,TP=TP,FP=FP)

ggplot(data=dt.ROC)+
  geom_line(mapping=aes(x=FP,
                        y=TP,
                        group=times,
                        color=times),
            size=1) +
  ggtitle("Time dependent ROC curve (veteran data)") +
  xlab("False Positive Rate") + ylab("True Positive Rate") +
  theme(plot.title = element_text(hjust = 0.5))


# Make AUC data

# uncensored
uncensored.density <- density(So[,1][So[,2]==1])
data5_1 <- data.frame(time=uncensored.density[["x"]], 
                      density=uncensored.density[["y"]])
data5_1$density.adjusted <- (max(na.omit(ROC.td[["AUC"]]))-min(na.omit(ROC.td[["AUC"]])))*
  (uncensored.density[["y"]]-min(uncensored.density[["y"]]))/
  (max(uncensored.density[["y"]])-min(uncensored.density[["y"]]))+min(na.omit(ROC.td[["AUC"]]))

# censored
censored.density <- density(So[,1][So[,2]==0])
data5_2 <- data.frame(time=censored.density[["x"]], 
                      density=censored.density[["y"]])
data5_2$density.adjusted <- (max(na.omit(ROC.td[["AUC"]]))-min(na.omit(ROC.td[["AUC"]])))*
  (censored.density[["y"]]-min(censored.density[["y"]]))/
  (max(censored.density[["y"]])-min(censored.density[["y"]]))+min(na.omit(ROC.td[["AUC"]]))

res.AUC <- data.frame(time=ROC.td[["times"]],
                      AUC=ROC.td[["AUC"]])
res.AUC <- melt(data = res.AUC, 
                id.vars = "time", 
                measure.vars = c("AUC"))
res.AUC$status <- sort(So)[,2]


# Time dependent AUC plot with ggplot2

res.AUC$censor <- ifelse(So[,2]==1,
                         1,min(na.omit(ROC.td[["AUC"]])))
ggplot() +
  geom_line(data=res.AUC, aes(x=time, y=value,
                              group=variable, color=variable)) +
  geom_point(data=res.AUC, aes(x=time, y=censor, color=as.factor(So[,2]))) +
  geom_line(data=data5_1, aes(x=time, y=density.adjusted), 
            size=1, colour='black') +
  geom_line(data=data5_2, aes(x=time, y=density.adjusted), 
            size=1, colour='orange') +
  scale_x_continuous(limits = c(0, max(So[,1]))) +
  geom_vline(xintercept = min(So[,1][So[,2]==1]),
             linetype="dotted", color = "blue", size=1) +
  ggtitle("Time dependent AUC with risk score (alpha=0.01, lambda=0.1, k=5, time.sort : 2)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))



#----------------#







library(survival)
So <- with(veteran, Surv(time,status==1))
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior), data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
risk.score_cox <- predict(object=coxph(So ~ X), newdata=veteran, type="risk") # risk score
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5, ncol(X)))
time.sort <- sort(So[,1]) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation



# Assume parametric basis distribution is Exponential distribution... 
lambda1 <- 0.003 # parameter for generate Y  
lambda2 <- 0.2/lambda1 # parameter for generate C 
beta0 <- 1:5 # true coefficients
Y <- c() # real survival time
C <- c() # censoring time
t <- c() # observed survival time
censor <- c() # censoring indicator
X <- matrix(, nrow=100, ncol=length(beta0)) # explanatory variables
U1 <- c() # random variable from U(0,1) for generate Y
U2 <- c() # random variable from U(0,1) for generate C 

for (i in 1:100) {
  U1[i] <- runif(1)
  U2[i] <- runif(1)
  X[i,] <- rnorm(5,mean=0,sd=1)
  Y[i] <- -(1/lambda1)*log(U1[i])*exp(-X[i,]%*%beta0)
  C[i] <- (lambda2)*(1-U2[i])
  t[i] <- min(Y[i],C[i])
  censor[i] <- ifelse(Y[i]<=C[i],1,0)
}

mean(1-censor) # censoring rate

dt5 <- data.frame(X=X, Y=Y, C=C, time=t, censor=censor)
head(dt5,10)

library(survival)
So <- with(dt5, Surv(time,censor==1))
X <- model.matrix(~X.1+X.2+X.3+X.4+X.5, data=dt5)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
risk.score_cox <- predict(object=coxph(So ~ X), newdata=dt5, type="risk") # risk score
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5, ncol(X)))
time.sort <- sort(So[,1]) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation



# alpha=0.01, lambda=0, k=5, number of sorting time=100 (case01)
time.sort_01 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res01 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0, k=5, time.sort_01)
risk.score_01 <- exp(X%*%res01$beta.new)
res01_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0, time.sort_01)
risk.score_01_1 <- exp(X%*%res01_1$beta.new)

# alpha=0.01, lambda=1, k=5, number of sorting time=100 (case02)
time.sort_02 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res02 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=1, k=5, time.sort_02)
risk.score_02 <- exp(X%*%res02$beta.new)
res02_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=1, time.sort_02)
risk.score_02_1 <- exp(X%*%res02_1$beta.new)

# alpha=0.01, lambda=10, k=5, number of sorting time=100 (case03)
time.sort_03 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res03 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=10, k=5, time.sort_03)
risk.score_03 <- exp(X%*%res03$beta.new)
res03_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=10, time.sort_03)
risk.score_03_1 <- exp(X%*%res03_1$beta.new)

# alpha=0.01, lambda=20, k=5, number of sorting time=100 (case04)
time.sort_04 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res04 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=20, k=5, time.sort_04)
risk.score_04 <- exp(X%*%res04$beta.new)
res04_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=20, time.sort_04)
risk.score_04_1 <- exp(X%*%res04_1$beta.new)



## Time-dependent AUC


# Loading packages

library(survival) 
library(timeROC)
library(timereg)
library(ggplot2)
library(reshape)
library(plyr)
library(tidyverse)
library(coxed)
library(gridExtra)
library(mgcv)


# Time dependent AUC (Cumulative/Dynamic time-dependent ROC curve)

# mini-batch gradient descent

ROC.00 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_cox,
                    cause=1,weighting="marginal",times=So[,1])
case00 <- ROC.00[["AUC"]]

ROC.td01 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_01,
                    cause=1,weighting="marginal",times=So[,1])
case01 <- ROC.td01[["AUC"]]

ROC.td02 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_02,
                    cause=1,weighting="marginal",times=So[,1])
case02 <- ROC.td02[["AUC"]]

ROC.td03 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_03,
                    cause=1,weighting="marginal",times=So[,1])
case03 <- ROC.td03[["AUC"]]

ROC.td04 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_04,
                    cause=1,weighting="marginal",times=So[,1])
case04 <- ROC.td04[["AUC"]]


dt <- data.frame(time=ROC.00[["times"]], case00=case00,
                 case01=case01, case02=case02, case03=case03,
                 case04=case04)
dt1 <- data.frame(time=rep(ROC.00[["times"]],5),
                  case=c(rep("case00",100), rep("case01",100), rep("case02",100),
                         rep("case03",100), rep("case04",100)),
                  AUC=c(case00, case01, case02, case03, case04))

dt2 <- na.omit(dt1)

aggregate(dt2$AUC, list(dt2$case), FUN=mean, na.action = na.omit)

ggplot(dt1, aes(x=time, y=AUC, group=case, color=case)) +
  geom_line() +
  ggtitle("Time dependent AUC with risk score (mini-batch gradient descent)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))

# gradient descent

ROC.td00_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_cox,
                  cause=1,weighting="marginal",times=So[,1])
case00_1 <- ROC.td00_1[["AUC"]]

ROC.td01_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_01_1,
                    cause=1,weighting="marginal",times=So[,1])
case01_1 <- ROC.td01_1[["AUC"]]

ROC.td02_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_02_1,
                    cause=1,weighting="marginal",times=So[,1])
case02_1 <- ROC.td02_1[["AUC"]]

ROC.td03_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_03_1,
                    cause=1,weighting="marginal",times=So[,1])
case03_1 <- ROC.td03_1[["AUC"]]

ROC.td04_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_04_1,
                      cause=1,weighting="marginal",times=So[,1])
case04_1 <- ROC.td04_1[["AUC"]]


dt_1 <- data.frame(time=ROC.td00_1[["times"]], case00_1=case00_1,
                 case01=case01_1, case02=case02_1, case03=case03_1,
                 case04=case04_1)
dt1_1 <- data.frame(time=rep(ROC.td00_1[["times"]],5),
                  case=c(rep("case00_1",100), rep("case01_1",100), rep("case02_1",100),
                         rep("case03_1",100), rep("case04_1",100)),
                  AUC=c(case00_1, case01_1, case02_1, case03_1, case04_1))

dt2_1 <- na.omit(dt1_1)

aggregate(dt2_1$AUC, list(dt2_1$case), FUN=mean, na.action = na.omit)

ggplot(dt1_1, aes(x=time, y=AUC, group=case, color=case)) +
  geom_line() +
  ggtitle("Time dependent AUC with risk score (gradient descent)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))


# ROC curve

TP1 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=11"]
FP1 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=11"]
TP2 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=25"]
FP2 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=25"]
TP3 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=36"]
FP3 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=36"]
TP4 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=411"]
FP4 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=411"]
TP5 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=587"]
FP5 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=587"]
TP <- c(TP1,TP2,TP3,TP4,TP5)
FP <- c(FP1,FP2,FP3,FP4,FP5)

time.survival <- c(rep("t=11",138),rep("t=25",552),rep("t=36",138),
                   rep("t=411",138),rep("t=587",138))
dt.ROC <- data.frame(times=time.survival,TP=TP,FP=FP)

ggplot(data=dt.ROC)+
  geom_line(mapping=aes(x=FP,
                        y=TP,
                        group=times,
                        color=times),
            size=1) +
  ggtitle("Time dependent ROC curve (veteran data)") +
  xlab("False Positive Rate") + ylab("True Positive Rate") +
  theme(plot.title = element_text(hjust = 0.5))


# Make AUC data

# uncensored
uncensored.density <- density(So[,1][So[,2]==1])
data5_1 <- data.frame(time=uncensored.density[["x"]], 
                      density=uncensored.density[["y"]])
data5_1$density.adjusted <- (max(na.omit(ROC.td[["AUC"]]))-min(na.omit(ROC.td[["AUC"]])))*
  (uncensored.density[["y"]]-min(uncensored.density[["y"]]))/
  (max(uncensored.density[["y"]])-min(uncensored.density[["y"]]))+min(na.omit(ROC.td[["AUC"]]))

# censored
censored.density <- density(So[,1][So[,2]==0])
data5_2 <- data.frame(time=censored.density[["x"]], 
                      density=censored.density[["y"]])
data5_2$density.adjusted <- (max(na.omit(ROC.td[["AUC"]]))-min(na.omit(ROC.td[["AUC"]])))*
  (censored.density[["y"]]-min(censored.density[["y"]]))/
  (max(censored.density[["y"]])-min(censored.density[["y"]]))+min(na.omit(ROC.td[["AUC"]]))

res.AUC <- data.frame(time=ROC.td[["times"]],
                      AUC=ROC.td[["AUC"]])
res.AUC <- melt(data = res.AUC, 
                id.vars = "time", 
                measure.vars = c("AUC"))
res.AUC$status <- sort(So)[,2]


# Time dependent AUC plot with ggplot2

res.AUC$censor <- ifelse(So[,2]==1,
                      1,min(na.omit(ROC.td[["AUC"]])))
ggplot() +
  geom_line(data=res.AUC, aes(x=time, y=value,
                           group=variable, color=variable)) +
  geom_point(data=res.AUC, aes(x=time, y=censor, color=as.factor(So[,2]))) +
  geom_line(data=data5_1, aes(x=time, y=density.adjusted), 
            size=1, colour='black') +
  geom_line(data=data5_2, aes(x=time, y=density.adjusted), 
            size=1, colour='orange') +
  scale_x_continuous(limits = c(0, max(So[,1]))) +
  geom_vline(xintercept = min(So[,1][So[,2]==1]),
             linetype="dotted", color = "blue", size=1) +
  ggtitle("Time dependent AUC with risk score (alpha=0.01, lambda=0.1, k=5, time.sort : 2)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))
