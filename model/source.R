standarized_errors_SLR <- function(y, x) {
  n <- length(y)
  linearModelVar <- lm(y ~ x)
  a <- coef(linearModelVar)[1] 
  b <- coef(linearModelVar)[2]
  se <- sigma(linearModelVar) 
  x_mean <- mean(x)
  x_var <- var(x)
  h <- vector("double", length = n)
  y_p <- vector("double", length = n)
  Error <- vector("double", length = n)
  Standard_E <- vector("double", length = n)
  for(i in 1:n) {
    h[i] <- (1/n) + (x[i] - x_mean)^2 / (x_var * (n - 1))
    y_p[i] <- a + b * x[i]
    Error[i] <- y[i] - y_p[i]
    Standard_E[i] <- Error[i] / (se * sqrt(1 - h[i]))
  }
  SEh <-cbind(Standard_E,h)
  return(SEh)    #return a matrix: standardized_residual & h
}


Run_Test <- function(x) {  
  size <- length(x)   
  Run_Count <- vector("double", length = size)   
  x_median <- median(x)   
  x_cat <- x > x_median   
  Run_Count[1] <- 1   
  for(j in 2:size) {      
    if (x_cat[j] == x_cat[j-1]) {          
      Run_Count[j] <- Run_Count[j-1]      
    } else {          
      Run_Count[j] <- Run_Count[j-1] + 1
    }   
  }
  n1 <- sum(x_cat == 'TRUE')   
  n2 <- sum(x_cat == 'FALSE')   
  Run_Mean <- 2*n1*n2/(n1+n2) + 1   
  Run_S <- sqrt(2*n1*n2*(2*n1*n2-n1-n2)/(n1+n2)^2/(n1+n2-1))   
  Run_Z = (Run_Count[size] - Run_Mean) / Run_S   
  if (Run_Z < 0) {       
    Run_Pvalue <- pnorm(Run_Count[size], mean=Run_Mean, sd=Run_S, lower.tail=TRUE)*2   
  } else {       
    Run_Pvalue <- pnorm(Run_Count[size], mean=Run_Mean, sd=Run_S, lower.tail=FALSE)*2 
  }    
  Run_Return_V <- vector("double", length = 3)   #small sample
  Run_Return_V[1] <- n1   
  Run_Return_V[2] <- n2   
  Run_Return_V[3] <- Run_Count[size]   
  if (n1 > 20 || n2 > 20 || (n1+n2) > 40)         #large sample
    return (Run_Pvalue)   
  return(Run_Return_V) 
}


standarized_errors_MR <- function(y, x) {
  n <- length(y)
  k <- ncol(x)
  x1 <- vector("double", length = n)
  for(j in 1:n) {
    x1[j] <- 1
  }
  Matrix_X <- cbind(x1, x)
  Matrix_Y <- cbind(y)
  Matrix_H <- Matrix_X%*%solve(t(Matrix_X)%*%Matrix_X)%*%t(Matrix_X)
  Matrix_YP <- Matrix_H%*%Matrix_Y
  y_p <- Matrix_YP[,1]
  h <- vector("double", length = n)
  Error <- vector("double", length = n)
  Standard_E <- vector("double", length = n)
  D <- vector("double", length = n)
  se <- sqrt(sum((y - y_p)^2)/(n - k - 1))
  for(i in 1:n) {
    h[i] <- Matrix_H[i,i]
    Error[i] <- y[i] - y_p[i]
    Standard_E[i] <- Error[i] / (se * sqrt(1 - h[i]))
    D[i] <- (y[i] - y_p[i])^2 * h[i] / ((k -1) * se^2 * (1 - h[i])^2)
  }
  SEhD <- cbind(Standard_E, h, D)
  return(SEhD)
}


Durbin_Watson_Test <- function(x) {
  x_square_sum <- sum(x*x)
  size <- length(x)
  x_d <- vector("double", length = n)
  x_d[1] = 0
  for(j in 2:n) {
    x_d[j] <- x[j] - x[j - 1]
  }
  d <- sum(x_d*x_d) / x_square_sum
  return(d)
} 



Center_Moving_Average <- function(x, q, s) {
  k <- length(x)
  y <- vector("double", length = k)
  SE <- vector("double", length = k)
  SI <- vector("double", length = s)
  if(s %% 2 == 0){
    lb = s / 2
    lb2 = lb + 1
    ub = k - (s / 2)
    cmv_y <- vector("double", length = k)
    for(i in lb:ub) {
      sum_x = 0
      lb_cmv = i - lb + 1
      ub_cmv = i + lb
      for(j in lb_cmv:ub_cmv){
        sum_x <- sum_x + x[j]
      }
      cmv_y[i] <- sum_x / s
    }
    for(i in lb2:ub) {
      y[i] <- (cmv_y[i-1] + cmv_y[i]) / 2
      SE[i] <- x[i] / y[i]
      SI[q[i]] <- SI[q[i]] + SE[i]
    }
    SI <- SI / ((k / s) - 1)
    sum_SI <- sum(SI)
    SI <- SI * s / sum_SI
  } else {
    lb = ceiling(s / 2) 
    ub = k - floor(s / 2)
    for(i in lb:ub) {
      sum_x = 0
      lb_cmv = i - lb + 1
      ub_cmv = i + lb - 1
      for(j in lb_cmv:ub_cmv){
        sum_x <- sum_x + x[j]
      }
      y[i] <- sum_x / s
      SE[i] <- x[i] / y[i]
      SI[q[i]] <- SI[q[i]] + SE[i]
    }
    SI <- SI / ((k / s) - 1)
    sum_SI <- sum(SI)        #normalize
    SI <- SI * s / sum_SI
  }
  return(SI)
}

Linear_Regression_Seasonal_Index <- function(x, q, s) {
  k <- length(x)
  SE <- vector("double", length = k)
  SI <- vector("double", length = s)
  t <- 0:(k - 1)
  linearModelVar <- lm(x ~ t)      #先做一條回歸線，用來算seasonal index
  b0 <- coef(linearModelVar)[1]
  b1 <- coef(linearModelVar)[2]
  y <- b0 + b1 * t
  for(i in 1:k) {
    SE[i] <- x[i] / y[i]
    SI[q[i]] <- SI[q[i]] + SE[i]
  }
  SI <- SI / ((k / s) - 1)
  sum_SI <- sum(SI)
  SI <- SI * s / sum_SI
  return(SI)
}


###
Forecast_by_SI <- function(x, q, s, SI, t) {   # t: 0~n-1
  #Deseasonalizing
  k <- length(x)
  Des_x <- vector("double", length = k)
  for(i in 1:k) {
    Des_x[i] <- x[i] / SI[q[i]]
  } 
  linearModelVar <- lm(Des_x ~ t)
  b0 <- coef(linearModelVar)[1]
  b1 <- coef(linearModelVar)[2]
  
  
  cat("### Linear Model ###")
  print(summary(linearModelVar))
  
  #residual analysis 
  Standard_E <- standarized_errors_SLR(Des_x, t)[ ,1]
  h <- standarized_errors_SLR(Des_x, t)[ ,2]
  print(shapiro.test(Standard_E) )                # normality
  y_p <- vector("double", length = k)      # Homoscedasticity & Heteroscedasticity 
  for(i in 1:k) {
    y_p[i] <- b0 + b1 * t[i]
  }
  plot(x = y_p, y = Standard_E, xlab="Predicted Des_x", ylab="Standardized Error", main="Predicted Des_x vs Error") 
  cat("run test:")
  print(Run_Test(Standard_E))        #errors are independent

  cat("The Outliers")
  Outliers <- abs(Standard_E) > 2
  if(sum(Outliers)!=0)
      print(which(Outliers))
  else
    print("no outliers")
  
  cat("influential observations:")
  Inf_Obs <- h > 3 * (12+1) / n
  if(sum(Inf_Obs)!=0)
     print(which(Inf_Obs))
  else
    print("no influential observations")
  
  
  #assesment
  cat("mean_y:",mean(Des_x),"\n\n\n")
  
  #forecast
  new_k <- k + s
  new_lb <- k + 1
  Q_t <- vector("double", length = new_k)
  SI_Q_t <- vector("double", length = new_k)
  for(i in 1:k) {
    Q_t[i] <- q[i]
  }
  for(j in new_lb:new_k) {
    Q_t[j] <- j - k
  }
  for(i in 1:new_k) {
    SI_Q_t[i] <- SI[Q_t[i]]
  }
  F_t <- 0:(new_k - 1)
  des_y <- b0 + b1 * F_t
  y <- des_y * SI_Q_t
  
  
  cat("\ndeseasonalized data:\n")
  print(Des_x)
  cat("\n\n\n")
  
  
  #plot the deseasonalized data
  plot(t,Des_x, col = "green" , main = "Visitors against Period", xlab = "Period", ylab = "Vistors", xlim=range(t, t+13), ylim=range(max(y)+Des_x/8,0))
  lines(t,Des_x, type="o", col = "green")  
  abline(lm(Des_x ~ t))

  
  return (y)
}



#Error Metrics
Mean_Absolute_Deviation <- function(x, y, s) {
  k <- length(x)
  Sum_E <- 0
  for(i in s:k) {
    Sum_E <- Sum_E + abs(x[i] - y[i])
  }
  MAD <- Sum_E / (k - s + 1)
  return(MAD)
}

Mean_Square_Error <- function(x, y, s) {
  k <- length(x)
  Sum_E <- 0
  for(i in s:k) {
    Sum_E <- Sum_E + (x[i] - y[i])^2
  }
  MSE <- Sum_E / (k - s + 1)
  return(MSE)
}

Mean_Absolute_Percentage_Error <- function(x, y, s) {
  k <- length(x)
  Sum_E <- 0
  for(i in s:k) {
    if(x[i] != 0) {
      Sum_E <- Sum_E + abs(x[i] - y[i])/x[i]
    } else {
      Sum_E <- Sum_E + abs(x[i] - y[i])/mean(x)
    }  
  }
  MAPE <- Sum_E / (k - s + 1) * 100
  return(MAPE)
}