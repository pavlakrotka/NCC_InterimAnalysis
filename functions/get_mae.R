# Function for computing the mean adjusted estimator for a given data set. It outputs the MAEs using theta_1 estimated from both periods, only period 1, only period 2, as well as using the true value of theta_1.

# Arguments:

# Data: Data frame. Contains the trial data from the evaluated platform trial. Must contain columns `response`, `treatment`, and `period_2`.
# theta: Double vector of length 2 with treatment effects for arms 1 and 2 in terms of difference of means compared to control.
# sigma: Double. Standard deviation of the responses.
# futility_bound: Double. Futility bound alpha_1 used in the interim analysis for arm 1.

get_mae <- function(Data, theta, sigma, futility_bound){
  
  n_22 <- nrow(Data[Data$treatment==2 & Data$period_2==2,])
  n_02 <- nrow(Data[Data$treatment==0 & Data$period_2==2,])
  n_01 <- nrow(Data[Data$treatment==0 & Data$period_2==1,])
  n_12 <- nrow(Data[Data$treatment==1 & Data$period_2==2,])
  n_11 <- nrow(Data[Data$treatment==1 & Data$period_2==1,])
  
  c_1 <- qnorm(1-futility_bound)
  
  y_22 <- mean(Data[Data$treatment==2 & Data$period_2==2,]$response)
  y_02 <- mean(Data[Data$treatment==0 & Data$period_2==2,]$response)
  y_01 <- mean(Data[Data$treatment==0 & Data$period_2==1,]$response)
  y_12 <- mean(Data[Data$treatment==1 & Data$period_2==2,]$response)
  y_11 <- mean(Data[Data$treatment==1 & Data$period_2==1,]$response)
  
  get_cumvue <- function(Data){
    l <- c_1
    u <- Inf
    
    z <- (mean(Data[Data$treatment==1,]$response) - mean(Data[Data$treatment==0,]$response))/sqrt(sigma^2/(n_01+n_02) + sigma^2/(n_11+n_12))
    
    I1 <- 1/(sigma^2*(1/n_11 + 1/n_01))
    
    I2 <- 1/(sigma^2*(1/(n_11 + n_12) + 1/(n_01 + n_02)))
    
    mle <- z/sqrt(I2)
    
    corr_term <- (I2 - I1) / (I2 * sqrt(I1)) *
      ( dnorm(u, mean = z * sqrt(I1 / I2), sd = sqrt((I2 - I1) / I2)) -
          dnorm(l, mean = z * sqrt(I1 / I2), sd = sqrt((I2 - I1) / I2)) ) /
      ( pnorm(u, mean = z * sqrt(I1 / I2), sd = sqrt((I2 - I1) / I2)) -
          pnorm(l, mean = z * sqrt(I1 / I2), sd = sqrt((I2 - I1) / I2)) )
    
    UMVUE <- mle - corr_term
    
    CUMVUE <- (z*sqrt(I2) - I1*UMVUE)/(I2-I1)
    CUMVUE
  }
  
  theta_1_true <- theta[1]
  theta_1_per12 <- mean(Data[Data$treatment==1,]$response) - mean(Data[Data$treatment==0,]$response)
  theta_1_per1 <- y_11 - y_01
  theta_1_per2 <- y_12 - y_02
  theta_1_rao <- get_cumvue(Data)
  
  rho <- (1/n_02)/(1/n_01 + 1/n_02 + 1/n_11 + 1/n_12)
  
  theta2_tilde_unadj <- (y_22 - ((1-rho)*y_02 + rho*(y_01+y_12-y_11)))
  
  theta2_tilde_MAE <- function(theta_1){
    
    # sd_per1 <- sqrt(sigma^2/n_11 + sigma^2/n_01)
    
    sd_per1 <- sigma*sqrt(1/n_11 + 1/n_01)
    
    # bias <- rho*pnorm(c_1 - theta_1)*sd_per1*(dnorm(c_1 - (theta_1/sd_per1))/(pnorm(-c_1 + (theta_1/sd_per1))*pnorm(c_1 - (theta_1/sd_per1))))
    
    # bias <- (rho*sd_per1*dnorm(c_1 - (theta_1/sd_per1)))/(1-pnorm(c_1 - (theta_1/sd_per1)))
    
    # equivalent to the formula in paper, but using x=exp(log(x)) to avoid computational problems:
    bias <- exp(log(rho*sd_per1)+dnorm((c_1 - (theta_1/sd_per1)), log = T) - (pnorm((-c_1 + (theta_1/sd_per1)), log.p = T)))
    
    return(list(est = theta2_tilde_unadj - bias, 
                bias = bias))
  }
  
  theta2_tilde_adj_theta1_true <- theta2_tilde_MAE(theta_1_true)$est
  
  theta2_tilde_adj_theta1_per12 <- theta2_tilde_MAE(theta_1_per12)$est
  
  theta2_tilde_adj_theta1_per1 <- theta2_tilde_MAE(theta_1_per1)$est
  
  theta2_tilde_adj_theta1_per2 <- theta2_tilde_MAE(theta_1_per2)$est
  
  theta2_tilde_adj_theta1_rao <- theta2_tilde_MAE(theta_1_rao)$est
  
  
  bias_est_theta1_true <- theta2_tilde_MAE(theta_1_true)$bias
  
  bias_est_theta1_per12 <- theta2_tilde_MAE(theta_1_per12)$bias
  
  bias_est_theta1_per1 <- theta2_tilde_MAE(theta_1_per1)$bias
  
  bias_est_theta1_per2 <- theta2_tilde_MAE(theta_1_per2)$bias
  
  bias_est_theta1_rao <- theta2_tilde_MAE(theta_1_rao)$bias
  
  return(list(theta2_tilde_unadj = theta2_tilde_unadj,
              theta2_tilde_adj_theta1_true = theta2_tilde_adj_theta1_true,
              theta2_tilde_adj_theta1_per12 = theta2_tilde_adj_theta1_per12,
              theta2_tilde_adj_theta1_per1 = theta2_tilde_adj_theta1_per1,
              theta2_tilde_adj_theta1_per2 = theta2_tilde_adj_theta1_per2,
              theta2_tilde_adj_theta1_rao = theta2_tilde_adj_theta1_rao,
              
              bias_est_theta1_true = bias_est_theta1_true,
              bias_est_theta1_per12 = bias_est_theta1_per12,
              bias_est_theta1_per1 = bias_est_theta1_per1,
              bias_est_theta1_per2 = bias_est_theta1_per2,
              bias_est_theta1_rao = bias_est_theta1_rao))
}


