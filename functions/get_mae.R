# Function computing the mean adjusted estimator for a given data set. It outputs the MAEs using theta_1 estimated from both periods, only period 1, only period 2, as well as using the true value of theta_1.

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
  
  theta_1_true <- theta[1]
  theta_1_per12 <- mean(Data[Data$treatment==1,]$response) - mean(Data[Data$treatment==0,]$response)
  theta_1_per1 <- y_11 - y_01
  theta_1_per2 <- y_12 - y_02
  
  rho <- (1/n_02)/(1/n_01 + 1/n_02 + 1/n_11 + 1/n_12)
  
  theta2_tilde_unadj <- (y_22 - ((1-rho)*y_02 + rho*(y_01+y_12-y_11)))
  
  theta2_tilde_MAE <- function(theta_1){
    
    sd_per1 <- sqrt(sigma^2/n_11 + sigma^2/n_01)
    
    bias <- rho*pnorm(c_1 - theta_1)*sd_per1*(dnorm(c_1 - (theta_1/sd_per1))/(pnorm(-c_1 + (theta_1/sd_per1))*pnorm(c_1 - (theta_1/sd_per1))))
    
    return(list(est = theta2_tilde_unadj - bias, 
                bias = bias))
  }
  
  theta2_tilde_adj_theta1_true <- theta2_tilde_MAE(theta_1_true)$est
  
  theta2_tilde_adj_theta1_per12 <- theta2_tilde_MAE(theta_1_per12)$est
  
  theta2_tilde_adj_theta1_per1 <- theta2_tilde_MAE(theta_1_per1)$est
  
  theta2_tilde_adj_theta1_per2 <- theta2_tilde_MAE(theta_1_per2)$est
  
  bias_est_theta1_true <- theta2_tilde_MAE(theta_1_true)$bias
  
  bias_est_theta1_per12 <- theta2_tilde_MAE(theta_1_per12)$bias
  
  bias_est_theta1_per1 <- theta2_tilde_MAE(theta_1_per1)$bias
  
  bias_est_theta1_per2 <- theta2_tilde_MAE(theta_1_per2)$bias
  
  return(list(theta2_tilde_unadj = theta2_tilde_unadj,
              theta2_tilde_adj_theta1_true = theta2_tilde_adj_theta1_true,
              theta2_tilde_adj_theta1_per12 = theta2_tilde_adj_theta1_per12,
              theta2_tilde_adj_theta1_per1 = theta2_tilde_adj_theta1_per1,
              theta2_tilde_adj_theta1_per2 = theta2_tilde_adj_theta1_per2,
              bias_est_theta1_true = bias_est_theta1_true,
              bias_est_theta1_per12 = bias_est_theta1_per12,
              bias_est_theta1_per1 = bias_est_theta1_per1,
              bias_est_theta1_per2 = bias_est_theta1_per2))
}


