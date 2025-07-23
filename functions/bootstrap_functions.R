# Function for performing stratified bootstrap with interim analysis to compute the variance of the mean adjusted estimator

# Arguments:

# Data: Data frame. Contains the trial data from the evaluated platform trial. Must contain columns `response`, `treatment`, and `period_2`.
# theta: Double vector of length 2 with treatment effects for arms 1 and 2 in terms of difference of means compared to control.
# sigma: Double. Standard deviation of the responses.
# futility_bound: Double. Futility bound alpha_1 used in the interim analysis for arm 1.
# B_boot: Integer. Number of bootstrap samples.

bootstrap_stratified_withIA <- function(Data, theta, sigma, futility_bound, B_boot = 1000){
  
  B_boot_aux <- B_boot
  mae_boot_unadj <- rep(NA, B_boot) 
  mae_boot_adj_theta1_true <- rep(NA, B_boot) 
  mae_boot_adj_theta1_per12 <- rep(NA, B_boot) 
  mae_boot_adj_theta1_per1 <- rep(NA, B_boot)
  mae_boot_adj_theta1_per2 <- rep(NA, B_boot) 
  
  boot_stop <- 0
  
  while (B_boot_aux>0) {
    
    Data_01 <- Data[Data$treatment==0 & Data$period_2==1,]
    Data_11 <- Data[Data$treatment==1 & Data$period_2==1,]
    
    n_01 <- nrow(Data_01)
    n_11 <- nrow(Data_11)
    
    Data_stage_1_boot <- rbind(Data_01[sample(1:nrow(Data_01), nrow(Data_01), replace = T),],
                               Data_11[sample(1:nrow(Data_11), nrow(Data_11), replace = T),])
    
    Z_1_boot <- (mean(Data_stage_1_boot[Data_stage_1_boot$treatment==1,]$response) - mean(Data_stage_1_boot[Data_stage_1_boot$treatment==0,]$response))/sqrt(sigma^2/n_11 + sigma^2/n_01)
    
    c_1_boot <- qnorm(1-futility_bound)
    
    if (Z_1_boot<c_1_boot) { # arm 1 stops for futility
      boot_stop <- boot_stop+1
      
    } else {
      
      Data_02 <- Data[Data$treatment==0 & Data$period_2==2,]
      Data_12 <- Data[Data$treatment==1 & Data$period_2==2,]
      Data_22 <- Data[Data$treatment==2 & Data$period_2==2,]
      
      Data_period_2_boot <- rbind(Data_02[sample(1:nrow(Data_02), nrow(Data_02), replace = T),],
                                  Data_12[sample(1:nrow(Data_12), nrow(Data_12), replace = T),],
                                  Data_22[sample(1:nrow(Data_22), nrow(Data_22), replace = T),])
      
      Data_stage_2_boot <- rbind(Data_stage_1_boot, Data_period_2_boot)
      
      MAEs_boot <- get_mae(Data_stage_2_boot, theta, sigma, futility_bound)
      
      mae_boot_unadj[B_boot - B_boot_aux +1] <- MAEs_boot$theta2_tilde_unadj
      mae_boot_adj_theta1_true[B_boot - B_boot_aux +1] <- MAEs_boot$theta2_tilde_adj_theta1_true
      mae_boot_adj_theta1_per12[B_boot - B_boot_aux +1] <- MAEs_boot$theta2_tilde_adj_theta1_per12
      mae_boot_adj_theta1_per1[B_boot - B_boot_aux +1] <- MAEs_boot$theta2_tilde_adj_theta1_per1
      mae_boot_adj_theta1_per2[B_boot - B_boot_aux +1] <- MAEs_boot$theta2_tilde_adj_theta1_per2
      
      B_boot_aux <- B_boot_aux-1 # for while loop
      
    }
    
  }
  
  var_boot_est_unadj <- 1/B_boot * sum((mae_boot_unadj-mean(mae_boot_unadj))^2)    
  var_boot_est_adj_theta1_true <- 1/B_boot * sum((mae_boot_adj_theta1_true-mean(mae_boot_adj_theta1_true))^2)
  var_boot_est_adj_theta1_per12 <- 1/B_boot * sum((mae_boot_adj_theta1_per12-mean(mae_boot_adj_theta1_per12))^2)
  var_boot_est_adj_theta1_per1 <- 1/B_boot * sum((mae_boot_adj_theta1_per1-mean(mae_boot_adj_theta1_per1))^2)
  var_boot_est_adj_theta1_per2 <- 1/B_boot * sum((mae_boot_adj_theta1_per2-mean(mae_boot_adj_theta1_per2))^2)
  
  return(list(var_boot_est_unadj = var_boot_est_unadj,
              var_boot_est_adj_theta1_true = var_boot_est_adj_theta1_true,
              var_boot_est_adj_theta1_per12 = var_boot_est_adj_theta1_per12,
              var_boot_est_adj_theta1_per1 = var_boot_est_adj_theta1_per1,
              var_boot_est_adj_theta1_per2 = var_boot_est_adj_theta1_per2,
              check = length(mae_boot_adj_theta1_per2),
              boot_stop = boot_stop))
}

