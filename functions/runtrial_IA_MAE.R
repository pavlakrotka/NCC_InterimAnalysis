# Function for running a platform trial with 2 treatment arms and an interim analysis of arm 1 at the time of adding the arm 2, where arm 2 is evaluated using the mean adjusted estimator
# For the data generation, user-defined sample sizes per arm per period and block randomization are assumed
# The sample sizes has to be chosen such that n_11/n_01 or n_01/n_11, as well as n_12/n_02 or n_02/n_12, and n_22/n_02 or n_02/n_22 are natural numbers in order to enable the block randomization.
# Linear and stepwise time trends of strength lambda for all arms can be present in the trial
# Interim analysis is performed for arm 1 with the possibility to stop for futility

# Arguments:

# n_01: Integer. Sample size in the control arm in period 1.
# n_11: Integer. Sample size in arm 1 in period 1.
# n_02: Integer. Sample size in the control arm in period 2.
# n_12: Integer. Sample size in arm 1 in period 2.
# n_22: Integer. Sample size in arm 2 in period 2.
# mu0: Double. Response in the control arm.
# theta: Double vector of length 2 with treatment effects for arms 1 and 2 in terms of difference of means compared to control.
# sigma: Double. Standard deviation of the responses.
# futility_bound: Double. Futility bound alpha_1 used in the interim analysis for arm 1
# alpha: Double. Significance level (one-sided) for the interim analysis of arm 1 and final analysis of arm 2.
# lambda: Integer vector of length 3 containing the strength of the time trend in each arm (control, arm 1, arm 2).
# trend_pattern: String indicating the time trend pattern ("linear" or "stepwise")
# period_blocks: Integer. Number to define the size of the blocks for the block randomization. The block size in each period equals `period_blocks` times the smallest block size that would be possible for a given allocation ratio. For example, if the allocation ratio is 3:1, the block size in this period is `period_blocks`*(3+1)=8.
# B_boot: Integer. Number of bootstrap samples. Default=1000.

library(rlang)

runtrial_IA_MAE <- function(n_01 = 100,
                            n_11 = 100,
                            n_02 = 150,
                            n_12 = 150,
                            n_22 = 150,
                            
                            mu0 = 0,
                            theta = c(0, 0),
                            sigma = 1,
                            
                            futility_bound = 0.3,
                            alpha = 0.025,
                            
                            lambda = rep(0, 3),
                            trend_pattern = "stepwise",
                            period_blocks = 2,
                            B_boot = 1000){  
  
  
  num_arms <- 2
  
  SS_matrix <- matrix(c(n_01, n_11, NA,
                        n_02, n_12, n_22), ncol = 2, byrow = F)
  
  if(sum(SS_matrix[1,]>SS_matrix[2,])==0){ # if sample size in arm 1 larger than control
    
    alloc_ratios <- matrix(NA, nrow = 3, ncol = 2)
    
    for (i in 1:nrow(SS_matrix)) {
      alloc_ratios[i,] <- SS_matrix[i,]/SS_matrix[1,]
    }
    
    alloc_ratios <- (ifelse(!is.na(alloc_ratios), alloc_ratios, 0))
  }
  
  
  if(sum(SS_matrix[1,]>SS_matrix[2,])==2){ # if sample size in control larger than arm 1
    
    alloc_ratios <- matrix(NA, nrow = 3, ncol = 2)
    
    for (i in 1:nrow(SS_matrix)) {
      alloc_ratios[i,] <- SS_matrix[i,]/SS_matrix[2,]
    }
    
    alloc_ratios <- (ifelse(!is.na(alloc_ratios), alloc_ratios, 0))
  }
  
  is_natural_matrix <- function(mat) {
    is.numeric(mat) && all(mat == floor(mat)) && all(mat >= 0)
  }
  
  if(is_natural_matrix(alloc_ratios)==FALSE){
    stop("Sample sizes has to be chosen such that allocation ratios are natural numbers!")
  }
  
  
  num_periods <- ncol(alloc_ratios) # total number of periods
  
  N_period <- colSums(SS_matrix, na.rm = T) # sample sizes per period
  N_arm <- rowSums(SS_matrix, na.rm = T) # sample sizes per arm
  n_total <- sum(SS_matrix, na.rm = T) # total sample size
  active_arms <- colSums(apply(SS_matrix, 2, is.na)==0) # active arms per period
  
  block_sizes <- period_blocks*colSums(alloc_ratios) # block sizes per period
  
  linear_trend <- function(j, lambda){
    lambda*(j-1)/(n_total-1)
  }
  
  sw_trend <- function(cj, lambda){
    as.numeric(lambda)*(cj-1)
  }
  
  # Generate data
  
  t <- c()
  
  for (i in 1:num_periods){
    
    m_i <- t(replicate(trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i]),
                       sample(rep(rep(c(0:(num_arms)), alloc_ratios[,i]), block_sizes[i]/length(rep(c(0:(num_arms)), alloc_ratios[,i]))))))
    
    rest <- sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])
    
    t_i <- c(t(m_i), sample(rep(c(0:(num_arms)), alloc_ratios[,i]*ceiling(rest/active_arms[i])),
                            size = sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])))
    
    t <- c(t, t_i)
  }
  
  j <- c(1:n_total)
  
  
  for (i in 0:num_arms) {
    assign(paste0("j", i), which(t==i)) # j0, j1, j2 ... position in time (order) of allocated patients in every arm
  }
  
  cj <- rep(1:num_periods, N_period)
  
  if(trend_pattern=="linear"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                  lambda = lambda[i+1]))
      
      assign(paste0("all_trend", i), linear_trend(j=j,
                                                  lambda = lambda[i+1]))
    }
  }
  
  if(trend_pattern=="stepwise"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), sw_trend(cj=cj[eval(sym(paste0("j", i)))],
                                              lambda = lambda[i+1]))
      
      assign(paste0("all_trend", i), sw_trend(cj=cj,
                                              lambda = lambda[i+1]))
    }
  }
  
  
  # Simulation of continuous endpoint
  
  means <- c()
  means[j0] <- ind_trend0
  
  for (i in 1:num_arms) {
    means[eval(sym(paste0("j", i)))] <- eval(sym(paste0("ind_trend", i))) + theta[i]
  }
  
  X <- rnorm(n=n_total, mean=mu0+means, sd=sigma)
  
  Data <- data.frame(j = c(1:n_total),
                     response = X,
                     treatment = t,
                     period_2 = rep(1:num_periods, N_period),
                     means = mu0+means)
  
  for (i in 0:num_arms) {
    Data[ ,paste0("lambda", i)] <- lambda[i+1]
  }
  
  
  Data_stage_1 <- Data[Data$period_2==1,]
  
  Z_1 <- (mean(Data_stage_1[Data_stage_1$treatment==1,]$response) - mean(Data_stage_1[Data_stage_1$treatment==0,]$response))/sqrt(sigma^2/n_01 + sigma^2/n_11)
  
  c_1 <- qnorm(1-futility_bound)
  
  if (Z_1<c_1) { # arm 1 stops for futility
    Data_stage_2 <- Data[Data$period_2==2,]
    
    theta2_tilde_unadj <- theta2_tilde_adj_theta1_true <- theta2_tilde_adj_theta1_per12 <- theta2_tilde_adj_theta1_per1 <- theta2_tilde_adj_theta1_per2 <- mean(Data_stage_2[Data_stage_2$treatment==2,]$response) - mean(Data_stage_2[Data_stage_2$treatment==0,]$response)
    
    bias_arm2_unadj <- bias_arm2_adj_theta1_true <- bias_arm2_adj_theta1_per12 <- bias_arm2_adj_theta1_per1 <- bias_arm2_adj_theta1_per2 <- theta2_tilde_unadj - theta[2]
    
    theta2_tilde_unadj_stop <- mean(Data_stage_2[Data_stage_2$treatment==2,]$response) - mean(Data_stage_2[Data_stage_2$treatment==0,]$response)
    
    bias_arm2_unadj_stop <- theta2_tilde_unadj_stop - theta[2]
    
    # Z-test for arm 2
    Z_stat <- theta2_tilde_unadj / sqrt(sigma^2/n_02 + sigma^2/n_22)
    
    reject_h02_unadj <- reject_h02_adj_theta1_true <- reject_h02_adj_theta1_per12 <- reject_h02_adj_theta1_per1 <- reject_h02_adj_theta1_per2 <- Z_stat > qnorm(1-alpha)
    
    reject_h02_unadj_stop <- Z_stat > qnorm(1-alpha)
    
    boot_stop <- NA
    
    var_boot_est_unadj <- NA
    var_boot_est_adj_theta1_true <- NA
    var_boot_est_adj_theta1_per12 <- NA
    var_boot_est_adj_theta1_per1 <- NA
    var_boot_est_adj_theta1_per2 <- NA
    
    bias_est_theta1_true <- NA
    bias_est_theta1_per12 <- NA
    bias_est_theta1_per1 <- NA
    bias_est_theta1_per2 <- NA
    
    theta2_tilde_unadj_cond <- NA
    theta2_tilde_adj_theta1_true_cond <- NA
    theta2_tilde_adj_theta1_per12_cond <- NA
    theta2_tilde_adj_theta1_per1_cond <- NA
    theta2_tilde_adj_theta1_per2_cond <- NA
    
    bias_arm2_unadj_cond <- NA
    bias_arm2_adj_theta1_true_cond <- NA
    bias_arm2_adj_theta1_per12_cond <- NA
    bias_arm2_adj_theta1_per1_cond <- NA
    bias_arm2_adj_theta1_per2_cond <- NA
    
    reject_h02_unadj_cond <- NA
    reject_h02_adj_theta1_true_cond <- NA
    reject_h02_adj_theta1_per12_cond <- NA
    reject_h02_adj_theta1_per1_cond <- NA
    reject_h02_adj_theta1_per2_cond <- NA
    
  } else { # arm 1 continues
    
    theta2_tilde_unadj_stop <- NA
    bias_arm2_unadj_stop <- NA
    reject_h02_unadj_stop <- NA
    
    Data_stage_2 <- Data
    
    MAEs <- get_mae(Data_stage_2, theta, sigma, futility_bound)
    
    theta2_tilde_unadj <- MAEs$theta2_tilde_unadj
    theta2_tilde_adj_theta1_true <- MAEs$theta2_tilde_adj_theta1_true
    theta2_tilde_adj_theta1_per12 <- MAEs$theta2_tilde_adj_theta1_per12
    theta2_tilde_adj_theta1_per1 <- MAEs$theta2_tilde_adj_theta1_per1
    theta2_tilde_adj_theta1_per2 <- MAEs$theta2_tilde_adj_theta1_per2
    
    bias_est_theta1_true <- MAEs$bias_est_theta1_true
    bias_est_theta1_per12 <- MAEs$bias_est_theta1_per12
    bias_est_theta1_per1 <- MAEs$bias_est_theta1_per1
    bias_est_theta1_per2 <- MAEs$bias_est_theta1_per2
    
    bias_arm2_unadj <- theta2_tilde_unadj - theta[2] 
    bias_arm2_adj_theta1_true <- theta2_tilde_adj_theta1_true - theta[2] 
    bias_arm2_adj_theta1_per12 <- theta2_tilde_adj_theta1_per12 - theta[2] 
    bias_arm2_adj_theta1_per1 <- theta2_tilde_adj_theta1_per1 - theta[2] 
    bias_arm2_adj_theta1_per2 <- theta2_tilde_adj_theta1_per2 - theta[2] 
    
    theta2_tilde_unadj_cond <- MAEs$theta2_tilde_unadj
    theta2_tilde_adj_theta1_true_cond <- MAEs$theta2_tilde_adj_theta1_true
    theta2_tilde_adj_theta1_per12_cond <- MAEs$theta2_tilde_adj_theta1_per12
    theta2_tilde_adj_theta1_per1_cond <- MAEs$theta2_tilde_adj_theta1_per1
    theta2_tilde_adj_theta1_per2_cond <- MAEs$theta2_tilde_adj_theta1_per2
    
    bias_arm2_unadj_cond <- theta2_tilde_unadj - theta[2] 
    bias_arm2_adj_theta1_true_cond <- theta2_tilde_adj_theta1_true - theta[2] 
    bias_arm2_adj_theta1_per12_cond <- theta2_tilde_adj_theta1_per12 - theta[2] 
    bias_arm2_adj_theta1_per1_cond <- theta2_tilde_adj_theta1_per1 - theta[2] 
    bias_arm2_adj_theta1_per2_cond <- theta2_tilde_adj_theta1_per2 - theta[2] 
    
    
    # Bootstrap variance of MAE
    
    res_bootstrap <- bootstrap_stratified_withIA(Data_stage_2, theta, sigma, futility_bound, B_boot)
    
    var_boot_est_unadj <- res_bootstrap$var_boot_est_unadj
    var_boot_est_adj_theta1_true <- res_bootstrap$var_boot_est_adj_theta1_true
    var_boot_est_adj_theta1_per12 <- res_bootstrap$var_boot_est_adj_theta1_per12
    var_boot_est_adj_theta1_per1 <- res_bootstrap$var_boot_est_adj_theta1_per1
    var_boot_est_adj_theta1_per2 <- res_bootstrap$var_boot_est_adj_theta1_per2
    boot_stop <- res_bootstrap$boot_stop
    
    
    
    
    # Wald test for arm 2 (unadjusted)
    test_stat_unadj <- theta2_tilde_unadj/sqrt(var_boot_est_unadj)
    reject_h02_unadj <- test_stat_unadj > qnorm(1-alpha)
    reject_h02_unadj_cond <- test_stat_unadj > qnorm(1-alpha)
    
    
    # Wald test for arm 2 (adjusted, true theta1)
    test_stat_adj_theta1_true <- theta2_tilde_adj_theta1_true/sqrt(var_boot_est_adj_theta1_true)
    reject_h02_adj_theta1_true <- test_stat_adj_theta1_true > qnorm(1-alpha)
    reject_h02_adj_theta1_true_cond <- test_stat_adj_theta1_true > qnorm(1-alpha)
    
    # Wald test for arm 2 (adjusted, theta1 from periods 1 and 2)
    test_stat_adj_theta1_per12 <- theta2_tilde_adj_theta1_per12/sqrt(var_boot_est_adj_theta1_per12)
    reject_h02_adj_theta1_per12 <- test_stat_adj_theta1_per12 > qnorm(1-alpha)
    reject_h02_adj_theta1_per12_cond <- test_stat_adj_theta1_per12 > qnorm(1-alpha)
    
    # Wald test for arm 2 (adjusted, theta1 from period 1)
    test_stat_adj_theta1_per1 <- theta2_tilde_adj_theta1_per1/sqrt(var_boot_est_adj_theta1_per1)
    reject_h02_adj_theta1_per1 <- test_stat_adj_theta1_per1 > qnorm(1-alpha)
    reject_h02_adj_theta1_per1_cond <- test_stat_adj_theta1_per1 > qnorm(1-alpha)
    
    # Wald test for arm 2 (adjusted, theta1 from period 2)
    test_stat_adj_theta1_per2 <- theta2_tilde_adj_theta1_per2/sqrt(var_boot_est_adj_theta1_per2)
    reject_h02_adj_theta1_per2 <- test_stat_adj_theta1_per2 > qnorm(1-alpha)
    reject_h02_adj_theta1_per2_cond <- test_stat_adj_theta1_per2 > qnorm(1-alpha)
    
  }
  
  
  
  # Separate analysis for comparison
  
  x0 <- Data[Data$treatment==0 & Data$period_2==2,]$response
  x1 <- Data[Data$treatment==2 & Data$period_2==2,]$response
  
  # calculate the z-statistic
  theta2_sep <- mean(x1) - mean(x0)
  
  z_stat <- (theta2_sep - mu0) / 
    sqrt(sigma^2/n_22 + sigma^2/n_02)
  
  # metrics
  bias_arm2_sep <- (mean(x1) - mean(x0)) - theta[2]
  reject_h02_sep <- (z_stat > qnorm(1-alpha))
  
  
  
  return(list(Z_1 = Z_1,
              arm_1_stop = Z_1<c_1,
              arm_1_cont = Z_1>=c_1,
              
              boot_stop = boot_stop,
              
              theta2_tilde_unadj = theta2_tilde_unadj,
              theta2_tilde_adj_theta1_true = theta2_tilde_adj_theta1_true,
              theta2_tilde_adj_theta1_per12 = theta2_tilde_adj_theta1_per12,
              theta2_tilde_adj_theta1_per1 = theta2_tilde_adj_theta1_per1,
              theta2_tilde_adj_theta1_per2 = theta2_tilde_adj_theta1_per2,
              
              theta2_tilde_unadj_cond = theta2_tilde_unadj_cond,
              theta2_tilde_adj_theta1_true_cond = theta2_tilde_adj_theta1_true_cond,
              theta2_tilde_adj_theta1_per12_cond = theta2_tilde_adj_theta1_per12_cond,
              theta2_tilde_adj_theta1_per1_cond = theta2_tilde_adj_theta1_per1_cond,
              theta2_tilde_adj_theta1_per2_cond = theta2_tilde_adj_theta1_per2_cond,
              
              bias_est_theta1_true = bias_est_theta1_true,
              bias_est_theta1_per12 = bias_est_theta1_per12,
              bias_est_theta1_per1 = bias_est_theta1_per1,
              bias_est_theta1_per2 = bias_est_theta1_per2,
              
              var_boot_est_unadj = var_boot_est_unadj,
              var_boot_est_adj_theta1_true = var_boot_est_adj_theta1_true,
              var_boot_est_adj_theta1_per12 = var_boot_est_adj_theta1_per12,
              var_boot_est_adj_theta1_per1 = var_boot_est_adj_theta1_per1,
              var_boot_est_adj_theta1_per2 = var_boot_est_adj_theta1_per2,
              
              bias_arm2_unadj = bias_arm2_unadj,
              bias_arm2_adj_theta1_true = bias_arm2_adj_theta1_true,
              bias_arm2_adj_theta1_per12 = bias_arm2_adj_theta1_per12,
              bias_arm2_adj_theta1_per1 = bias_arm2_adj_theta1_per1,
              bias_arm2_adj_theta1_per2 = bias_arm2_adj_theta1_per2,
              
              bias_arm2_unadj_cond = bias_arm2_unadj_cond,
              bias_arm2_adj_theta1_true_cond = bias_arm2_adj_theta1_true_cond,
              bias_arm2_adj_theta1_per12_cond = bias_arm2_adj_theta1_per12_cond,
              bias_arm2_adj_theta1_per1_cond = bias_arm2_adj_theta1_per1_cond,
              bias_arm2_adj_theta1_per2_cond = bias_arm2_adj_theta1_per2_cond,
              
              MSE_arm2_unadj = bias_arm2_unadj^2,
              MSE_arm2_adj_theta1_true = bias_arm2_adj_theta1_true^2,
              MSE_arm2_adj_theta1_per12 = bias_arm2_adj_theta1_per12^2,
              MSE_arm2_adj_theta1_per1 = bias_arm2_adj_theta1_per1^2,
              MSE_arm2_adj_theta1_per2 = bias_arm2_adj_theta1_per2^2,
              
              MSE_arm2_unadj_cond = bias_arm2_unadj_cond^2,
              MSE_arm2_adj_theta1_true_cond = bias_arm2_adj_theta1_true_cond^2,
              MSE_arm2_adj_theta1_per12_cond = bias_arm2_adj_theta1_per12_cond^2,
              MSE_arm2_adj_theta1_per1_cond = bias_arm2_adj_theta1_per1_cond^2,
              MSE_arm2_adj_theta1_per2_cond = bias_arm2_adj_theta1_per2_cond^2,
              
              reject_h02_unadj = reject_h02_unadj, 
              reject_h02_adj_theta1_true = reject_h02_adj_theta1_true,
              reject_h02_adj_theta1_per12 = reject_h02_adj_theta1_per12,
              reject_h02_adj_theta1_per1 = reject_h02_adj_theta1_per1,
              reject_h02_adj_theta1_per2 = reject_h02_adj_theta1_per2,
              
              reject_h02_unadj_cond = reject_h02_unadj_cond, 
              reject_h02_adj_theta1_true_cond = reject_h02_adj_theta1_true_cond,
              reject_h02_adj_theta1_per12_cond = reject_h02_adj_theta1_per12_cond,
              reject_h02_adj_theta1_per1_cond = reject_h02_adj_theta1_per1_cond,
              reject_h02_adj_theta1_per2_cond = reject_h02_adj_theta1_per2_cond,
              
              theta2_tilde_unadj_stop = theta2_tilde_unadj_stop,
              bias_arm2_unadj_stop = bias_arm2_unadj_stop,
              MSE_arm2_unadj_stop = bias_arm2_unadj_stop^2,
              reject_h02_unadj_stop = reject_h02_unadj_stop,
              
              theta2_sep = theta2_sep,
              bias_arm2_sep = bias_arm2_sep,
              MSE_arm2_sep = bias_arm2_sep^2,
              reject_h02_sep = reject_h02_sep))
}

