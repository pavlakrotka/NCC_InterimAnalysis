# Wrapper function for performing simulation studies for a given set of scenarios (parallelized on replication level)

# Arguments:

# scenario_df: Data frame containing all parameters for scenarios that should be simulated.
# nsim: Integer. Number of replications.
# perc_cores: Double. What percentage of available cores should be used for the simulations.

library(parallelly)
library(foreach)
library(iterators)
library(doFuture)
library(future)

sim_function_MAE_par <- function(scenarios_df, nsim, perc_cores = 0.99){
  
  cores <- availableCores()
  n_cores <- ifelse(floor(unname(cores)*perc_cores)<=1, 1, floor(unname(cores)*perc_cores)) # always use at least one core
  plan(strategy = "cluster", workers = n_cores)
  
  res_sim <- data.frame()
  
  # Create log file
  
  write(paste0("Starting the simulations. Time: ", Sys.time()), file = paste0("simlog_", deparse(substitute(scenarios_df)), ".txt"), append = FALSE)
  
  for (i in 1:nrow(scenarios_df)) {
    
    theta_i <- as.numeric(scenarios_df[i, grepl("^theta\\d", names(scenarios_df))])
    
    lambda_i <- as.numeric(scenarios_df[i, grepl("^lambda\\d", names(scenarios_df))])
    
    
    
    res <- foreach(icount(nsim), .combine = rbind, .options.future = list(seed = TRUE)) %dofuture% {
      runtrial_IA_MAE(n_01 = scenarios_df$n_01[i],
                      n_11 = scenarios_df$n_11[i],
                      n_02 = scenarios_df$n_02[i],
                      n_12 = scenarios_df$n_12[i],
                      n_22 = scenarios_df$n_22[i],
                      
                      mu0 = scenarios_df$mu[i],
                      theta = theta_i,
                      sigma = scenarios_df$sigma[i],
                      
                      futility_bound = scenarios_df$futility_bound[i],
                      alpha = scenarios_df$alpha[i],
                      
                      lambda = lambda_i,
                      trend_pattern = scenarios_df$trend_pattern[i],
                      period_blocks = scenarios_df$period_blocks[i],
                      
                      B_boot = scenarios_df$B_boot[i])
    }
    
    res_sim <- rbind(res_sim,
                     cbind(t(colMeans(matrix(as.numeric(res), nrow = nsim), na.rm = T)),
                           var(unlist(as.data.frame(res)$theta2_tilde_unadj)),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_true)),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_per12)),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_per1)),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_per2)),
                           var(unlist(as.data.frame(res)$theta2_tilde_unadj_cond), na.rm = T),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_true_cond), na.rm = T),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_per12_cond), na.rm = T),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_per1_cond), na.rm = T),
                           var(unlist(as.data.frame(res)$theta2_tilde_adj_theta1_per2_cond), na.rm = T),
                           var(unlist(as.data.frame(res)$theta2_tilde_unadj_stop), na.rm = T),
                           var(unlist(as.data.frame(res)$bias_est_theta1_true), na.rm = T),
                           var(unlist(as.data.frame(res)$bias_est_theta1_per12), na.rm = T),
                           var(unlist(as.data.frame(res)$bias_est_theta1_per1), na.rm = T),
                           var(unlist(as.data.frame(res)$bias_est_theta1_per2), na.rm = T),
                           
                           
                           sd(unlist(as.data.frame(res)$bias_arm2_unadj), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_true), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_per12), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_per1), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_per2), na.rm = T),
                           
                           sd(unlist(as.data.frame(res)$bias_arm2_unadj_cond), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_true_cond), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_per12_cond), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_per1_cond), na.rm = T),
                           sd(unlist(as.data.frame(res)$bias_arm2_adj_theta1_per2_cond), na.rm = T)))
    
    print(paste0("Scenario ", i, "/", nrow(scenarios_df), " done. Time: ", Sys.time()))
    
    line = paste0("Scenario ", i, "/", nrow(scenarios_df), " done. Time: ", Sys.time())
    write(line, file = paste0("simlog_", deparse(substitute(scenarios_df)), ".txt"), append = TRUE)
  }
  
  
  colnames(res_sim) <- c(colnames(res),
                         "var_theta2_tilde_unadj",
                         "var_theta2_tilde_adj_theta1_true",
                         "var_theta2_tilde_adj_theta1_per12",
                         "var_theta2_tilde_adj_theta1_per1",
                         "var_theta2_tilde_adj_theta1_per2",
                         "var_theta2_tilde_unadj_cond",
                         "var_theta2_tilde_adj_theta1_true_cond",
                         "var_theta2_tilde_adj_theta1_per12_cond",
                         "var_theta2_tilde_adj_theta1_per1_cond",
                         "var_theta2_tilde_adj_theta1_per2_cond",
                         "var_theta2_tilde_unadj_stop",
                         "var_bias_est_theta1_true",
                         "var_bias_est_theta1_per12",
                         "var_bias_est_theta1_per1",
                         "var_bias_est_theta1_per2",
                         "sd_bias_arm2_unadj" ,
                         "sd_bias_arm2_adj_theta1_true",
                         "sd_bias_arm2_adj_theta1_per12",
                         "sd_bias_arm2_adj_theta1_per1",
                         "sd_bias_arm2_adj_theta1_per2",
                         "sd_bias_arm2_unadj_cond",
                         "sd_bias_arm2_adj_theta1_true_cond",
                         "sd_bias_arm2_adj_theta1_per12_cond",
                         "sd_bias_arm2_adj_theta1_per1_cond",
                         "sd_bias_arm2_adj_theta1_per2_cond")
  
  
  plan(strategy = "sequential")
  gc()
  
  res_final <- cbind(scenarios_df, res_sim)
  res_final
}
