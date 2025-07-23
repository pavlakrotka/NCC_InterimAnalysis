library(tidyverse)

source("../functions/runtrial_IA_MAE.R")
source("../functions/bootstrap_functions.R")
source("../functions/get_mae.R")
source("../functions/sim_function_MAE_par.R")

n_sim <- 100000 # number of replications of each scenario

# NULL HYPOTHESIS - MAIN PAPER ----

## Varying alpha_1 ----

set.seed(2)
df_scenarios_t1e_alpha1 <- data.frame(n_01 = 150,
                                      n_11 = 150,
                                      n_02 = 150,
                                      n_12 = 150,
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = c(0.1, 0.15, 0.2, 0.25, 0.35, 0.50, 0.65, 0.75, 0.95),
                                      
                                      lambda0 = 0,
                                      lambda1 = 0,
                                      lambda2 = 0,
                                      trend_pattern = "linear",
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_t1e_alpha1 <- sim_function_MAE_par(df_scenarios_t1e_alpha1, nsim = n_sim)
write_csv(results_scenarios_t1e_alpha1, "results_final/results_scenarios_t1e_alpha1.csv")


## Varying r ----

set.seed(4)
df_scenarios_t1e_r <- data.frame(n_01 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_02 = 150,
                                 n_12 = 150,
                                 n_22 = 150,
                                 
                                 mu0 = 0,
                                 theta1 = 0,
                                 theta2 = 0,
                                 sigma = 1,
                                 
                                 alpha = 0.025,
                                 
                                 futility_bound = 0.5,
                                 
                                 lambda0 = 0,
                                 lambda1 = 0,
                                 lambda2 = 0,
                                 trend_pattern = "linear",
                                 
                                 period_blocks = 2,
                                 B_boot = 1000)

results_scenarios_t1e_r <- sim_function_MAE_par(df_scenarios_t1e_r, nsim = n_sim)
write_csv(results_scenarios_t1e_r, "results_final/results_scenarios_t1e_r.csv")


## Varying a ----

set.seed(9)
df_scenarios_t1e_a <- data.frame(n_01 = 150,
                                 n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_02 = 150,
                                 n_12 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_22 = 150,
                                 
                                 mu0 = 0,
                                 theta1 = 0,
                                 theta2 = 0,
                                 sigma = 1,
                                 
                                 alpha = 0.025,
                                 
                                 futility_bound = 0.5,
                                 
                                 lambda0 = 0,
                                 lambda1 = 0,
                                 lambda2 = 0,
                                 trend_pattern = "linear",
                                 
                                 period_blocks = 2,
                                 B_boot = 1000)

results_scenarios_t1e_a <- sim_function_MAE_par(df_scenarios_t1e_a, nsim = n_sim)
write_csv(results_scenarios_t1e_a, "results_final/results_scenarios_t1e_a.csv")


## Varying lambda ----

set.seed(5)
df_scenarios_t1e_lambda <- data.frame(n_01 = 150,
                                      n_11 = 150,
                                      n_02 = 150,
                                      n_12 = 150,
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = 0.5,
                                      
                                      lambda0 = seq(-0.15, 0.15, length.out=5),
                                      lambda1 = seq(-0.15, 0.15, length.out=5),
                                      lambda2 = seq(-0.15, 0.15, length.out=5),
                                      trend_pattern = "linear",
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_t1e_lambda <- sim_function_MAE_par(df_scenarios_t1e_lambda, nsim = n_sim)
write_csv(results_scenarios_t1e_lambda, "results_final/results_scenarios_t1e_lambda.csv")





# NULL HYPOTHESIS - SUPPLEMENTARY MATERIAL ----


## Varying alpha_1 ----

set.seed(16)
df_scenarios_t1e_alpha1_supp <- data.frame(n_01 = 150,
                                           n_11 = 150,
                                           n_02 = 150,
                                           n_12 = 150,
                                           n_22 = 150,
                                           
                                           mu0 = 0,
                                           theta1 = 0,
                                           theta2 = 0,
                                           sigma = 1,
                                           
                                           alpha = 0.025,
                                           
                                           futility_bound = c(0.1, 0.15, 0.2, 0.25, 0.50, 0.65, 0.75, 0.95),
                                           
                                           lambda0 = 0.15,
                                           lambda1 = 0.15,
                                           lambda2 = 0.15,
                                           trend_pattern = c(rep("linear", 8), rep("stepwise", 8)),
                                           
                                           period_blocks = 2,
                                           B_boot = 1000)

results_scenarios_t1e_alpha1_supp <- sim_function_MAE_par(df_scenarios_t1e_alpha1_supp, nsim = n_sim)
write_csv(results_scenarios_t1e_alpha1_supp, "results_final/results_scenarios_t1e_alpha1_supp.csv")


## Varying r ----

set.seed(42)
df_scenarios_t1e_r_supp <- data.frame(n_01 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_02 = 150,
                                      n_12 = 150,
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = 0.5,
                                      
                                      lambda0 = 0.15,
                                      lambda1 = 0.15,
                                      lambda2 = 0.15,
                                      trend_pattern = c(rep("linear", 7), rep("stepwise", 7)),
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_t1e_r_supp <- sim_function_MAE_par(df_scenarios_t1e_r_supp, nsim = n_sim)
write_csv(results_scenarios_t1e_r_supp, "results_final/results_scenarios_t1e_r_supp.csv")


## Varying a ----

set.seed(59)
df_scenarios_t1e_a_supp <- data.frame(n_01 = 150,
                                      n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_02 = 150,
                                      n_12 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = 0.5,
                                      
                                      lambda0 = 0.15,
                                      lambda1 = 0.15,
                                      lambda2 = 0.15,
                                      trend_pattern = c(rep("linear", 7), rep("stepwise", 7)),
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_t1e_a_supp <- sim_function_MAE_par(df_scenarios_t1e_a_supp, nsim = n_sim)
write_csv(results_scenarios_t1e_a_supp, "results_final/results_scenarios_t1e_a_supp.csv")


## Varying lambda ----

set.seed(67)
df_scenarios_t1e_lambda_supp <- data.frame(n_01 = 150,
                                           n_11 = 150,
                                           n_02 = 150,
                                           n_12 = 150,
                                           n_22 = 150,
                                           
                                           mu0 = 0,
                                           theta1 = 0,
                                           theta2 = 0,
                                           sigma = 1,
                                           
                                           alpha = 0.025,
                                           
                                           futility_bound = 0.5,
                                           
                                           lambda0 = seq(-0.15, 0.15, length.out=5),
                                           lambda1 = seq(-0.15, 0.15, length.out=5),
                                           lambda2 = seq(-0.15, 0.15, length.out=5),
                                           trend_pattern = rep(c("linear", "stepwise"), each = 5),
                                           
                                           period_blocks = 2,
                                           B_boot = 1000)

results_scenarios_t1e_lambda_supp <- sim_function_MAE_par(df_scenarios_t1e_lambda_supp, nsim = n_sim)
write_csv(results_scenarios_t1e_lambda_supp, "results_final/results_scenarios_t1e_lambda_supp.csv")





# ALTERNATIVE HYPOTHESIS  - MAIN PAPER ----

## Varying alpha_1 ----

set.seed(63)
df_scenarios_pow_alpha1 <- data.frame(n_01 = 150,
                                      n_11 = 150,
                                      n_02 = 150,
                                      n_12 = 150,
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0.32,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = c(0.1, 0.15, 0.2, 0.25, 0.50, 0.65, 0.75, 0.95),
                                      
                                      lambda0 = 0,
                                      lambda1 = 0,
                                      lambda2 = 0,
                                      trend_pattern = "linear",
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_pow_alpha1 <- sim_function_MAE_par(df_scenarios_pow_alpha1, nsim = n_sim)
write_csv(results_scenarios_pow_alpha1, "results_final/results_scenarios_pow_alpha1.csv")


## Varying r ----

set.seed(42)
df_scenarios_pow_r <- data.frame(n_01 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_02 = 150,
                                 n_12 = 150,
                                 n_22 = 150,
                                 
                                 mu0 = 0,
                                 theta1 = 0,
                                 theta2 = 0.32,
                                 sigma = 1,
                                 
                                 alpha = 0.025,
                                 
                                 futility_bound = 0.5,
                                 
                                 lambda0 = 0,
                                 lambda1 = 0,
                                 lambda2 = 0,
                                 trend_pattern = "linear",
                                 
                                 period_blocks = 2,
                                 B_boot = 1000)


results_scenarios_pow_r <- sim_function_MAE_par(df_scenarios_pow_r, nsim = n_sim)
write_csv(results_scenarios_pow_r, "results_final/results_scenarios_pow_r.csv")


## Varying a ----

set.seed(95)
df_scenarios_pow_a <- data.frame(n_01 = 150,
                                 n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_02 = 150,
                                 n_12 = c(10, 50, 150, 300, 600, 1050, 1500),
                                 n_22 = 150,
                                 
                                 mu0 = 0,
                                 theta1 = 0,
                                 theta2 = 0.32,
                                 sigma = 1,
                                 
                                 alpha = 0.025,
                                 
                                 futility_bound = 0.5,
                                 
                                 lambda0 = 0,
                                 lambda1 = 0,
                                 lambda2 = 0,
                                 trend_pattern = "linear",
                                 
                                 period_blocks = 2,
                                 B_boot = 1000)

results_scenarios_pow_a <- sim_function_MAE_par(df_scenarios_pow_a, nsim = n_sim)
write_csv(results_scenarios_pow_a, "results_final/results_scenarios_pow_a.csv")


## Varying lambda ----

set.seed(57)
df_scenarios_pow_lambda <- data.frame(n_01 = 150,
                                      n_11 = 150,
                                      n_02 = 150,
                                      n_12 = 150,
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0.32,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = 0.5,
                                      
                                      lambda0 = seq(-0.15, 0.15, length.out=5),
                                      lambda1 = seq(-0.15, 0.15, length.out=5),
                                      lambda2 = seq(-0.15, 0.15, length.out=5),
                                      trend_pattern = "linear",
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_pow_lambda <- sim_function_MAE_par(df_scenarios_pow_lambda, nsim = n_sim)
write_csv(results_scenarios_pow_lambda, "results_final/results_scenarios_pow_lambda.csv")





# ALTERNATIVE HYPOTHESIS - SUPPLEMENTARY MATERIAL ----

## Varying alpha_1 ----

set.seed(653)
df_scenarios_pow_alpha1_supp <- data.frame(n_01 = 150,
                                           n_11 = 150,
                                           n_02 = 150,
                                           n_12 = 150,
                                           n_22 = 150,
                                           
                                           mu0 = 0,
                                           theta1 = 0,
                                           theta2 = 0.32,
                                           sigma = 1,
                                           
                                           alpha = 0.025,
                                           
                                           futility_bound = c(0.1, 0.15, 0.2, 0.25, 0.50, 0.65, 0.75, 0.95),
                                           
                                           lambda0 = 0.15,
                                           lambda1 = 0.15,
                                           lambda2 = 0.15,
                                           trend_pattern = c(rep("linear", 8), rep("stepwise", 8)),
                                           
                                           period_blocks = 2,
                                           B_boot = 1000)

results_scenarios_pow_alpha1_supp <- sim_function_MAE_par(df_scenarios_pow_alpha1_supp, nsim = n_sim)
write_csv(results_scenarios_pow_alpha1_supp, "results_final/results_scenarios_pow_alpha1_supp.csv")


## Varying r ----

set.seed(346)
df_scenarios_pow_r_supp <- data.frame(n_01 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_02 = 150,
                                      n_12 = 150,
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0.32,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = 0.5,
                                      
                                      lambda0 = 0.15,
                                      lambda1 = 0.15,
                                      lambda2 = 0.15,
                                      trend_pattern = c(rep("linear", 7), rep("stepwise", 7)),
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_pow_r_supp <- sim_function_MAE_par(df_scenarios_pow_r_supp, nsim = n_sim)
write_csv(results_scenarios_pow_r_supp, "results_final/results_scenarios_pow_r_supp.csv")


## Varying a ----

set.seed(395)
df_scenarios_pow_a_supp <- data.frame(n_01 = 150,
                                      n_11 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_02 = 150,
                                      n_12 = c(10, 50, 150, 300, 600, 1050, 1500),
                                      n_22 = 150,
                                      
                                      mu0 = 0,
                                      theta1 = 0,
                                      theta2 = 0.32,
                                      sigma = 1,
                                      
                                      alpha = 0.025,
                                      
                                      futility_bound = 0.5,
                                      
                                      lambda0 = 0.15,
                                      lambda1 = 0.15,
                                      lambda2 = 0.15,
                                      trend_pattern = c(rep("linear", 7), rep("stepwise", 7)),
                                      
                                      period_blocks = 2,
                                      B_boot = 1000)

results_scenarios_pow_a_supp <- sim_function_MAE_par(df_scenarios_pow_a_supp, nsim = n_sim)
write_csv(results_scenarios_pow_a_supp, "results_final/results_scenarios_pow_a_supp.csv")


## Varying lambda ----

set.seed(857)
df_scenarios_pow_lambda_supp <- data.frame(n_01 = 150,
                                           n_11 = 150,
                                           n_02 = 150,
                                           n_12 = 150,
                                           n_22 = 150,
                                           
                                           mu0 = 0,
                                           theta1 = 0,
                                           theta2 = 0.32,
                                           sigma = 1,
                                           
                                           alpha = 0.025,
                                           
                                           futility_bound = 0.5,
                                           
                                           lambda0 = seq(-0.15, 0.15, length.out=5),
                                           lambda1 = seq(-0.15, 0.15, length.out=5),
                                           lambda2 = seq(-0.15, 0.15, length.out=5),
                                           trend_pattern = rep(c("linear", "stepwise"), each = 5),
                                           
                                           period_blocks = 2,
                                           B_boot = 1000)

results_scenarios_pow_lambda_supp <- sim_function_MAE_par(df_scenarios_pow_lambda_supp, nsim = n_sim)
write_csv(results_scenarios_pow_lambda_supp, "results_final/results_scenarios_pow_lambda_supp.csv")



































