source("Report/Functions/kernel+mcmc_multimix.R")
source("Report/Functions/random_walk_MH.R")
source("Report/Functions/log gamma for simulation.R")

set.seed(447)

sample_sim1 <- mcmc_multmix(x0=list(c(-5,5)), N=10000, log_gamma=log_gamma_sim1, 
                            Gamma=list(diag(2)), beta=0.5, alpha=0.4, 
                            x_sigma=list(diag(2)), dat=rnorm(n=100, mean=0, sd=1))
saveRDS(sample_sim1, "Report/RData/sample_sim1.rds")



sms_data = c(13,24,8,24,7,35,14,11,15,11,22,22,11,57,11,19,29,6,19,
             12,22,12,18,72,32,9,7,13,19,23,27,20,6,17,13,10,14,6,
             16,15,7,2,15,15,19,70,49,7,53,22,21,31,19,11,18,20,12,
             35,17,23,17,4,2,31,30,13,27,0,39,37,5,14,13,22)
x0 <- c(rexp(1, rate = 1/100), 
        rexp(1, rate = 1/100), 
        extraDistr::rdunif(1, min = 1, max = length(sms_data)))
sample_sim2 <- mcmc_multmix(x0=list(c(70,70), 70), N=30000, 
                            log_gamma=log_gamma_sim2, Gamma=list(diag(2),30),                              
                            beta=0.5, alpha=0.4, x_sigma=list(diag(2),30), 
                            dat=sms_data)
saveRDS(sample_sim2, "Report/RData/sample_sim2.rds")
MH_sample_sim2 <- random_walk_MH(rates = list(70,70),
                                 change_point = 70,
                                 y=sms_data,
                                 n_iterations = 30000)
saveRDS(MH_sample_sim2, "Report/RData/MH_sample_sim2.rds")

ESS_summary <- matrix(0, nrow=3, ncol=4)
MCMC_result <- list()

for(j in 1:4){
  if(j==1){
    MCMC_result[[j]] <- mcmc_multmix(x0=list(70,70,70), N=30000, 
                                     log_gamma=log_gamma_sim2, Gamma=list(1,1,30), 
                                     beta=0.5, alpha=0.4, x_sigma=list(1,1,30), 
                                     dat=sms_data)
    ESS_1_rate1 <- MCMC_result[[j]]$samples[[1]][10001:30000]
    ESS_1_rate2 <- MCMC_result[[j]]$samples[[2]][10001:30000]
    ESS_1_CP <- MCMC_result[[j]]$samples[[3]][10001:30000]
    seconds <- MCMC_result[[j]]$total_time
    ESS_dataframe <- tibble(ESS_1_rate1,ESS_1_rate2, ESS_1_CP)
    ESS_summary[,j] = effectiveSize(as.mcmc(ESS_dataframe))/as.numeric(seconds)
  }
  if(j==2){
    MCMC_result[[j]] <- mcmc_multmix(x0=list(c(70,70),70), N=30000, 
                                     log_gamma=log_gamma_sim2, Gamma=list(diag(2),30), 
                                     beta=0.5, alpha=0.4, x_sigma=list(diag(2),30), 
                                     dat=sms_data)
    ESS_2_rate1 <- MCMC_result[[j]]$samples[[1]][10001:30000,1]
    ESS_2_rate2 <- MCMC_result[[j]]$samples[[1]][10001:30000,2]
    ESS_2_CP <- MCMC_result[[j]]$samples[[2]][10001:30000]
    seconds <- MCMC_result[[j]]$total_time
    ESS_dataframe <- tibble(ESS_2_rate1,ESS_2_rate2, ESS_2_CP)
    ESS_summary[,j] = effectiveSize(as.mcmc(ESS_dataframe))/as.numeric(seconds)
  }
  if(j==3){
    MCMC_result[[j]] <- mcmc_multmix(x0=list(c(70,70,70)), N=30000, 
                                     log_gamma=log_gamma_sim2,
                                     Gamma=list(diag(c(1,1,30))), 
                                     beta=0.5, alpha=0.4,
                                     x_sigma=list(diag(c(1,1,30))), 
                                     dat=sms_data)
    ESS_3_rate1 <- MCMC_result[[j]]$samples[[1]][10001:30000,1]
    ESS_3_rate2 <- MCMC_result[[j]]$samples[[1]][10001:30000,2]
    ESS_3_CP <- MCMC_result[[j]]$samples[[1]][10001:30000,3]
    seconds <- MCMC_result[[j]]$total_time
    ESS_dataframe <- tibble(ESS_3_rate1,ESS_3_rate2, ESS_3_CP)
    ESS_summary[,j] =
      effectiveSize(as.mcmc(ESS_dataframe))/as.numeric(seconds)
  }
  if(j==4){
    MCMC_result[[j]] <- random_walk_MH(rates = list(70,70),
                                       change_point = 70,
                                       y=sms_data,
                                       n_iterations = 30000)
    ESS_4_rate1 <- MCMC_result[[j]]$rates_trace[10001:30000,1]
    ESS_4_rate2 <- MCMC_result[[j]]$rates_trace[10001:30000,2]
    ESS_4_CP <- MCMC_result[[j]]$change_point_trace[10001:30000]
    seconds <- MCMC_result[[j]]$total_time
    ESS_dataframe <- tibble(ESS_4_rate1,ESS_4_rate2, ESS_4_CP)
    ESS_summary[,j] =
      effectiveSize(as.mcmc(ESS_dataframe))/as.numeric(seconds)
  }
}
colnames(ESS_summary) = c("No Covariance",
                          "Covariance between Two Exponential Parameterss", 
                          "Covariance among all Parameters",
                          "Traditional MH")
rownames(ESS_summary) = c("Exponential Rate before Change Point",
                          "Exponential Rate after Change Point",
                          "Change Point")

saveRDS(MCMC_result, "Report/RData/MCMC_results_sim2.rds")
saveRDS(ESS_summary, "Report/RData/ESS_summary.rds")
