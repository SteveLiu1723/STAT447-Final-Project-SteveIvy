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

MCMC_list1 <- list()
MCMC_list2 <- list()
MCMC_list3 <- list()
MCMC_list4 <- list()
ESS_summary <- matrix(nrow=3, ncol=4)
for(i in 1:10){
  MCMC_list1[[i]] <- mcmc_multmix(x0=list(70,70,70), N=30000, 
                                  log_gamma=log_gamma_sim2, Gamma=list(1,1,30), 
                                  beta=0.5, alpha=0.4, x_sigma=list(1,1,30), 
                                  dat=sms_data)
  MCMC_list2[[i]] <- mcmc_multmix(x0=list(c(70,70),70), N=30000, 
                                      log_gamma=log_gamma_sim2, Gamma=list(diag(2),30), 
                                      beta=0.5, alpha=0.4, x_sigma=list(diag(2),30), 
                                      dat=sms_data)
  MCMC_list3[[i]] <- mcmc_multmix(x0=list(c(70,70,70)), N=30000, 
                                  log_gamma=log_gamma_sim2,
                                  Gamma=list(diag(c(1,1,30))), 
                                  beta=0.5, alpha=0.4,
                                  x_sigma=list(diag(c(1,1,30))), 
                                  dat=sms_data)
  MCMC_list4[[i]] <- random_walk_MH(rates = list(70,70),
                                    change_point = 70,
                                    y=sms_data,
                                    n_iterations = 30000)
}

saveRDS(MCMC_list1, "Report/RData/MCMC_list1.rds")
saveRDS(MCMC_list2, "Report/RData/MCMC_list2.rds")
saveRDS(MCMC_list3, "Report/RData/MCMC_list3.rds")
saveRDS(MCMC_list4, "Report/RData/MCMC_list4.rds")

ESS_samples1 <- list()
time_1 <- list()
ESS_samples2 <- list()
time_2 <- list()
ESS_samples3 <- list()
time_3 <- list()
for(i in 1:10){
  ESS_samples1[[i]] <- effectiveSize(as.mcmc(tibble(MCMC_list1[[i]]$samples[[1]][10001:30000],
                           MCMC_list1[[i]]$samples[[2]][10001:30000],
                           MCMC_list1[[i]]$samples[[3]][10001:30000])))
  time_1[[i]] <- as.numeric(MCMC_list1[[i]]$total_time)
  ESS_samples2[[i]]<- effectiveSize(as.mcmc(tibble(MCMC_list2[[i]]$samples[[1]][10001:30000,1],
                           MCMC_list2[[i]]$samples[[1]][10001:30000,2],
                           MCMC_list2[[i]]$samples[[2]][10001:30000])))
  time_2[[i]] <- as.numeric(MCMC_list2[[i]]$total_time)
  ESS_samples3[[i]] <- effectiveSize(as.mcmc(tibble(MCMC_list3[[i]]$samples[[1]][10001:30000,1],
                           MCMC_list3[[i]]$samples[[1]][10001:30000,2],
                           MCMC_list3[[i]]$samples[[1]][10001:30000,3])))
  time_3[[i]] <- as.numeric(MCMC_list3[[i]]$total_time)
}
ratioList1 <- vector("list", length(ESS_samples1))
ratioList2 <- vector("list", length(ESS_samples2))
ratioList3 <- vector("list", length(ESS_samples3))
# Compute the ratio for each element
for (i in 1:10) {
  ratioList1[[i]] <- ESS_samples1[[i]] / time_1[[i]]
  ratioList2[[i]] <- ESS_samples2[[i]] / time_2[[i]]
  ratioList3[[i]] <- ESS_samples3[[i]] / time_3[[i]]
}
ESS_summary <- matrix(nrow=3,ncol=3)
mean_vector_1 <- sapply(1:length(ratioList1[[1]]), function(i) {
  mean(sapply(ratioList1, `[`, i))
})
mean_vector_2 <- sapply(1:length(ratioList2[[1]]), function(i) {
  mean(sapply(ratioList2, `[`, i))
})
mean_vector_3 <- sapply(1:length(ratioList3[[1]]), function(i) {
  mean(sapply(ratioList3, `[`, i))
})
colnames(ESS_summary) <- c("Rate before Change Point", 
                           "Rate after Change Point", 
                           "Change Point")
rownames(ESS_summary) <- c("BloRA: No Covariance",
                           "BloRA: Covariance between Rates",
                           "BloRA: Covariance among all Parameters")
ESS_summary[,1] <- mean_vector_1
ESS_summary[,2] <- mean_vector_2
ESS_summary[,3] <- mean_vector_3
ESS_summary <- t(ESS_summary)
saveRDS(ESS_summary, "Report/RData/ESS_summary.rds")

source("Report/Functions/compute_parameter_mcse.R")

list_of_mcmc_lists <- list(mcmc1 = MCMC_list1, mcmc2 = MCMC_list2, 
                           mcmc3 = MCMC_list3, mcmc4 = MCMC_list4)
results <- lapply(list_of_mcmc_lists, compute_parameter_mcse)

param_names <- c("Rate before Change Point", 
                 "Rate after Change Point", 
                 "Change Point")

results_df <- do.call(rbind, lapply(results, function(x) setNames(as.data.frame(t(x)), param_names)))
rownames(results_df) <- c("BloRA: No Covariance",
                          "BloRA: Covariance between Rates",
                          "BloRA: Covariance among all Parameters",
                          "Traditional MH")
saveRDS(results_df, "Report/RData/MCSE_table.rds")
