random_walk_MH = function(rates, change_point, y, n_iterations) {
  change_point_trace = rep(-1, n_iterations)
  rates_trace = matrix(0, nrow=n_iterations, ncol=2)
  # define a kernel function
  kernel = function(gamma, rates, change_point, kernel_index_choice) {
    # update parameters based on K_1
    if(kernel_index_choice==1){
      proposal_1 <- list(rnorm(n=1,mean=rates[[1]],
                               sd=1), 
                         rnorm(n=1,mean=rates[[2]],
                               sd=1))
      proposal_2 <- change_point
    }
    # update parameters based on K_2
    if(kernel_index_choice==2){
      proposal_1 <- rates
      proposal_2 <- sample(x=seq(change_point-15,
                                 change_point+15, by=1), size=1)
    }
    # compute difference in log_joint
    log_ratio = gamma(proposal_1, proposal_2, y)-gamma(rates, change_point, y) 
    if (runif(1) < exp(log_ratio)) {
      return(list(proposal_1, proposal_2))
    } else {
      return(list(rates, change_point))
    }
  }
  start_time <- Sys.time()
  # iteratively generate sampled change point
  for (i in 1:n_iterations) {
    kernel_index_choice = if (runif(1) < 0.5) 1 else 2
    next_state <- kernel(log_joint, rates, change_point, kernel_index_choice)
    rates <- next_state[[1]]
    rates_trace[i,] <- unlist(next_state[[1]])
    change_point <- next_state[[2]]
    change_point_trace[i] <- next_state[[2]]
    if(i==10000){
      start_time <- Sys.time()
    }
  }
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "secs")
  
  # Return:
  # - the trace of the change points (for question 1) 
  # - the rates at the last iteration (for question 2)
  return(
    list(
      change_point_trace = change_point_trace, 
      rates_trace = rates_trace,
      total_time=total_time
    )
  )
}