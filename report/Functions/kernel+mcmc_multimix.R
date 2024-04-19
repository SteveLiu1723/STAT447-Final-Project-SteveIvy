# a kernel function
kernel <- function(Gamma, x, x_sigma, beta, d){
  kernel_choice <- rbinom(1,1,prob=beta)
  if(kernel_choice==0){
    z <- MASS::mvrnorm(n = 1, mu=rep(0,d), Sigma = (2.38)^2 * x_sigma/d)
  }
  else{
    z <- MASS::mvrnorm(n = 1, mu=rep(0,d), Sigma = 0.1 * Gamma/d)
  }
  x_prime = x + z
  # return the next state
  return(x_prime)
}
# a MCMC algorithm which can take multiple sets of parameters
mcmc_multmix <- function(x0, N, log_gamma, Gamma, beta, alpha, x_sigma, dat){
  d_list <- sapply(x0, length) # dim of each parameter
  d <- sum(d_list)
  xbar <- x0
  x <- x0
  samples <- lapply(d_list, function(d){matrix(data=0, nrow=N, ncol=d)})
  start_time <- Sys.time()
  # start the iteration
  for(t in 1:N){
    Wn <- (alpha - 0.4)/sqrt(t)
    Gamma <- lapply(Gamma, function(gamm){(1 + Wn) * gamm})
    x_prime <- list()
    
    # get new samples
    for(i in 1:length(d_list)){
      x_prime[[i]] <- kernel(Gamma[[i]], x[[i]], x_sigma[[i]], beta, d_list[i])
    }
    
    #MH accept/reject
    alpha <- min(c(1, exp(log_gamma(unlist(x_prime),dat)-
                            log_gamma(unlist(x), dat))))
    if(runif(1) < alpha){
      x = x_prime
    }
    
    # store new sample
    for(i in 1:length(d_list)){
      samples[[i]][t,] <- x[[i]]
    }
    
    # update covariance matrix using Sherman-Morrison formula
    if(t == 2){
      x_sigma = lapply(samples, function(sample){var(sample[1:2, ])})
    } else if(t > 2){
      x_sigma <- mapply(function(sigma, mu, y) {
        (t-1)/t * sigma + (t-1)/t^2 * (mu - y) %*% t(mu - y)
      }, x_sigma, xbar, x, SIMPLIFY = FALSE)
    }
    
    # update xbar
    xbar = mapply(function(mu, y) {((t-1)*mu + y)/t}, 
                  xbar, x, SIMPLIFY = FALSE)
    
    if(t==10000){
      start_time <- Sys.time()
    }
  }
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "secs")
  return(list(samples=samples,
              total_time=total_time))
}