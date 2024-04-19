log_gamma_sim1 = function(para, data) {
  # Return log(0.0) if parameters are outside of the support
  if (para[2]<0) 
    return(-Inf)
  log_prior = 
    dnorm(para[1], mean=0, sd=1, log = TRUE) + 
    dexp(para[2], 1/100, log = TRUE)
  log_likelihood = 0.0
  for (j in 1:length(data)) {
    log_likelihood = log_likelihood + dnorm(data[j], mean=para[1],
                                            sd = para[2], log = TRUE)
  }
  return(log_prior + log_likelihood)
}

log_joint = function(rates, change_point, y) {
  # Return log(0.0) if parameters are outside of the support
  if (rates[[1]] < 0 | rates[[2]] < 0 | change_point < 1 | change_point > length(y)) 
    return(-Inf)
  
  log_prior = 
    dexp(rates[[1]], 1/100, log = TRUE) + 
    dexp(rates[[2]], 1/100, log = TRUE)
  
  log_likelihood = 0.0
  for (i in 1:length(y)) {
    rate = if (i < change_point) rates[[1]] else rates[[2]]
    log_likelihood = log_likelihood + dpois(y[[i]], rate, log = TRUE)
  }
  
  return(log_prior + log_likelihood)
}

log_gamma_sim2 = function(para, data) {
  # Return log(0.0) if parameters are outside of the support
  if (para[1] < 0 | para[2] < 0 | para[3] < 1 | para[3] > length(data)) 
    return(-Inf)
  log_prior = 
    dexp(para[1], 1/100, log = TRUE) + dexp(para[2], 1/100, log = TRUE)
  log_likelihood = 0.0
  for (j in 1:length(data)) {
    rate = if (j < para[3]) para[1] else para[2]
    log_likelihood = log_likelihood + dpois(data[j], rate, log = TRUE)
  }
  
  return(log_prior + log_likelihood)
}