gegenbauer_sim <- function(n = 1500, u = 0.7, d = 0.35, sigma = 1, M = 50) {

# White noise;  M: truncation lag for filter
e <- rnorm(n + M, sd = sigma)
  
  # Compute Gegenbauer coefficients recursively
  psi <- numeric(M)
  psi[1] <- 1
  psi[2] <- 2*u*d
  
  for (k in 3:M) {
    psi[k] <- (2*u*(k-2+d)*psi[k-1] - (k-3+2*d)*psi[k-2]) / (k-1)
  }
  
  # Apply filter (convolution)
  y <- filter(e, filter = psi, method = "convolution", sides = 1)
  
  # Remove burn-in
  y <- y[(M+1):(M+n)]
  return(as.numeric(y))
}

set.seed(123)
y <- gegenbauer_sim_spec(n = 1500, u = 0.7, d = 0.35)