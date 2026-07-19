gegenbauer_sim_k2 <- function(n=1500, u = c(0.7,0.3), d = c(0.35,0.25), sigma = 1, M = 100) {
  
# Compute Gegenbauer coefficients for one factor; M: truncating lag for filter
compute_psi <- function(u, d, M) {
  psi <- numeric(M)
  psi[1] <- 1
  psi[2] <- 2 * u * d
    
for (k in 3:M) {
  psi[k] <- (2 * u * (k - 2 + d) * psi[k - 1] - (k - 3 + 2 * d) * psi[k - 2]) / (k - 1)
  }
return(psi)
  }
  
# Compute psi for both factors
psi1 <- compute_psi(u[1], d[1], M)
psi2 <- compute_psi(u[2], d[2], M)
  
# Convolution of filters
psi_total <- convolve(psi1, rev(psi2), type = "open")
  
# Truncate to M (important for stability)
psi_total <- psi_total[1:M]
  
# White noise
e <- rnorm(n + M, sd = sigma)
  
# Apply filter
y <- filter(e, filter = psi_total, method = "convolution", sides = 1)
  
# Remove burn-in
y <- y[(M + 1):(M + n)]
return(as.numeric(y))
}

set.seed(123)
y <- gegenbauer_sim_k2(n = 2000, u = c(0.2, 0.7), d = c(0.4, 0.3), M = 200)