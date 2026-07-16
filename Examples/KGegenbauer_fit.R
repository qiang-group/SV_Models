library(rstan)
library(quantmod)
library(latex2exp)
set.seed(100)

#function to compute Cauchy product of series
cauchy_prod <- function(a) {
  if (dim(a)[2] == 2) {
    lambda = numeric(dim(a)[1])
    for (k in 1:(dim(a)[1])) {
      lambda[k] = sum(a[1:k, 1] * a[k:1, 2])
    }
    return(lambda)
  }
  if (dim(a)[2] > 2) {
    return(cauchy_prod(cbind(cauchy_prod(a[, -1]), a[, 1])))
  }
  else{
    stop('Only one series is provided')
  }
}

#function to calculate coefficients of analytical expantion of \prod (1-2u_iB + B^2)^-d_i up tp order J
g_coef <- function(u, d, J) {
  if (length(u) != length(d)) {
    stop('The length of u(',
         length(u),
         ') does not match length of d(',
         length(d),
         ').')
  }
  lambda = matrix(0, nrow = J, ncol = length(d))
  lambda[1, ] = 1
  lambda[2, ] = 2 * u * d
  for (j in 3:J) {
    lambda[j, ] = 2 * u * ((d - 1) / (j-1) + 1) * lambda[j - 1, ] - (2 * (d - 1) /  (j-1) + 1 ) * lambda[j - 2, ]
    
  }
  if (length(u) == 1) {
    return(as.vector(lambda))
  }
  else{
    return(cauchy_prod(lambda))
  }
}


#Create model specification in Stan
code = code = r'(
functions{
  vector cauchy_prod(matrix a) {
    if (dims(a)[2] == 2) {
        vector[dims(a)[1]] lambda;
        lambda = to_vector(rep_array(0, dims(a)[1]));
      for (k in 1:(dims(a)[1])) {
        for (i in 1:k){
          lambda[k] += a[i, 1] * a[k-i+1, 2];
        }
      }
      return lambda;
    } else {
      return cauchy_prod(append_col(cauchy_prod(a[, 2:dims(a)[2]]), a[, 1]));
    }
  }
  
  vector g_coef(row_vector u, row_vector d, int J) {
    matrix[J + 1 ,size(u)] lambda;
    lambda[1, ] = to_row_vector(rep_array(1, size(u)));
    lambda[2, ] = 2 .* u .* d;
    for (j in 3:(J+1)) {
      lambda[j, ] = 2 * u .* ((d - 1) / (j-1) + 1) .* lambda[j - 1, ] - (2 * (d - 1) /  (j-1) + 1 ) .* lambda[j - 2, ];
    }
    if (size(u) == 1) {
      return to_vector(lambda);
    }
    else{
      return cauchy_prod(lambda);
    }
  }
}

data {
  int<lower=0> N;
  int<lower=1> k;
  vector[N] y;
  int<lower=3> J;
}

parameters {
  ordered[k] lu;
  row_vector<lower=0, upper=0.5>[k] d;
  real mu;
  real<lower=0> sigma;
  real eta[N];
}

transformed parameters{
  vector[J+1] lambda;
  real h[N-J];
  row_vector<lower=0, upper=1>[k] a;
  row_vector<lower=-1, upper=1>[k] u;
  u = 2/(1+exp(-lu'))-1;
  lambda = g_coef(u ,d , J);
  a = 2*d;
  for (i in 1:(N-J)){
    h[i] = mu;
    for (j in 0:J){
      h[i] += eta[i+J-j]*lambda[j+1];
    }
  }
}

model {
  a ~ beta(2,2);
  lu ~ normal(0,5);
  mu ~ normal(-5,2);
  sigma ~ gamma(5,10);
  eta ~ normal(0, sigma);
  for (i in 1:(N-J)){
    y[i+J] ~ normal(0, exp(h[i]/2));
  }
}

)'

#compile the model
mod = stan_model(model_code = code)

###############One-Factor GSV
# Simulate a time series y_t from Gegenbauer SV model
# N: sample size
# mu: mean of h
# sigma: standard error of eta_t
# u,d: coefficients in Gegenbauer polynomial
# J: Truncation value when approximating infinite long memory 
N = 2000
mu = log(0.0001)
sigma = 2
u = 0.8
d = 0.3
J = 50

y = numeric(N)
h = numeric(N)
lambda = numeric(J)

eta = rnorm(N+J,0,sigma)

lambda = g_coef(u, d, J+1)

for (t in 1:N){
  h[t] = sum(eta[(t+J):t]*lambda) + mu
  y[t] = rnorm(1,0,exp(h[t]/2))
}
#pdf(file = 'plot/GSV.pdf')
plot(h, type = 'l', cex.main = 2, cex.lab = 1.5)
acf(h, main = 'ACF: h', cex.main = 2, cex.lab = 1.5)
spec.pgram(h, main ='Periodogram: h', cex.main = 2, cex.lab = 1.5)
plot(y, type = 'l', cex.main = 1.5, cex.lab = 1.5)
acf(y, main = 'ACF: y', cex.main = 1.5, cex.lab = 1.5)
spec.pgram(y, main = 'Periodogram: y', cex.main = 1.5, cex.lab = 1.5)
plot(log(y^2), type = 'l', cex.main = 1.5, cex.lab = 1.5)
acf(log(y^2), main = TeX(r'(ACF: $\log(y^2)$)'), cex.main = 1.5, cex.lab = 1.5)
spec.pgram(log(y^2), main = TeX(r'(Periodogram: $\log(y^2)$)'), cex.main = 1.5, cex.lab = 1.5)
#dev.off()

test_data = list(
  N = length(y),
  k = 1,
  y = y,
  J = 30
)

#MCMC sampling
fit = sampling(mod, data = test_data, iter = 10000, warmup = 5000,chains = 1)
#compute posterior means of u, d and sigma
mean(extract(fit,'u')$u)
mean(extract(fit,'d')$d)
mean(extract(fit,'sigma')$sigma)
sd(extract(fit,'u')$u)
sd(extract(fit,'d')$d)
sd(extract(fit,'sigma')$sigma)

################Two-Factor GSV
# Simulate a time series y_t from Gegenbauer SV model
# N: sample size
# mu: mean of h
# sigma: standard error of eta_t
# u_i,d_i: coefficients in Gegenbauer polynomial
# J: Truncation value when approximating infinite long memory 
N = 2000
mu = log(0.0001)
sigma = 2
u = c(0.8, -0.5)
d = c(0.3,0.4)
J = 50

y = numeric(N)
h = numeric(N)
lambda = numeric(J)

eta = rnorm(N+J,0,sigma)

lambda = g_coef(u, d, J+1)

for (t in 1:N){
  h[t] = sum(eta[(t+J):t]*lambda) + mu
  y[t] = rnorm(1,0,exp(h[t]/2))
}

#pdf(file = 'plot/GSV2.pdf')
plot(h, type = 'l', cex.main = 2, cex.lab = 1.5)
acf(h, main = 'ACF: h', cex.main = 2, cex.lab = 1.5)
spec.pgram(h, main ='Periodogram: h', cex.main = 2, cex.lab = 1.5)
plot(y, type = 'l', cex.main = 1.5, cex.lab = 1.5)
acf(y, main = 'ACF: y', cex.main = 1.5, cex.lab = 1.5)
spec.pgram(y, main = 'Periodogram: y', cex.main = 1.5, cex.lab = 1.5)
plot(log(y^2), type = 'l', cex.main = 1.5, cex.lab = 1.5)
acf(log(y^2), main = TeX(r'(ACF: $\log(y^2)$)'), cex.main = 1.5, cex.lab = 1.5)
spec.pgram(log(y^2), main = TeX(r'(Periodogram: $\log(y^2)$)'), cex.main = 1.5, cex.lab = 1.5)
#dev.off()

test_data = list(
  N = length(y),
  k = 2,
  y = y,
  J = 30
)

fit2 = sampling(mod, data = test_data, iter = 10000, warmup = 5000,chains = 1)
apply(extract(fit2,'u')$u,2,mean)
apply(extract(fit2,'d')$d,2,mean)
mean(extract(fit2,'sigma')$sigma)
apply(extract(fit2,'u')$u,2,sd)
apply(extract(fit2,'d')$d,2,sd)
sd(extract(fit2,'sigma')$sigma)
