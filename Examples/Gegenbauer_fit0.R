library(rstan)
set.seed(100)

# Simulate a time series y_t from Gegenbauer SV model
# N: sample size
# mu: mean of h
# sigma: standard error of eta_t
# u,j: coefficients in Gegenbauer polynomial
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

lambda[1] = 2*u*d
lambda[2] = 2*u*((d-1)/2+1)*lambda[1] - (2*(d-1)/2+1)
for (j in 3:J){
  lambda[j] = 2*u*((d-1)/j+1)*lambda[j-1] - (2*(d-1)/j+1)*lambda[j-2];
}

for (t in 1:N){
  h[t] = sum(eta[(t+J):t]*c(1,lambda)) + mu
  y[t] = rnorm(1,0,exp(h[t]/2))
}
acf(h)
pacf(h)

#Create model specification in Stan
code = '
data {
  int<lower=0> N;
  vector[N] y;
  int<lower=3> J;
}

parameters {
  real<lower=-1,upper=1> u;
  real<lower=0, upper=0.5> d;
  real mu;
  real<lower=0> sigma;
  real eta[N];
}

transformed parameters{
  real lambda[J];
  real h[N-J];
  real<lower=0, upper=1> a;
  real<lower=0, upper=1> b;
  lambda[1] = 2*u*d;
  lambda[2] = 2*u*((d-1)/2+1)*lambda[1] - (2*(d-1)/2+1);
  a = (1+u)/2;
  b = 2*d;
  for (j in 3:J){
    lambda[j] = 2*u*((d-1)/j+1)*lambda[j-1] - (2*(d-1)/j+1)*lambda[j-2];
  }
  for (i in 1:(N-J)){
    h[i] = eta[i+J]+mu;
    for (j in 1:J)
      h[i] += eta[i+J-j]*lambda[j];
  }
}

model {
  a ~ beta(2,2);
  b ~ beta(2,2);
  mu ~ normal(-5,2);
  sigma ~ gamma(5,10);
  eta ~ normal(0, sigma);
  for (i in 1:(N-J)){
    y[i+J] ~ normal(0, exp(h[i]/2));
  }
}'

#compile the model
mod = stan_model(model_code = code)

test_data = list(
  N = length(y),
  y = y,
  J = 30
)

#MCMC sampling
fit = sampling(mod, data = test_data, iter = 10000, warmup = 5000,chains = 1)
#compute posterior means of u, d and sigma
mean(extract(fit,'u')$u)
mean(extract(fit,'d')$d)
mean(extract(fit,'sigma')$sigma)
