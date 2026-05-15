library(rstan)
library(quantmod)
library(latex2exp)
set.seed(100)

#Create model specification in Stan
code = r'(
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

mod = stan_model(model_code = code)

code = '
data {
  int<lower=0> N;
  vector[N] y;
  int<lower=3> J;
}

parameters {
  real<lower=0, upper=0.5> d;
  real mu;
  real<lower=0> sigma;
  real eta[N];
}

transformed parameters{
  real h[N-J];
  vector[J+1] lambda;
  lambda[1] = 1;
  real<lower=0, upper=1> b;
  b = 2*d;
  for (j in 1:J){
      lambda[j+1] = lambda[j]*(-d-j+1)/j;
  }
  for (i in 1:(N-J)){
    h[i] = mu + eta[i+J];
    for (j in 1:J){
      lambda[j+1] = -lambda[j]*(-d-j+1)/j;
      h[i] += eta[i+J-j]*lambda[j+1];
    }
  }
}

model {
  b ~ beta(2,2);
  mu ~ normal(-5,2);
  sigma ~ gamma(5,10);
  eta ~ normal(0, sigma);
  for (i in 1:(N-J)){
    y[i+J] ~ normal(0, exp(h[i]/2));
  }
}'

mod2 = stan_model(model_code = code)

#Case study
getSymbols("^GSPC", src = "yahoo", from = "2000-01-01", to = "2025-12-31")

# Daily log returns
ret <- dailyReturn(Cl(GSPC), type = "log")
ret <- na.omit(ret)
ret <- ret[ret!=0]
t = index(ret)

pdf(file = 'plot/sp500.pdf')
plot(t, ret, main = "Daily return of S&P 500", xlab = 'Date', ylab = 'r', type = 'l', cex.main = 2, cex.lab = 1.5)
acf(ret, lag.max=100, main = "ACF: S&P 500", cex.main = 2, cex.lab = 1.5)
spec.pgram(ret, main = "Periodogram: S&P 500", cex.main = 2, cex.lab = 1.5)
plot(t, ret^2, main = "Squared daily return of S&P 500", xlab = 'Date', ylab = TeX(r'($r^2$)'), type = 'l', cex.main = 2, cex.lab = 1.5)
acf(ret^2, lag.max=100, main =  TeX(r'(ACF: $($S&P 500$)^2$)'), cex.main = 2, cex.lab = 1.5)
spec.pgram(ret^2, main = TeX(r'(Periodogram: $($S&P 500$)^2$)'), cex.main = 2, cex.lab = 1.5)
plot(t, log(ret^2), main = "Log of squared daily return of S&P 500", xlab = 'Date', ylab = TeX(r'($y$)'), type = 'l', cex.main = 2, cex.lab = 1.5)
acf(log(ret^2), lag.max=100, main = "ACF: y", cex.main = 2, cex.lab = 1.5)
spec.pgram(log(ret^2), main = "Periodogram: y", cex.main = 2, cex.lab = 1.5)
dev.off()

# Realized volatility proxy and log transform
RV    <- ret^2
logRV <- log(RV)
logRV <- na.omit(logRV)

# Construct HAR regressors and estimate a simple HAR model
y    <- as.numeric(logRV)
n    <- length(y)
y_w  <- filter(y, rep(1/5, 5),  sides = 1)
y_m  <- filter(y, rep(1/22,22), sides = 1)

har_df <- data.frame(
  y      = y,
  y_lag1 = c(NA, y[-n]),
  y_w    = c(NA,as.numeric(y_w)[-n]),
  y_m    = c(NA,as.numeric(y_m)[-n])
)
har_df <- na.omit(har_df)
har_fit <- lm(y ~ y_lag1 + y_w + y_m, data = har_df)
summary(har_fit)
MSE = mean(har_fit$residuals^2)

# 
J = 30
r = as.numeric(ret)

test_data = list(
  N = length(r),
  k = 1,
  y = r,
  J = J
)

fit1 = sampling(mod, data = test_data, iter = 10000, warmup = 5000,chains = 1)
h1 = apply(extract(fit1,'h')$h,2,mean)
MSE1 = mean((h1 - 1.27 - log(ret[-(1:J)]^2))^2)

test_data = list(
  N = length(r),
  k = 2,
  y = r,
  J = J
)

fit2 = sampling(mod, data = test_data, iter = 10000, warmup = 5000,chains = 1)
h2 = apply(extract(fit2,'h')$h,2,mean)
MSE2 = mean((h2 - 1.27 - log(ret[-(1:J)]^2))^2)



fit3 = sampling(mod2, data = test_data, iter = 10000, warmup = 5000,chains = 1)
h3 = apply(extract(fit3,'h')$h,2,mean)
MSE3 = mean((h3 - 1.27 - log(ret[-(1:J)]^2))^2)
pdf(file = 'plot/sp500_fitted.pdf',height=5, width=8)
layout(matrix(1:2, nrow=2, byrow=TRUE),
       heights = c(0.95, 0.05))
par(mar=c(4,2.5,1,1))
plot(t,y,type = 'l',col = 'gray70', xlab = 'Date', ylab = 'y')
lines(t[-(1:22)],har_fit$fitted.values, col = 'gray40')
lines(t[-(1:J)],h1 - 1.27, col = 'gray30')
lines(t[-(1:J)],h2 - 1.27, col = 'gray20')
lines(t[-(1:J)],h3 - 1.27, col = 'gray10')
par(mar=c(0,0,0,0)) 
plot(1, type = "n", axes=FALSE, xlab="", ylab="") 
legend('top',c('Observed','HAR','1-factor','2-factor', 'ARFIMA'), 
       horiz=TRUE,col = c('gray70','gray40','gray30','gray20','gray10'),lty =c(1,1,1,1,1),
       lwd = c(1,1,1,1,1),cex = 0.8, bty = "n")
dev.off()