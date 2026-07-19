# Load required packages
library(quantmod)
library(xts)

# Download S&P 500 data
getSymbols("^GSPC", src = "yahoo", from = "2000-01-01", to = "2025-12-31")

# Daily log returns
y <- dailyReturn(Cl(GSPC), type = "log")
y <- na.omit(y)
y <- y[y!=0]
y <- index(y)
y <- as.numeric(y)

# Realized volatility proxy and log transform
z <- log(y^2)

# Construct HAR regressors and estimate a simple HAR model
z    <- as.numeric(logRV)
n    <- length(z)
z_w  <- filter(z, rep(1/5, 5),  sides = 1)
z_m  <- filter(z, rep(1/22,22), sides = 1)

har_df <- data.frame(
  z      = z,
  z_lag1 = c(NA, z[-n]),
  z_w    = c(NA,as.numeric(z_w)[-n]),
  z_m    = c(NA,as.numeric(z_m)[-n])
)
har_df <- na.omit(har_df)
har_fit <- lm(z ~ z_lag1 + z_w + z_m, data = har_df)
summary(har_fit)
MSE = mean(har_fit$residuals^2); MSE