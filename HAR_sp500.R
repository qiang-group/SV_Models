# Load required packages
library(quantmod)
library(xts)

# Download S&P 500 data
getSymbols("^GSPC", src = "yahoo", from = "2000-01-01")
spx <- GSPC

# Daily log returns
ret <- dailyReturn(Cl(spx), type = "log")
ret <- na.omit(ret)

# Realized volatility proxy and log transform
RV    <- ret^2
logRV <- log(RV+0.001)
logRV <- na.omit(logRV)

# Construct HAR regressors and estimate a simple HAR model
y    <- as.numeric(logRV)
n    <- length(y)
y_w  <- filter(y, rep(1/5, 5),  sides = 1)
y_m  <- filter(y, rep(1/22,22), sides = 1)

har_df <- data.frame(
  y      = y,
  y_lag1 = c(NA, y[-n]),
  y_w    = as.numeric(y_w),
  y_m    = as.numeric(y_m)
)
har_df <- na.omit(har_df)
har_fit <- lm(y ~ y_lag1 + y_w + y_m, data = har_df)
summary(har_fit)