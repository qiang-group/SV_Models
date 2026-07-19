install.packages("garma");  install.packages("forecast");  install.packages("tseries")
library(garma); library(forecast); library(tseries)

# Example 1: AirPassengers (built-in)
# Load data Example 1: AirPassengers (built-in)
data("AirPassengers")
y <- AirPassengers   # stabilize variance and remove trend

# Remove deterministic seasonality
fit_season <- tslm(y ~ trend + season)
y_detrend <- residuals(fit_season)

# Plot data, ACF and Periodogram 
plot(y_detrend, main = "detrend(AirPassengers)")
acf(y_detrend, lag.max=100, main = "ACF: detrend(diff log AirPassengers)")
spec.pgram(as.numeric(y_detrend), main = "Periodogram: detrend(AirPassengers)")

# Fit GARMA (1 Gegenbauer factor)
fit_garma <- garma(y_detrend, order = c(0,0,0), k = 1)
summary(fit_garma)