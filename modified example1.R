# Install required packages (run once)
#install.packages("garma")
#install.packages("forecast")
#install.packages("tseries")

# Load libraries
library(garma)
library(forecast)
library(tseries)

###########################################################
# Example 1: AirPassengers (built-in)   y=detrend(log(AirPassengers)) 
# Load data Example 1: AirPassengers (built-in)
data("AirPassengers")
y <- AirPassengers  
# Plot
plot(y, main = "AirPassengers")

# Periodogram
spec.pgram(y, main = "Periodogram: AirPassengers")
abline(v=1, col="red", lwd=2)
acf(y, lag.max=100, main = "ACF: AirPassengers")

# Remove deterministic trend and seasonality
fit_season <- tslm(y ~ trend + season)
y_detrend <- residuals(fit_season)
plot(y_detrend, main = "AirPassengers")
spec.pgram(y_detrend, main = "Periodogram: detrend(AirPassengers)")
acf(y_detrend, lag.max=100, main = "ACF: detrend(AirPassengers)")


# Fit GARMA (1 Gegenbauer factor)
fit_garma <- garma(y_detrend, order = c(0,0,0), k = 1)

summary(fit_garma)

# Compare with ARIMA
auto.arima(y_detrend)
