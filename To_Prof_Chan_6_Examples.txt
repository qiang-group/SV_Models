###############################################################
# Time Series Analysis using Gegenbauer Processes
#
# This project analyses several time series datasets to investigate
# cyclical long-memory behaviour using Gegenbauer models. The datasets
# include both real-world time series and a simulated Gegenbauer process.
#
# Datasets used:
# 1) AirPassengers – Monthly airline passengers (1949–1960).
#    Shows strong trend and seasonality.
#
# 2) NottinghamTemperatures – Monthly air temperatures in Nottingham, UK.
#    Displays strong seasonal cycles.
#
# 3) IndustrialProductionIndex – U.S. Industrial Production Index
#    obtained from FRED (Federal Reserve Economic Data).
#
# 4) SouthernOscillationIndex – Climate index related to El Niño /
#    La Niña atmospheric pressure cycles.
#
# 5) Sunspot – Monthly sunspot numbers reflecting solar activity
#    (~11-year solar cycle).
#
# 6) SimulatedGegenbauer – Artificial series generated with
#    parameters u = 0.7 and d = 0.35.
#
###############################################################

###############################################################
# 2. Key Steps in the Analysis
#
# Step 1: Data Preparation
# All datasets are stored as time series (ts objects) to preserve
# temporal ordering and seasonal frequency.
#
# Step 2: Data Transformations
# For each dataset three versions of the data are analysed:
#
#   - Raw data
#   - Log-transformed data (if values are positive)
#   - De-trended data (seasonal components removed using tslm)
#
# Step 3: Exploratory Analysis
# For each transformation we generate:
#   - Time series plot
#   - Autocorrelation Function (ACF)
#   - Periodogram
#
# These diagnostics help identify persistence, seasonality,
# and dominant cyclical frequencies in the data.
#
# Step 4: Gegenbauer Parameter Estimation
# Semi-parametric estimation is used to estimate:
#
#   u : cyclical frequency parameter
#   d : long-memory parameter
#
# The parameter u identifies the location of spectral peaks,
# while d measures persistence around that frequency.
#
# Step 5: Short-term Dynamics
# An ARMA(1,1) model is fitted to capture short-run dynamics,
# producing estimates of:
#
#   AR1
#   MA1
#   innovation variance (sigma^2)
#
# Step 6: Forecasting
# A one-step-ahead forecast is generated for each dataset
# based on the fitted model.
#
###############################################################


###############################################################
# 3. Main Output
#
# The final output of the analysis is a summary table containing
# parameter estimates and forecasts for each dataset and
# transformation.
#
# The table includes:
#
#   Example        : dataset index
#   Dataset        : dataset name
#   Transformation : raw / log / detrended
#   u              : Gegenbauer frequency estimate
#   d              : long memory parameter
#   AR1            : AR(1) coefficient
#   MA1            : MA(1) coefficient
#   Sigma2         : innovation variance
#   Forecast_t1    : one-step-ahead forecast
#
# This summary allows comparison of cyclical long-memory
# behaviour across the different datasets.
#
###############################################################
# Load Required Packages
###############################################################
# if needed below
install.packages(
  c("quantmod","astsa","ggbr","garma","forecast","tseries","ggplot2","gridExtra"),
  dependencies = TRUE
)

library(quantmod)
library(astsa)
#library(ggbr)
library(garma)
library(forecast)
library(tseries)
library(ggplot2)
library(gridExtra)


###############################################################
# Gegenbauer Estimation Function
###############################################################

estimate_gegenbauer <- function(y, log_transform = FALSE, diff_order = 0){
  
  y <- as.numeric(y)
  
  if(log_transform && min(y) > 0) y <- log(y)
  if(diff_order == 1) y <- diff(y)
  if(diff_order == 12) y <- diff(y, lag = 12)
  
  sp <- tryCatch(ggbr_semipara(y), error=function(e) NULL)
  
  if(is.null(sp)){
    return(list(u=NA,d=NA,AR=NA,MA=NA,Sigma2=NA,fit=NULL))
  }
  
  gf <- sp[[1]]
  
  epsilon <- 1e-3
  u <- max(min(gf$freq,1-epsilon),0.01)
  d <- gf$fd
  
  fit <- tryCatch(
    garma(y, order=c(1,0,1), k=0, method="CSS", include.mean=FALSE),
    error=function(e) NULL
  )
  
  if(is.null(fit)){
    AR <- NA; MA <- NA; Sigma2 <- NA
  } else {
    AR <- round(fit$coef[1],4)
    MA <- round(fit$coef[2],4)
    Sigma2 <- round(fit$sigma2,4)
  }
  
  list(
    u = round(u,4),
    d = round(d,4),
    AR = AR,
    MA = MA,
    Sigma2 = Sigma2,
    fit = fit
  )
}


###############################################################
# Simulated Gegenbauer Series Generator
###############################################################

simulate_gegenbauer <- function(n=1500,u=0.7,d=0.35,sigma=1){
  
  j <- 0:(n-1)
  lambda <- 2*pi*j/n
  
  f <- (1 - 2*u*cos(lambda) + u^2)^(-2*d)
  
  z <- complex(real=rnorm(n),imaginary=rnorm(n))
  y_fft <- sqrt(f)*z
  
  y <- Re(fft(y_fft,inverse=TRUE))/sqrt(n)
  y <- sigma*scale(y)
  
  as.numeric(y)
}


###############################################################
# Detrending Function
###############################################################

detrend_season <- function(x){
  fit <- tslm(x ~ season)
  residuals(fit)
}


###############################################################
# Plot Function (3x2 ACF + Periodogram)
###############################################################

plot_series_acf_periodogram <- function(y,dataset_name){
  
  series_list <- list(
    "Raw" = y,
    "Log" = if(min(y)>0) log(y) else y,
    "De-trended" = detrend_season(y)
  )
  
  plot_list <- list()
  
  for(i in seq_along(series_list)){
    
    sname <- names(series_list)[i]
    sdata <- series_list[[i]]
    
    acf_obj <- acf(sdata,plot=FALSE)
    
    acf_plot <- autoplot(acf_obj) +
      ggtitle(paste(dataset_name,"-",sname,"ACF")) +
      theme_minimal()
    
    spec <- spec.pgram(sdata,plot=FALSE)
    
    df <- data.frame(freq=spec$freq,spec=spec$spec)
    
    per_plot <- ggplot(df,aes(freq,spec)) +
      geom_line(color="blue") +
      ggtitle(paste(dataset_name,"-",sname,"Periodogram")) +
      xlab("Frequency") +
      ylab("Spectral Density") +
      theme_minimal()
    
    plot_list[[length(plot_list)+1]] <- acf_plot
    plot_list[[length(plot_list)+1]] <- per_plot
  }
  
  grid.arrange(grobs=plot_list,nrow=3,ncol=2)
}


###############################################################
# Load Datasets
###############################################################

datasets <- list()

data(AirPassengers)
datasets[["AirPassengers"]] <- AirPassengers

datasets[["NottinghamTemperatures"]] <- nottem

getSymbols("INDPRO",src="FRED")
datasets[["IndustrialProductionIndex"]] <-
  ts(as.numeric(INDPRO),start=c(1919,1),frequency=12)

data(soi)
soi_vec <- as.numeric(t(soi))
soi_vec <- soi_vec[!is.na(soi_vec)]

datasets[["SouthernOscillationIndex"]] <-
  ts(soi_vec,start=c(1950,1),frequency=12)

data("sunspot.month")
datasets[["Sunspot"]] <- sunspot.month

datasets[["SimulatedGegenbauer"]] <-
  ts(simulate_gegenbauer(),start=c(1950,1),frequency=12)


###############################################################
# Transformations
###############################################################

transformations <- list(
  
  "Raw" = list(log=FALSE,diff=0),
  "Log" = list(log=TRUE,diff=0),
  "First Diff" = list(log=FALSE,diff=1),
  "Log Diff" = list(log=TRUE,diff=1),
  "12 Mth Diff" = list(log=FALSE,diff=12),
  "Log 12 Mth Diff" = list(log=TRUE,diff=12)
  
)


###############################################################
# Main Loop
###############################################################

results <- data.frame()

i <- 1

for(ds_name in names(datasets)){
  
  y <- datasets[[ds_name]]
  
  # Plot diagnostics
  plot_series_acf_periodogram(y,ds_name)
  
  for(tr_name in names(transformations)){
    
    tr <- transformations[[tr_name]]
    
    est <- estimate_gegenbauer(
      y,
      log_transform = tr$log,
      diff_order = tr$diff
    )
    
    fc1 <- if(!is.null(est$fit)){
      tryCatch(forecast(est$fit,h=1)$mean[1],
               error=function(e) NA)
    } else NA
    
    row <- data.frame(
      Example = paste0("Example ",i),
      Dataset = ds_name,
      Transformation = tr_name,
      u = est$u,
      d = est$d,
      AR1 = est$AR,
      MA1 = est$MA,
      Sigma2 = est$Sigma2,
      Forecast_t1 = round(fc1,4)
    )
    
    results <- rbind(results,row)
    
    i <- i + 1
  }
}


###############################################################
# Final Results Table
###############################################################

print(results,row.names=FALSE)