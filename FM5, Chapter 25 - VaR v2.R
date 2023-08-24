# We thank Sagi Haim for developing this script

##############################
# CHAPTER 25 - VALUE AT RISK
##############################

### 25.1 Overview
## A Simple Example
# Inputs:
mean <- 0.20
sigma <- 0.30
init_invst <- 100
cutoff <- 80

# (Cumulative) Probability that portfolio worth less than cutoff
pnorm(cutoff, mean = init_invst * (1 + mean),
      sd = sigma * init_invst, log = FALSE) # similar to Excel's NORM.DIST

# Probability distribution
port_value <- c(0:25)*10  #bins
prob <- sapply(port_value, function(x) pnorm(x, mean = init_invst * (1 + mean),
                                     sd = sigma * init_invst, log = FALSE))

# Plot (comulative distribution)
plot(port_value, prob, type = "l", main = "End-of-Year Portfolio Value", 
     xlab = "Portfolio Value", ylab = "Cumulative Probability")
lines(port_value,prob, type = "p")

### 25.2  The three types of VaR models
##  Analytic VaR 
# in this example: 10M (weekly)			
# Inputs: 
mean_annual <- 0.02
mean_weekly	<- mean_annual*1/52
sigma_annual	<- 0.10
sigma_week	<- sigma_annual * sqrt(1/52)
port_size 	<- 10000000
spot_rate	<- 1.15 #(Dollar/Euro spot rate)

# Standard deviations
prob <- c(0.9, 0.95, 0.99)
st_dev <- qnorm(prob, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) 
names(st_dev) <- prob

# VaR (1 week, $)
var_analytic <- st_dev * sigma_week * port_size * spot_rate
names(var_analytic) <- prob
var_analytic

## Calculating the quantiles
mean <- 0.20
sigma <- 0.30
init_invst <- 100
cutoff <- qnorm(0.01, mean = init_invst * (1 + mean),
      sd = sigma * init_invst, lower.tail = TRUE, log.p = FALSE)

# Var at 1.00% Level
var_analytic <- init_invst - cutoff

## Quantiles for Log-Normal distribution
init_invst <- 100
mean <- 0.10
sigma <- 0.30
t <- 1.0

#Parameters of normal distribution of ln(VT)
mean <- log(init_invst) + (mean - 0.5*sigma^2) * t
sigma <- sigma * sqrt(t)
cutoff <- qlnorm(0.01, mean = mean, sd = sigma)

# Var at 1.00% Level
var_analytic <- init_invst - cutoff

##  Monte-Carlo VaR of 10M (weekly)
# Simulate weekly returns
n <- 10000 #number of simulations
Z <- rnorm(n, mean = 0, sd = 1) 
week_ret <- (mean_weekly - 0.5 * sigma_week ^ 2) + sigma_week * Z

# Calculate Percentiles
percentile <- quantile(week_ret, probs = 1-prob) # same as PERCENTILE.INC in Excel

# Calculate Value at Risk
Var_monte_carlo <- percentile * port_size * spot_rate

##  Historical VaR (weekly)
# Set working directory
workdir <- readline(prompt="working directory?")
setwd(workdir)

# Read data
data <- read.csv("Chapter_25_Data_1.csv", row.names = 1)

# Sort the data by date
data <- data[order(as.Date(row.names(data), format="%d-%m-%y"),decreasing = FALSE),]

# Calculate returns
week_ret <- diff(log(data))
port_change <- port_size * exp(week_ret) - port_size

# Calculate Value at Risk
Var_historical <- quantile(port_change, probs = 1 - prob) # same as PERCENTILE.INC in Excel

### 25.3  VaR of n-asset portfolio
## Analytic VaR of 3-asset portfolio
# Inputs:
mean_ret <- c(0.1, 0.12, 0.13)
cov_mat <- matrix(c(0.10,	0.04,	0.03,
                    0.04,	0.20,	-0.04,
                    0.03,	-0.04,	0.60), 
                  nrow = 3)
port_prop <- c(0.3, 0.25, 0.45)
Init_invst	<- 100

# Portfolio Stats
port_mean_ret <- mean_ret %*% port_prop
port_sigma	<- sqrt(port_prop %*% cov_mat %*% port_prop)
mean_invs_value <- Init_invst * ( 1 + port_mean_ret )
invs_value_sigma <- port_sigma * Init_invst

# Analytic VaR
sig_level <- 0.01
cutoff <- qnorm(p = sig_level, mean = mean_invs_value, sd = invs_value_sigma)

# VaR at 1.00%  level
var_analytic <- cutoff - Init_invst 

## Monte-Carlo simulation VaR of the 3-asset portfolio

# Assets Variance
asset_var <- diag(cov_mat)

# Lower-Triangular Cholesky Decomposition 
chol_decompos <- t(chol(cov_mat)) 

n <- 10000 # number of simulations

# Simulate correlated random variables
Z <- t( sapply(1:n, function(x) chol_decompos %*% rnorm(n = 3, mean = 0, sd = 1)))

# Simulate Monthly Returns Of 3 Correlated Stocks
sim_asset_ret <- apply(Z, 1, function(z) ( (mean_ret - 0.5 * asset_var) 
                                           + sqrt(asset_var) * z) )

# Simulate Portfolio Return
port_ret <- t(sim_asset_ret) %*% port_prop

# Calculate the 1% Percentile Return 
percentile <- quantile(port_ret, probs = 0.01)

# VaR at 1.00%  level
var_monte_carlo <- percentile * Init_invst 

## Historical VaR of the 3-asset portfolio
# Read the return data from the csv file
asset_ret <- as.matrix(read.csv("Chapter_25_Data_2.csv"))
port_ret <- asset_ret %*% port_prop

# Calculate the 1% Percentile Return 
percentile <- quantile(port_ret, probs = 0.01)

# VaR at 1.00%  level
Var_historical <- percentile * Init_invst 

### 25.4  Backtesting
## The Standard Coverage Test
# Inputs
VaR_level	<- 0.01
n	<- 750
Signif_level <- 0.05

# Confidence interval
upper_bound <- qbinom(1-Signif_level/2, n, VaR_level) # Similar to Excel's BINOM.INV
lower_bound <- n - qbinom(1-Signif_level/2, n, 1-VaR_level) 
non_rejection_values <-  paste(lower_bound, "-", upper_bound, sep = "")

# Data Table variables
VaR_levels <- c(0.005, 0.01, 0.025, 0.05, 0.1)
n <- c(50, 100, 200, 250, 500, 750, 1000, 1250, 1500) #Sample size

#  define function for the non-rjection values
non_rej_val <- function(VaR_level, n){
  upper_bound <- qbinom(1-Signif_level/2, n, VaR_level) # Similar to Excel's BINOM.INV
  lower_bound <- n - qbinom(1-Signif_level/2, n, 1-VaR_level) 
  non_rejection_values <-  paste(lower_bound, "-", upper_bound, sep = "")
  return(non_rejection_values)
}

# run function for each combination
results <- outer(VaR_levels, n, FUN = "non_rej_val")
colnames(results) <- n
rownames(results) <- VaR_levels
results
