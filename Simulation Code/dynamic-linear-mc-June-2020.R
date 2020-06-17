####################################################################
# MONTE CARLO CALIBRATED TO DYNAMIC LINEAR PANEL MODEL
####################################################################
# Cross-Over Jackknife of Nonlinear Panel
# Shuowen Chen, Victor Chernozhukov and Ivan Fernandez-Val
# Data source: Daron Acemoglu (MIT), N = 147 countries, T = 23 years (1987-2009)
# include data for four lags of dependent variable, balanced panel
#
# Description of the data: the sample selection and variable contruction follow ANRR
#
# The variables in the data set include:
#
# country_name    = Country name
# wbcode          = World Bank country code
# year            = Year 
# id              = Generated numeric country code
# dem             = Democracy measure by ANRR
# lgdp            = log of GDP per capita in 2000 USD from World Bank
####################################################################

library(foreign)
library(xtable)
library(plm)
library(gmm)
library(boot)
filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
setwd(filepath)

data2 <- read.dta("democracy-balanced-l4.dta")
data2 <- pdata.frame(data2, index = c("id", "year"))
N <- length(unique(data2$id))
T <- length(unique(data2$year))
T2 <- T - 4
################################
## Calibration by Fixed Effects
################################
# Regression Specification
form <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id) 
# Create lagged dependent variable
data2$l1lgdp <- lag(data2$lgdp, 1)
data2$l2lgdp <- lag(data2$lgdp, 2) 
data2$l3lgdp <- lag(data2$lgdp, 3) 
data2$l4lgdp <- lag(data2$lgdp, 4) 

##### 1. FE using data at hand ####
fit <- lm(form, data2, x = TRUE)
coefs <- coef(fit)[2:6]
sigma <- summary(fit)$sigma

##### 2. Compute X*beta_hat + time FE + ind FE (will be used in the data simulation) ####
# Drop the first four years at they contain NA due to lagging
speriod <- !((data2$year == 1987) | (data2$year == 1988) | (data2$year == 1989) | (data2$year == 1990))
sdata <- data2[speriod, ]
sdata$l1lgdp <- 0
sdata$l2lgdp <- 0
sdata$l3lgdp <- 0
sdata$l4lgdp <- 0
index0 <- predict(fit, sdata)
#index0 <- matrix(index0, nrow = T - 4, ncol = N)
rm(sdata)

## Auxiliary Functions: will be called in simulations for bias corrections
# Analytical Bias correction
abc <- function(data, formu, lags, N) {
  fit <- lm(formu, data, x = TRUE, na.action = na.omit)
  res <- fit$residuals
  jac <- solve(t(fit$x) %*% fit$x / length(res))[2:6, 2:6]
  indexes <- c(1:length(res))
  bscore <- rep(0, 5)
  T <- length(res)/N
  for (i in 1:lags) {
    indexes <- indexes[-c(1 + c(0:(N - 1))*T)]
    lindexes <- indexes - i
    bscore <- bscore + t(fit$x[indexes, 2:6]) %*% res[lindexes] / length(indexes)
  }
  bias <- -jac %*% bscore
  return(as.vector(bias/T))
}

abc2 <- function(data, formu, lags, N) {
  fit <- lm(formu, data, x = TRUE, na.action = na.omit)
  res <- fit$residuals
  jac <- solve(t(fit$x) %*% fit$x / length(res))
  indexes <- c(1:length(res))
  bscore <- rep(0, 5)
  T <- length(res)/N
  for (i in 1:lags) {
    indexes   <- indexes[-c(1 + c(0:(N - 1))*T)]
    lindexes  <- indexes - i
    bscore  <- bscore + t(fit$x[indexes, ]) %*% res[lindexes] / length(res)
  }
  bias <- -jac %*% bscore
  return(bias[2:6]/T)
}

# Split Sample Jackknife (F-V and Weidner, 2016)
sbc <- function(data, formu) {
  N <- length(unique(data$id))
  T <- length(unique(data$year)) - 4
  unit <- as.double(data$id)
  time <- as.double(data$year)
  # Full sample estimation
  whole <- plm(formu, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Split across time series dimension
  t_lower <- plm(formu, data = data[time <= ceiling(T/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  t_upper <- plm(formu, data = data[time >= floor(1 + T/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  bc_t <- 0.5*(coef(t_lower) + coef(t_upper))
  # Split across individual dimension
  n_lower <- plm(formu, data = data[unit <= ceiling(N/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  n_upper <- plm(formu, data = data[unit >= floor(1 + N/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  bc_n <- 0.5*(coef(n_lower) + coef(n_upper))
  # SBC bias correction
  coef_sbc <- 3*coef(whole) - bc_t - bc_n
  return(coef_sbc)
}

# Cross-over Bias Correction (New Proposal)
cbc <- function(data, formu) {
  unit <- as.double(data$id)
  time <- as.double(data$year)
  # Full sample estimation
  whole <- plm(formu, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Two cross-over subsamples
  subsample1 <- (unit <= (median(unit)) & time <= median(time)) | 
    (unit >= median(unit) & time >= median(time))
  subsample2 <- (unit <= median(unit) & time >= median(time)) | 
    (unit >= median(unit) & time <= median(time))
  cross1 <- plm(formu, data = data[subsample1, ], model = "within", effect = "twoways", index = c("id", "year"))
  cross2 <- plm(formu, data = data[subsample2, ], model = "within", effect = "twoways", index = c("id", "year"))
  coef_cbc <- 2*coef(whole) - 0.5*(coef(cross1) + coef(cross2))
  return(coef_cbc)
}

#### 3. Simulation ####
set.seed(88)
R <- 500
y <- matrix(0, nrow = T2, ncol = N)
l1y <- matrix(0, nrow = T2, ncol = N)
l2y <- matrix(0, nrow = T2, ncol = N)
l3y <- matrix(0, nrow = T2, ncol = N)
l4y <- matrix(0, nrow = T2, ncol = N)
# first non NA entry of l1lgdp for each country 
l1y[1, ] <- data2$l1lgdp[2 + c(0:(N - 1))*T]
l2y[1, ] <- data2$l2lgdp[3 + c(0:(N - 1))*T]
l3y[1, ] <- data2$l3lgdp[4 + c(0:(N - 1))*T]
l4y[1, ] <- data2$l4lgdp[5 + c(0:(N - 1))*T]
# placeholders
coefs_fe_b <- matrix(0, nrow = R, ncol = 5)
coefs_abc4 <- matrix(0, nrow = R, ncol = 5)
coefs_abc8 <- matrix(0, nrow = R, ncol = 5)
coefs_abc12 <- matrix(0, nrow = R, ncol = 5)
coefs_sbc <- matrix(0, nrow = R, ncol = 5)
coefs_cbc <- matrix(0, nrow = R, ncol = 5)
se_fe <- matrix(0, nrow = R, ncol = 5)
# loop over simulations
for (i in 1:R) {
  y[1, ] <- index0[1 + c(0:(N - 1))*T2] + coefs[2]*l1y[1, ] + coefs[3]*l2y[1, ] + 
    coefs[4]*l3y[1, ] + coefs[5]*l4y[1, ] + rnorm(N, mean = 0, sd = sigma)
  for (t in 2:T2) {
    # generate lagged variables
    l4y[t, ] <- l3y[t - 1, ]
    l3y[t, ] <- l2y[t - 1, ]
    l2y[t, ] <- l1y[t - 1, ]
    l1y[t, ] <- y[t - 1, ]
    # update dependent variable at time t
    y[t, ] <- index0[t + c(0:(N - 1))*T2] + coefs[2]*l1y[t, ] + coefs[3]*l2y[t, ] + 
      coefs[4]*l3y[t, ] + coefs[5]*l4y[t, ] + rnorm(N, mean = 0, sd = sigma)
  }
  # Construct the dataset
  ys <- matrix(y, ncol = 1)
  l1ys <- matrix(l1y, ncol = 1)
  l2ys <- matrix(l2y, ncol = 1)
  l3ys <- matrix(l3y, ncol = 1)
  l4ys <- matrix(l4y, ncol = 1)
  data_b <- data.frame(id = kronecker(c(1:N), rep(1, T2)), year = kronecker(rep(1, N), c(1:T2)), 
                       lgdp = ys, l1lgdp = l1ys, l2lgdp = l2ys, l3lgdp = l3ys, l4lgdp = l4ys, 
                       dem = fit$x[, "dem"])
  #slgdp <- matrix(data2$lgdp, nrow = T, ncol = N)
  # Generate simulated dependent variable from t=5 
  #for (t in 5:T) {
  #  slgdp[t, ] <- index0[t - 4, ] + coefs[2]*slgdp[t - 1, ] + coefs[3]*slgdp[t - 2, ] + coefs[4]*slgdp[t - 3, ] + coefs[5]*slgdp[t - 4, ] + sigma*rnorm(N)
  #}
  #data_b <- data2
  # Replace with simulated dependent variable
  #data_b$lgdp <- matrix(slgdp, ncol = 1)
  #data_b <- data.frame(data_b)
  #data_b <- pdata.frame(data_b, index = c("id", "year"))
  
  # Fixed Effects
  fit_b <- lm(form, data_b)
  coefs_fe_b[i, ] <- coef(fit_b)[2:6]
  # abc
  bias_l4 <- abc(data =  data_b, formu = form, lags = 4, N = length(unique(data_b$id)))
  bias_l8 <- abc(data =  data_b, formu = form, lags = 8, N = length(unique(data_b$id)))
  bias_l12 <- abc(data =  data_b, formu = form, lags = 12, N = length(unique(data_b$id)))
  coefs_abc4[i, ] <- coefs_fe_b[i, ] - bias_l4 
  coefs_abc8[i, ] <- coefs_fe_b[i, ] - bias_l8
  coefs_abc12[i, ] <- coefs_fe_b[i, ] - bias_l12
  # sbc
  coefs_sbc[i, ] <- sbc(data_b, formu = form)
  # cbc
  coefs_cbc[i, ] <- cbc(data_b, formu = form)
  # standard errors
  se_fe[i, ] <- coef(summary(fit))[2:6, 2]
}

# Post Simulation Summary

# A function that calculates bias, std, RMSE, SE/SD and pvalue at 0.95
sumstatslinear <- function(sim, real, se_sim) {
  output <- matrix(0, nrow = dim(sim)[2], ncol = 5)
  colnames(output) <- c("bias", "std", "RMSE", "SE/SD", "p95")
  rownames(output) <- names(real)
  for (i in 1:dim(sim)[2]) {
    # bias
    output[i, 1] <- 100*(mean(sim[, i], na.rm = TRUE)/real[i] - 1)
    # standard deviation
    output[i, 2] <- 100*sd(sim[, i]/real[i], na.rm = TRUE)
    # rmse
    output[i, 3] <- 100*sqrt(mean((sim[, i]/real[i] - 1)^2, na.rm = TRUE))
    # se_sd
    output[i, 4] <- (mean((se_sim[, i]/sd(sim[, i], na.rm = TRUE)), na.rm = TRUE))
    # pvalue
    output[i, 5] <- (mean((sim[, i] + qnorm(.05/2) * se_sim[, i] <= real[i]) &
                            (sim[, i] + qnorm(1 - .05/2) * se_sim[, i] >= real[i]), na.rm = TRUE))
  }
  return(output)
}

sum_nobc <- sumstatslinear(coefs_fe_b, coefs, se_fe)
sum_abc4 <- sumstatslinear(coefs_abc4, coefs, se_fe)
sum_abc8 <- sumstatslinear(coefs_abc8, coefs, se_fe)
sum_abc12 <- sumstatslinear(coefs_abc12, coefs, se_fe)
sum_sbc <- sumstatslinear(coefs_sbc, coefs, se_fe)
sum_cbc <- sumstatslinear(coefs_cbc, coefs, se_fe)
