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

data2 <- read.dta("democracy-balanced-l4.dta")
data2 <- pdata.frame(data2, index = c("id", "year"))
N <- length(unique(data2$id))
T <- length(unique(data2$year))
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
fit <- lm(form, data2) 
coefs <- coef(fit)[2:6]
sigma <- summary(fit)$sigma
# Long run effects of treatment
lr <- coefs[1]/(1 - sum(coefs[2:5]))

##### 2. Compute X*beta_hat + time FE + ind FE (will be used in the data simulation) ####
# Drop the first four years at they contain NA due to lagging
speriod <- !((data2$year == 1987) | (data2$year == 1988) | (data2$year == 1989) | (data2$year == 1990))
sdata <- data2[speriod, ]
sdata$l1lgdp <- 0
sdata$l2lgdp <- 0
sdata$l3lgdp <- 0
sdata$l4lgdp <- 0
index0 <- predict(fit, sdata)
index0 <- matrix(index0, nrow = T - 4, ncol = N)
rm(sdata)

## Auxiliary Functions 
# (a) Generate data sets
data.rg <- function(data, coef = coefs, sig = sigma, index = index0, mle) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  coef <- coef[-1]
  slgdp <- matrix(data$lgdp, nrow = T, ncol = N)
  # Generate simulated dependent variable from t=5 
  for (t in 5) {
    slgdp[t, ] <- index[t - 4, ] + coefs[1]*slgdp[t - 1, ] + coefs[2]*slgdp[t - 2, ] + coefs[3]*slgdp[t - 3, ] + coefs[4]*slgdp[t - 4, ] + sig*rnorm(N)
  }
  data_b <- data
  # Replace with simulated dependent variable
  data_b$lgdp <- matrix(slgdp, ncol = 1)
  data_b <- data.frame(data_b)
  data_b <- pdata.frame(data_b, index = c("id", "year"))  
  return(data_b)
}

data_rg <- function(data, fit_original = fit, index0 = index0, mle) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  coefs <- coef(fit_original)[3:6] # coefs on lagged dependent variable
  sigma <- summary(fit_original)$sigma 
  slgdp <- matrix(data$lgdp, nrow = T, ncol = N)
  index0 <- matrix(index0, nrow = T - 4, ncol = N)
  # Generate simulated dependent variable from t=5 
  for (t in 5:T) {
    slgdp[t, ] <- index0[t - 4, ] + coefs[1]*slgdp[t - 1, ] + 
      coefs[2]*slgdp[t - 2, ] + coefs[3]*slgdp[t - 3, ] + 
      coefs[4]*slgdp[t - 4, ] + sigma*rnorm(N)
  }
  data_b <- data
  # Replace with simulated dependent variable
  data_b$lgdp <- matrix(slgdp, ncol = 1)
  data_b <- data.frame(data_b)
  data_b <- pdata.frame(data_b, index = c("id", "year"))    
  return(data_b)
}

check <- function(data, form_fe, form_abc) {
  return(names(data))
}
# (b) Analytical Bias correction
abc <- function(data, form, lags, N) {
  data$l1lgdp <- lag(data$lgdp, 1)
  data$l2lgdp <- lag(data$lgdp, 2) 
  data$l3lgdp <- lag(data$lgdp, 3) 
  data$l4lgdp <- lag(data$lgdp, 4) 
  fit <- lm(form, data, x = TRUE, na.action = na.omit)
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

# (c) Split Sample Jackknife (F-V and Weidner, 2016)
sbc <- function(data, form) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  unit <- as.double(data$id)
  time <- as.double(data$year)
  # Full sample estimation
  whole <- plm(form, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Split across time series dimension
  t_lower <- plm(form, data = data[time <= ceiling(T/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  t_upper <- plm(form, data = data[time >= floor(1 + T/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  bc_t <- 0.5*(coef(t_lower) + coef(t_upper))
  # Split across individual dimension
  n_lower <- plm(form, data = data[unit <= ceiling(N/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  n_upper <- plm(form, data = data[unit >= floor(1 + N/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  bc_n <- 0.5*(coef(n_lower) + coef(n_upper))
  # SBC bias correction
  coef_sbc <- 3*coef(whole) - bc_t - bc_n
  return(coef_sbc)
}

# (d) Cross-over Bias Correction (New Proposal)
cbc <- function(data, form) {
  unit <- as.double(data$id)
  time <- as.double(data$year)
  # Full sample estimation
  whole <- plm(form, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Two cross-over subsamples
  subsample1 <- (unit <= (median(unit)) & time <= median(time)) | 
    (unit >= median(unit) & time >= median(time))
  subsample2 <- (unit <= median(unit) & time >= median(time)) | 
    (unit >= median(unit) & time <= median(time))
  cross1 <- plm(form, data = data[subsample1, ], model = "within", effect = "twoways", index = c("id", "year"))
  cross2 <- plm(form, data = data[subsample2, ], model = "within", effect = "twoways", index = c("id", "year"))
  coef_cbc <- 2*coef(whole) - 0.5*(coef(cross1) + coef(cross2))
  return(coef_cbc)
}
#### 3. Simulation ####
R <- 500

## statistics to be computed in each bootstrap draw ##
boot_fe <- function(data, form_fe, form_abc) {
  # Fixed Effects
  fe_fit <- plm(form_fe, data, model = "within", effect = "twoways", index = c("id", "year"))
  coefs_fe <- coef(fe_fit)
  # abc
  bias_l4 <- abc(data, form = form_abc, lags = 4, N = length(levels(data$id)))
  coefs_abc4 <- coefs_fe - bias_l4 
  # sbc
  coefs_sbc <- sbc(data, form = form_fe)
  # cbc
  coefs_cbc <- cbc(data, form = form_fe)
  # Return outputs
  return(c(coefs_fe, coefs_abc4, coefs_sbc, coefs_cbc))
}

set.seed(888)
form_fe <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id) - 1
form_abc <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id)

sim_fe <- boot(data = data2, statistic = boot_fe, sim = "parametric", ran.gen = data.rg, 
               mle = 0, form_fe = form_fe, form_abc = form_abc,  parallel = "multicore", 
               ncpus = 1, R = 5)

# robust estimator of std deviation based on IQR
rsd <- function(x) {
  return((quantile(x, .75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25)))
}