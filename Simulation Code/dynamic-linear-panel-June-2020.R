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
data_rg <- function(data, coef = coefs, sig = sigma, index = index0, mle) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  slgdp <- matrix(data$lgdp, nrow = T, ncol = N)
  # Generate simulated dependent variable from t=5 
  for (t in 5:T) {
    slgdp[t, ] <- index[t - 4, ] + coef[2]*slgdp[t - 1, ] + coef[3]*slgdp[t - 2, ] + coef[4]*slgdp[t - 3, ] + coef[5]*slgdp[t - 4, ] + sig*rnorm(N)
  }
  data_b <- data
  # Replace with simulated dependent variable
  data_b$lgdp <- matrix(slgdp, ncol = 1)
  data_b <- data.frame(data_b)
  data_b <- pdata.frame(data_b, index = c("id", "year"))  
  return(data_b)
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

## statistics to be computed in each bootstrap draw ##
boot_fe <- function(data, form_fe, form_abc) {
  # Fixed Effects
  fit <- lm(form_fe, data) # essentially the same as plm, but using richer summary output
  #fe_fit <- plm(form_fe, data, model = "within", effect = "twoways", index = c("id", "year"))
  coefs_fe <- coef(fit)[1:5]
  # abc
  bias_l4 <- abc(data, form = form_abc, lags = 4, N = length(levels(data$id)))
  coefs_abc4 <- coefs_fe - bias_l4 
  # sbc
  coefs_sbc <- sbc(data, form = form_fe)
  # cbc
  coefs_cbc <- cbc(data, form = form_fe)
  # standard errors
  se_fe <- coef(summary(fit))[2:6, 2]
  # Return outputs
  return(c(coefs_fe, coefs_abc4, coefs_sbc, coefs_cbc, se_fe))
}

set.seed(88)
R <- 500
coefs_fe_b <- matrix(0, nrow = R, ncol = 5)
coefs_abc4 <- matrix(0, nrow = R, ncol = 5)
coefs_sbc <- matrix(0, nrow = R, ncol = 5)
coefs_cbc <- matrix(0, nrow = R, ncol = 5)
se_fe <- matrix(0, nrow = R, ncol = 5)
for (i in 1:R) {
  slgdp <- matrix(data$lgdp, nrow = T, ncol = N)
  # Generate simulated dependent variable from t=5 
  for (t in 5:T) {
    slgdp[t, ] <- index[t - 4, ] + coefs[2]*slgdp[t - 1, ] + coefs[3]*slgdp[t - 2, ] + coefs[4]*slgdp[t - 3, ] + coefs[5]*slgdp[t - 4, ] + sigma*rnorm(N)
  }
  data_b <- data
  # Replace with simulated dependent variable
  data_b$lgdp <- matrix(slgdp, ncol = 1)
  data_b <- data.frame(data_b)
  data_b <- pdata.frame(data_b, index = c("id", "year"))
  # Fixed Effects
  fit_b <- lm(form, data_b) # essentially the same as plm, but using richer summary output
  #fe_fit <- plm(form_fe, data, model = "within", effect = "twoways", index = c("id", "year"))
  coefs_fe_b[i, ] <- coef(fit_b)[2:6]
  # abc
  bias_l4 <- abc(data_b, form, lags = 4, N = length(levels(data$id)))
  coefs_abc4[i, ] <- coefs_fe_b[i, ] - bias_l4 
  # sbc
  coefs_sbc[i, ] <- sbc(data_b, form)
  # cbc
  coefs_cbc[i, ] <- cbc(data_b, form)
  # standard errors
  se_fe[i, ] <- coef(summary(fit))[2:6, 2]
}


#form_fe <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id)
#form_abc <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id)
#sim_fe <- boot(data = data2, statistic = boot_fe, sim = "parametric", ran.gen = data_rg, 
#               mle = 0, form_fe = form_fe, form_abc = form_abc,  parallel = "multicore", 
#               ncpus = 1, R = 5)

# Post Simulation Summary

#sim_nobc <- sim_fe$t[, 1:5]
#sim_abc4 <- sim_fe$t[, 6:10]
#sim_sbc <- sim_fe$t[, 11:15]
#sim_cbc <- sim_fe$t[, 16:20]
#sim_se <- sim_fe$t[, 21:25]

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
sum_sbc <- sumstatslinear(coefs_sbc, coefs, se_fe)
sum_cbc <- sumstatslinear(coefs_cbc, coefs, se_fe)

# robust estimator of std deviation based on IQR
rsd <- function(x) {
  return((quantile(x, .75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25)))
}