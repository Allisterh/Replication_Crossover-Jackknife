library(foreign)
library(xtable)
library(plm)
library(boot)

filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
setwd(filepath)
###  read-in Democracy data
democracydata <- read.dta("democracy-balanced-l4.dta")
democracydata <- pdata.frame(democracydata, index = c("id", "year"))
# create lagged terms
democracydata$l1lgdp <- lag(democracydata$lgdp, 1)
democracydata$l2lgdp <- lag(democracydata$lgdp, 2) 
democracydata$l3lgdp <- lag(democracydata$lgdp, 3) 
democracydata$l4lgdp <- lag(democracydata$lgdp, 4) 

# Bias correction functions
# 1. Analytical Bias correction
abc <- function(data, form, lags, N) {
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
  bias <- -jac %*% bscore;
  return(as.vector(bias/T))
}

# a function that randomly splits s times along the cross-section
# and pool to get mega half panels for crossover jackknife bc
mega_data <- function(s, n, unit, time, dat, uc, fm) {
  mega1 <- list()
  mega2 <- list()
  for (i in 1:s) {
    index <- sample(n, ceiling(n/2), replace = FALSE)
    sub1 <- (unit %in% index & time <= median(time)) |
      (!(unit %in% index) & time >= median(time))
    sub2 <- (unit %in% index & time >= median(time)) | 
      (!(unit %in% index) & time <= median(time))
    mega1[[i]] <- dat[sub1, ]
    mega2[[i]] <- dat[sub2, ]
  }
  mega1 <- do.call(rbind, mega1)
  mega2 <- do.call(rbind, mega2)
  return(list(mega1, mega2))
}

# 2. Split Sample Jackknife (F-V and Weidner, 2016)
sbc <- function(data, formula, lag_dependent) {
  N <- length(unique(data$id))
  T <- length(unique(data$year)) - lag_dependent
  unit <- as.double(data$id)
  time <- as.double(data$year)
  # Full sample estimation
  whole <- plm(formula, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Split across time series dimension
  t_lower <- plm(formula, data = data[time <= ceiling(T/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  t_upper <- plm(formula, data = data[time >= floor(1 + T/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  bc_t <- 0.5*(coef(t_lower) + coef(t_upper))
  # Split across individual dimension
  n_lower <- plm(formula, data = data[unit <= ceiling(N/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  n_upper <- plm(formula, data = data[unit >= floor(1 + N/2), ], model = "within", 
                 effect = "twoways", index = c("id", "year"))
  bc_n <- 0.5*(coef(n_lower) + coef(n_upper))
  # SBC bias correction
  coef_sbc <- 3*coef(whole) - bc_t - bc_n
  return(coef_sbc)
}

# 3. Cross-over Bias Correction (New Proposal)
cbc <- function(data, formula) {
  unit <- as.double(data$id)
  time <- as.double(data$year)
  # Full sample estimation
  whole <- plm(formula, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Two cross-over subsamples
  subsample1 <- (unit <= (median(unit)) & time <= median(time)) | 
    (unit >= median(unit) & time >= median(time))
  subsample2 <- (unit <= median(unit) & time >= median(time)) | 
    (unit >= median(unit) & time <= median(time))
  cross1 <- plm(formula, data = data[subsample1, ], model = "within", effect = "twoways", index = c("id", "year"))
  cross2 <- plm(formula, data = data[subsample2, ], model = "within", effect = "twoways", index = c("id", "year"))
  coef_cbc <- 2*coef(whole) - 0.5*(coef(cross1) + coef(cross2))
  return(coef_cbc)
}


########## Calibration
##### 1. FE using data at hand ####
form <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id) 
fit <- lm(form, democracydata, x = TRUE)
coefs <- coef(fit)[2:6]
sigma <- summary(fit)$sigma

##### 2. Compute X*beta_hat + time FE + ind FE (will be used in the data simulation) ####
# Drop the first four years at they contain NA due to lagging
speriod <- !((democracydata$year == 1987) | (democracydata$year == 1988) | (democracydata$year == 1989) | (democracydata$year == 1990))
sdata <- democracydata[speriod, ]
sdata$l1lgdp <- 0
sdata$l2lgdp <- 0
sdata$l3lgdp <- 0
sdata$l4lgdp <- 0
index0 <- predict(fit, sdata)
#index0 <- matrix(index0, nrow = T - 4, ncol = N)
rm(sdata)

########## Simulation
# Simulate data function
data_sim <- function(data, index = index0, coeff = coefs, sig = sigma, lmfit = fit) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  T2 <- T - 4
  y <- matrix(0, nrow = T2, ncol = N)
  l1y <- matrix(0, nrow = T2, ncol = N)
  l2y <- matrix(0, nrow = T2, ncol = N)
  l3y <- matrix(0, nrow = T2, ncol = N)
  l4y <- matrix(0, nrow = T2, ncol = N)
  # first NA entry of l1lgdp for each country s
  l1y[1, ] <- data$l1lgdp[2 + c(0:(N - 1))*T]
  l2y[1, ] <- data$l2lgdp[3 + c(0:(N - 1))*T]
  l3y[1, ] <- data$l3lgdp[4 + c(0:(N - 1))*T]
  l4y[1, ] <- data$l4lgdp[5 + c(0:(N - 1))*T]
  # First observation for each country
  y[1, ] <- index[1 + c(0:(N - 1))*T2] + coeff[2]*l1y[1, ] + coeff[3]*l2y[1, ] + 
    coeff[4]*l3y[1, ] + coeff[5]*l4y[1, ] + rnorm(N, mean = 0, sd = sig)
  # generate lagged variables
  for (t in 2:T2) {
    l4y[t, ] <- l3y[t - 1, ]
    l3y[t, ] <- l2y[t - 1, ]
    l2y[t, ] <- l1y[t - 1, ]
    l1y[t, ] <- y[t - 1, ]
    # update dependent variable at time t
    y[t, ] <- index[t + c(0:(N - 1))*T2] + coeff[2]*l1y[t, ] + coeff[3]*l2y[t, ] + 
      coeff[4]*l3y[t, ] + coeff[5]*l4y[t, ] + rnorm(N, mean = 0, sd = sig)
  }
  # Construct the dataset
  ys <- matrix(y, ncol = 1)
  l1ys <- matrix(l1y, ncol = 1)
  l2ys <- matrix(l2y, ncol = 1)
  l3ys <- matrix(l3y, ncol = 1)
  l4ys <- matrix(l4y, ncol = 1)
  #data_b <- data.frame(id = kronecker(c(1:N), rep(1, T2)), year = kronecker(rep(1, N), c(1:T2)), 
  #                     lgdp = ys, l1lgdp = l1ys, l2lgdp = l2ys, l3lgdp = l3ys, l4lgdp = l4ys, 
  #                     dem = lmfit$x[, "dem"])
  return(data)
}

# bootstrap function
## statistics to be computed in each bootstrap draw ##
check <- function(data, fm) {
  return(data)
}
boot_fe <- function(data, fm) {
  # Fixed Effects
  fit_b <- lm(fm, data) 
  coefs_fe_b <- coef(fit_b)[2:6]
  # abc
  bias4 <- abc(data, form = fm, lags = 4, N = length(unique(data$id)))
  bias8 <- abc(data, form = fm, lags = 8, N = length(unique(data$id)))
  bias12 <- abc(data, form = fm, lags = 12, N = length(unique(data$id)))
  coefs_abc4_b <- coefs_fe_b - bias4
  coefs_abc8_b <- coefs_fe_b - bias8
  coefs_abc12_b <- coefs_fe_b - bias12
  # sbc
  coefs_sbc_b <- sbc(data, formula = fm, lag_dependent = 4)
  # cbc
  coefs_cbc_b <- cbc(data, formula = fm)
  # standard errors
  se_fe_b <- coef(summary(fit_b))[2:6, 2]
  # Return outputs
  return(c(coefs_fe_b, coefs_abc4_b, coefs_abc8_b, coefs_abc12_b, coefs_sbc_b, coefs_cbc_b, se_fe_b))
}

set.seed(88)
R <- 5
# use the boot library
fe_sim <- boot(data = democracydata, statistic = check, sim = "parametric", ran.gen = data_sim, mle = 0, 
               fm = form, parallel = "multicore", ncpus = 1, R = 2)







