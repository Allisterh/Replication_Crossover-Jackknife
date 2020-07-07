# 6/19/2020 Shuowen Chen
# Functions for bias correction in dynamic panel
# Note: the current functions are tailed to the dynamic linear panel with four lags
# Later on will revise the code so as to be more general

# 1. Analytical bias correction
abc <- function(data, formu, lags, N) {
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

# 2. Split-sample bias correction
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

# 3. Cross-over bias correction
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
