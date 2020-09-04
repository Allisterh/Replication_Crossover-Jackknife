# The calibration exercise use the DGP that features (1) a linear time trend (2) a quadratic time trend, with 1 lagged dependent variable
democracydata <- read.dta("democracy-balanced-l4.dta")
datanew <- data.table(democracydata)
#lag <- datanew[, shift(lgdp, 1:4), by = id][, 2:5]
# create a lagged dependent variable
lag <- datanew[, shift(lgdp, 1), by = id][, 2]
colnames(lag) <- c("l1lgdp")
# merge with original data set
democracydata <- data.frame(democracydata, lag)
data2 <- pdata.frame(datanew, index = c("id", "year"))
data2$l1lgdp <- lag(data2$lgdp, 1)
# Note: democracydata is not set to be panel while data2 is. 

N <- length(unique(democracydata$id)) # number of countries
T <- length(unique(democracydata$year)) # total number of years
T2 <- T - 1 # T - number of lags
timeindex <- rep(seq(T), N)
demdata <- cbind(democracydata, timeindex)

reg1 <- lgdp ~ dem + l1lgdp + factor(id) + timeindex
reg1_quad <- lgdp ~ dem + l1lgdp + factor(id) + I(timeindex^2)
#reg2 <- lgdp ~ dem + l1lgdp + factor(id) + factor(id):year
#reg3 <- lgdp ~ dem + l1lgdp + factor(id)*I(year^2)
#trend2_lm <- lgdp ~ dem + l1lgdp + factor(id)*I(year^2)
#trend_plm <- lgdp ~ dem + l1lgdp + id : as.integer(year)
#trend2_plm <- lgdp ~ dem + l1lgdp + id : I(as.integer(year)^2)

# Regression specifications
# unit fixed effect with a time trend
# reg2 <- lgdp ~ dem + l1lgdp + factor(id)*year
# f <- lm(reg, democracydata)
# quadratic time trend
# reg3 <- lgdp ~ dem + l1lgdp + factor(id)*I(year^2)

# Estimation using true data
fit_trend <- lm(reg1, demdata, x = TRUE)
fit_trend2 <- lm(reg1_quad, demdata, x = TRUE)
sigma_trend <- summary(fit_trend)$sigma
sigma_trend2 <- summary(fit_trend2)$sigma
coef_trend <- fit_trend$coefficients[2:3]
coef_trend2 <- fit_trend2$coefficients[2:3]

# Drop the first year
speriod <- !(demdata$year == 1987)
sdata <- demdata[speriod, ]
sdata$l1lgdp <- 0
index_trend <- predict(fit_trend, sdata)
index_trend2 <- predict(fit_trend2, sdata)
rm(sdata)

set.seed(88) # seed
R <- 500 # number of simulations

# function for cross-over bias correction
cbc <- function(data, formu) {
  unit <- as.double(data$id)
  time <- as.double(data$year)
  # Full sample estimation
  whole <- lm(formu, data)
  # Two cross-over subsamples
  subsample1 <- (unit <= (median(unit)) & time <= median(time)) | 
    (unit >= median(unit) & time >= median(time))
  subsample2 <- (unit <= median(unit) & time >= median(time)) | 
    (unit >= median(unit) & time <= median(time))
  # subsample regression
  cross1 <- lm(formu, data = data[subsample1, ])
  cross2 <- lm(formu, data = data[subsample2, ])
  coef_cbc <- 2*coef(whole)[2:3] - 0.5*(coef(cross1)[2:3] + coef(cross2)[2:3])
  return(coef_cbc)
}

# function for sample splitting bias correction
sbc <- function(data, formu) {
  # only need to split over the time series dimension as there is no time fixed effects
  T <- length(unique(data$year)) - 1
  time <- as.double(data$year)
  # Full sample estimation
  whole <- lm(formu, data)
  # Split across time series dimension
  t_lower <- lm(formu, data = data[time <= ceiling(T/2), ])
  t_upper <- lm(formu, data = data[time >= floor(1 + T/2), ])
  bc_t <- 0.5*(coef(t_lower)[2:3] + coef(t_upper)[2:3])
  # SBC bias correction
  coef_sbc <- 2*coef(whole)[2:3] - bc_t
  return(coef_sbc)
}

y <- matrix(0, nrow = T2, ncol = N)
l1y <- matrix(0, nrow = T2, ncol = N)
coef_trend_b <- matrix(0, nrow = R, ncol = 2)
coef_trendcbc_b <- matrix(0, nrow = R, ncol = 2)
coef_trendsbc_b <- matrix(0, nrow = R, ncol = 2)
se_trend_sim <- matrix(0, nrow = R, ncol = 2)
# first non NA entry of l1lgdp for each country 
l1y[1, ] <- democracydata$l1lgdp[2 + c(0:(N - 1))*T]

# For linear time trend DGP
for (i in 1:R) {
  # Data simulation
  y[1, ] <- index_trend[1 + c(0:(N - 1))*T2] + coef_trend[2]*l1y[1, ] + rnorm(N, mean = 0, sd = sigma_trend)
  for (t in 2:T2) {
    # generate lagged variables
    l1y[t, ] <- y[t - 1, ]
    # update dependent variable at time t
    y[t, ] <- index_trend[t + c(0:(N - 1))*T2] + coef_trend[2]*l1y[t, ] + rnorm(N, mean = 0, sd = sigma_trend)
  }
  # Construct the dataset
  ys <- matrix(y, ncol = 1)
  l1ys <- matrix(l1y, ncol = 1)
  data_b <- data.frame(id = kronecker(c(1:N), rep(1, T2)), year = kronecker(rep(1, N), c(1:T2)), 
                       lgdp = ys, l1lgdp = l1ys, dem = fit_trend$x[, "dem"], 
                       timeindex = fit_trend$x[, "timeindex"])
  # Estimation using simulated data
  fit_trend_b <- lm(reg1, data_b)
  coef_trend_b[i, ] <- coef(fit_trend_b)[2:3]
  coef_trendcbc_b[i, ] <- cbc(data_b, formu = lgdp ~ dem + l1lgdp + factor(id) + timeindex)
  coef_trendsbc_b[i, ] <- sbc(data_b, formu = lgdp ~ dem + l1lgdp + factor(id) + timeindex)
  se_trend_sim[i, ] <- coef(summary(fit_trend_b))[2:3, 2]
}

# For quadratic time trend DGP
y_2 <- matrix(0, nrow = T2, ncol = N)
l1y_2 <- matrix(0, nrow = T2, ncol = N)
coef_trend2_b <- matrix(0, nrow = R, ncol = 2)
coef_trend2cbc_b <- matrix(0, nrow = R, ncol = 2)
coef_trend2sbc_b <- matrix(0, nrow = R, ncol = 2)
se_trend2_sim <- matrix(0, nrow = R, ncol = 2)
# first non NA entry of l1lgdp for each country 
l1y_2[1, ] <- democracydata$l1lgdp[2 + c(0:(N - 1))*T]

for (i in 1:R) {
  # Data simulation
  y_2[1, ] <- index_trend2[1 + c(0:(N - 1))*T2] + coef_trend2[2]*l1y_2[1, ] + rnorm(N, mean = 0, sd = sigma_trend2)
  for (t in 2:T2) {
    # generate lagged variables
    l1y_2[t, ] <- y_2[t - 1, ]
    # update dependent variable at time t
    y_2[t, ] <- index_trend2[t + c(0:(N - 1))*T2] + coef_trend2[2]*l1y_2[t, ] + rnorm(N, mean = 0, sd = sigma_trend2)
  }
  # Construct the dataset
  ys_2 <- matrix(y_2, ncol = 1)
  l1ys_2 <- matrix(l1y_2, ncol = 1)
  data2_b <- data.frame(id = kronecker(c(1:N), rep(1, T2)), year = kronecker(rep(1, N), c(1:T2)), 
                       lgdp = ys_2, l1lgdp = l1ys_2, dem = fit_trend2$x[, "dem"], 
                       timeindex = fit_trend$x[, "timeindex"])
  # Estimation using simulated data
  fit_trend2_b <- lm(reg1_quad, data2_b)
  coef_trend2_b[i, ] <- coef(fit_trend2_b)[2:3]
  coef_trend2cbc_b[i, ] <- cbc(data2_b, formu = lgdp ~ dem + l1lgdp + factor(id) + I(timeindex^2))
  coef_trend2sbc_b[i, ] <- sbc(data2_b, formu = lgdp ~ dem + l1lgdp + factor(id) + I(timeindex^2))
  se_trend2_sim[i, ] <- coef(summary(fit_trend2_b))[2:3, 2]
}


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

trend_nobc <- sumstatslinear(coef_trend_b, coef_trend, se_trend_sim)
trend_sbc <- sumstatslinear(coef_trendsbc_b, coef_trend, se_trend_sim)
trend_cbc <- sumstatslinear(coef_trendcbc_b, coef_trend, se_trend_sim)

trend2_nobc <- sumstatslinear(coef_trend2_b, coef_trend2, se_trend2_sim)
trend2_sbc <- sumstatslinear(coef_trend2sbc_b, coef_trend2, se_trend2_sim)
trend2_cbc <- sumstatslinear(coef_trend2cbc_b, coef_trend2, se_trend2_sim)



