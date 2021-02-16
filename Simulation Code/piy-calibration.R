# 1/31/2021 Shuowen Chen
# This script conducts calibrated exercise using Covid data

library(lfe)
library(speedglm)
library(boot)
library(plm)

#filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
#setwd(filepath)
# load functions for regression, bias correction and bootstrap
#source('auxiliaryfunctions.R')
# load data (sdf)
load("covidsince0307.RData")

# add weight to sdf (for panel weighted bootstrap)
sdf$sweight <- 1
# create week factor to the data
sdf$week <- as.factor(strftime(sdf$date, format = "%V"))
N <- 51 # number of states
dates <- 89 # number of dates

#### 1. Define LHS and RHS variables
# dependent variable for pbiy and piy
yvar <- "dlogdc" # cases growth
yvar2 <- "dlogdd" # deaths growth
# policy variables
pols <- c("pmaskbus", "pk12", "pshelter", "pmovie", "prestaurant", "pnonessential") 
# behavioral variables
bvars <- c("workplaces", "retail", "grocery", "transit") 
# confounder from the SIR
tvars <- "dlogtests"
# information structure I and III 
# cases as information
infovars <- list(c("dlogdc", "logdc"),
                 c("dlogdc", "logdc", "dlogdc.national", "logdc.national")) 
# deaths as information
infovarsd <- list(c("dlogdd", "logdd"),
                  c("dlogdd", "logdd", "dlogdd.national", "logdd.national")) 

#### 2. Regression specification
# Notes: 
# (a) For calibration exercise, currently don't consider past cases/deaths at 
#     national level because I don't know how they were constructed
# (b) There are missing data in tvars. Total obs for estimation should be 3825,
#     but due to missing data there is only 3823. Therefore we first use deaths
#     as the dependent variable. 

# (1) Regression specification
#casereg <- mainregressions(sdf, yvar, pols, bvars, infovars[[1]], tvars, L = 14, 
#                           fixed = "state + week", spec = "piy", purpose = "cali")

#deathreg <- mainregressions(sdf, yvar2, pols, bvars, infovarsd[[1]], NULL, L = 21, 
#                            fixed = "state + week", spec = "piy", purpose = "cali")

# coefficients for calibration
#d_fit <- deathreg[[1]]$fit
#coef0 <- coef(d_fit) 

###############
## regression specifications
regressors <- sprintf("%s%d", c(pols, infovarsd[[1]]), 21)
rhs <- paste(regressors, collapse = " + ")
rhs_lm <- paste(rhs, "+", "state + week", collapse = " ")
# Regression formula
formula_piy <- as.formula(paste(yvar2, "~", rhs_lm, sep = " "))

# create lagged data
lagterm <- sapply(c(pols, infovarsd[[1]]), function(x) out <- lag(sdf[, x], 21))
colnames(lagterm) <- regressors

# obtain data for regression
df <- cbind(sdf[, c("state", "week", yvar2, pols, "logdd", "sweight", "date")], lagterm)

# regression 
d_fit <- lm(formula_piy, data = df, x = TRUE, weights = df$sweight)
coef0 <- coef(d_fit)[regressors] # coefficients
##############

# (2) Obtain the residuals and conduct an AR estimation
# Reason: allow for dependence in the errors given how the data are constructed
# Note: different from a pooled OLS in that we take panel structure into consideration
#       when creating the lagged obs
res <- d_fit$residuals
# construct a panel residual and take lags
res_panel <- data.frame(id = kronecker(seq(N), rep(1, dates - 21)),
                        time = kronecker(rep(1, N), seq(dates - 21)),
                        res = res)
res_panel <- pdata.frame(res_panel, index = c("id", "time"))
# create panel lags
placeholder <- matrix(0, nrow = length(res), ncol = 21)
colnames(placeholder) <- sapply(seq(21), function(x) sprintf("%s%d", "res", x))
for (l in 1:21) placeholder[, l] <- lag(res_panel$res, l)
# stack to a big panel with lagged obs
res_panel <- cbind(res_panel, placeholder)

# regression formula for AR residuals
rhs <- paste(colnames(placeholder), collapse = " + ")
res_fm <- as.formula(paste("res", "~", rhs, "-1", sep = " "))
fit_res <- lm(res_fm, res_panel, x = TRUE)
arcoef <- coef(fit_res) # ar coefficients
sigma_ar <- sqrt(sum(fit_res$residuals^2)/fit_res$df.residual) 

# compute the part that is unchanged in the DGP construction
d_fitted_value <- d_fit$fitted.values
# Note: the timing of unchanged is from t = 1 till dates-21 for each i
unchanged <- d_fitted_value - d_fit$x[, 8:9] %*% coef0[7:8]
unchanged <- matrix(unchanged, nrow = dates - 21, byrow = FALSE) # each column correpsonds to a state

#### 3. DGP Construction for each simulation
dgp_d <- function(data, mle) {
  # ar(21) structure
  res_sim <- matrix(0, nrow = dates, ncol = N) # placeholder
  u_it <- matrix(rnorm(N*dates, 0, sigma_ar), nrow = dates, ncol = N)
  res_sim[22:42, ] <- matrix(res, nrow = dates - 21, byrow = FALSE)[1:21, ] # initial obs
  for (t in 43:nrow(res_sim)) res_sim[t, ] <- arcoef %*% res_sim[(t - 1):(t - 21), ] + u_it[t, ]
  
  # synthetic logdd (an information variable)
  logdd_sim <- matrix(0, nrow = dates, ncol = N) # placeholder
  logdd_sim[1:21, ] <- matrix(data[, "logdd"], nrow = dates, byrow = FALSE)[1:21, ] # initial obs
  
  # synthetic dependent variable (dlogdd)
  y_sim <- matrix(0, nrow = dates, ncol = N) # placeholder
  y_sim[1:21, ] <- matrix(data[, "dlogdd"], nrow = dates, byrow = FALSE)[1:21, ] # initial obs
  u_iid <- matrix(rnorm(N*dates, 0, 1), nrow = dates, ncol = N)
  
  # Both logdd and dlogdd will be updated in the for loop
  for (t in 22:dates) {
    y_sim[t, ] <- unchanged[t - 21, ] + coef0[7]*y_sim[t - 21, ] + coef0[8]*logdd_sim[t - 21, ] + res_sim[t, ]
    # update logdd_sim
    logdd_sim[t, ] <- y_sim[t, ] + logdd_sim[t - 7, ]
  }
  # construct the synthetic panel data
  y_sim <- matrix(y_sim, ncol = 1)
  logdd_sim <- matrix(logdd_sim, ncol = 1)
  data_b <- data.frame(state = kronecker(seq(N), rep(1, dates)),
                       date = kronecker(rep(1, N), seq(dates)),
                       dlogdd = y_sim, logdd = logdd_sim)
  data_b <- pdata.frame(data_b, index = c("state", "date"))
  data_b <- cbind(data_b, data[, c(pols, "week", "sweight")])
  return(data_b)
}

# estimators to be computed in each simulation
estimator <- function(data, fm) {
  # create lagged data
  lagterms <- sapply(c(pols, infovarsd[[1]]), function(x) out <- lag(data[, x], 21))
  colnames(lagterms) <- regressors
  dat_s <- cbind(data, lagterms)
  environment(fm) <- environment() # so that weighted lm works
  d_fit_synthetic <- lm(fm, data = dat_s, weights = dat_s$sweight)
  
  # uncorrected
  coef_fe <- coef(d_fit_synthetic)[regressors]
  
  # crossover jackknife
  q <- 0.6
  unit <- as.double(dat_s$state)
  time <- as.double(dat_s$week)
  t_cutoff <- c(floor(quantile(time, q)), ceiling(quantile(time, 1 - q)))
  subsample1 <- (unit <= median(unit) & time <= t_cutoff[1]) | (unit >= median(unit) & time >= t_cutoff[2])
  subsample2 <- (unit <= median(unit) & time >= t_cutoff[2]) | (unit >= median(unit) & time <= t_cutoff[1])
  cross1 <- lm(fm, data = dat_s[subsample1, ], weights = dat_s[subsample1, ]$sweight)
  cross2 <- lm(fm, data = dat_s[subsample2, ], weights = dat_s[subsample2, ]$sweight)
  
  cbc <- coef_fe/(1 - q) - 0.5*(coef(cross1)[regressors] + coef(cross2)[regressors])*q/(1 - q)
  # return results
  return(c(coef_fe, cbc))
}

# function for bootstrap se
data_wb <- function(data, mle) {
  # number of states
  n <- length(unique(data$state))
  t <- length(unique(data$date))
  # Exponential weights
  multipliers <- rexp(n)
  # For each state, weight is the same for all dates
  weight <- rep(multipliers/sum(multipliers), each = t)
  # Add to the data and treat it as sampling weight
  data$sweight <- weight
  return(data)
}

bse <- function(data, fm, ncores, btimes, bseed) {
  # store seeds set in the global environment
  old <- .Random.seed
  # upon exiting the function environment, 
  # restore to the original seed
  on.exit({.Random.seed <<- old})
  # within this function environment, set
  # the new seed
  set.seed(bseed, kind = "L'Ecuyer-CMRG")
  result <- boot(data = data, statistic = estimator, sim = "parametric", 
                 ran.gen = data_wb, mle = 0, fm = fm, 
                 parallel = "multicore", ncpus = ncores, R = btimes)
  result <- structure(vapply(result$t, as.double, numeric(1)), 
                      dim = dim(result$t))
  se <- apply(result, 2, function(x) {
    return((quantile(x, .75, na.rm = TRUE) - 
              quantile(x, .25, na.rm = TRUE))/(qnorm(.75) - qnorm(.25)))
  })
  return(se)
}

# statistics to be computed in each bootstrap draw 
bootstat <- function(data, formula, ncores, btimes, bseed){
  bestimate <- estimator(data, formula)
  se <- bse(data, formula, ncores, btimes, bseed)
  return(c(bestimate, se))
}

time_cal <- 500 # number of simulations
core <- 28
seed_bse <- 13 # seed for bootstrap standard error estimation
set.seed(88, kind = "L'Ecuyer-CMRG") # seed for calibration
parallel::mc.reset.stream()
bd <- boot(data = df, statistic = bootstat, sim = "parametric", mle = 0, 
           formula = formula_piy, ncores = core, btimes = 200, bseed = seed_bse,
           ran.gen = dgp_d, R = time_cal, ncpus = core)


save.image(file = "piy.RData")

