# 2/11/2020 Shuowen Chen
# This script is a test file for the covid calibration

library(lfe)
library(speedglm)
library(boot)
library(plm)

load("C:/Users/Shuowen Chen/iCloudDrive/Desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code/covidsince0307.RData")

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
pols <- c("pmaskbus", "pk12", "pshelter", "pmovie", "prestaurant", "pnonessential") 
# deaths as information
infovarsd <- list(c("dlogdd", "logdd"),
                  c("dlogdd", "logdd", "dlogdd.national", "logdd.national")) 


###############
## regression specifications
vars <- c(pols, infovarsd[[1]])
regressors <- sprintf("%s%d", vars, 21)
rhs <- paste(regressors, collapse = " + ")
rhs_lm <- paste(rhs, "+", "state + week", collapse = " ")
# Regression formula
formula_piy <- as.formula(paste(yvar2, "~", rhs_lm, sep = " "))

# create lagged data
lagterm <- sapply(vars, function(x) out <- lag(sdf[, x], 21)) # note: lag is in the plm package
colnames(lagterm) <- regressors

# obtain data for regression
df <- cbind(sdf[, c("state", "week", yvar2, pols, "logdd", "sweight", "date")], lagterm)

# regression 
d_fit <- lm(formula_piy, data = df, x = TRUE, weights = df$sweight)
coef0 <- coef(d_fit)[regressors] # coefficients
# compute the part that is unchanged in the DGP construction
d_fitted_value <- d_fit$fitted.values
# Note: the timing of unchanged is from t = 1 till dates-21 for each i
unchanged <- d_fitted_value - d_fit$x[, c("dlogdd21", "logdd21")] %*% coef0[c("dlogdd21", "logdd21")]
# each column correpsonds to a state
unchanged <- matrix(unchanged, nrow = dates - 21, byrow = FALSE)
##############

# (2) Residuals 
res <- d_fit$residuals

# (a) ignore any structure 
sigma_ar <- sqrt(sum(d_fit$residuals^2)/d_fit$df.residual) 

# (b) Impose a panel AR structure
# Note: different from a pooled OLS in that we take panel structure into consideration
# when creating the lagged obs
# create lags
res_panel <- data.frame(id = kronecker(seq(N), rep(1, dates - 21)),
                        time = kronecker(rep(1, N), seq(dates - 21)), res = res)
res_panel <- pdata.frame(res_panel, index = c("id", "time"))
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
sigma_ar2 <- sqrt(sum(d_fit$residuals^2)/d_fit$df.residual) 

# (b) impose a panel MA structure: for each state, estimate MA coefficients
#res_mat <- matrix(res, nrow = dates - 21, byrow = FALSE)
#ar_coefs <- apply(res_mat, 2, function(x) coef(arma(x, order = c(21, 0), include.intercept = FALSE)))
#ma_coefs <- apply(res_mat, 2, function(x) coef(arma(x, order = c(0, 21), include.intercept = FALSE)))


hybrid_mbcovp <- function(s, n, unit, time, dat, uc, fm, q) {
  across <- 0 * uc
  t_cutoff <- c(floor(quantile(time, q)), 
                ceiling(quantile(time, 1 - q)))
  for (k in 1:s) {
    sample1 <- sample(n, ceiling(n/2), replace = FALSE)
    subsample1 <- (unit %in% sample1 & time <= t_cutoff[1]) | 
      (!(unit %in% sample1) & time >= t_cutoff[2])
    subsample2 <- (unit %in% sample1 & time >= t_cutoff[2]) | 
      (!(unit %in% sample1) & time <= t_cutoff[1])
    cross1 <- felm(fm, data = dat[subsample1, ], 
                   weights = dat[subsample1, ]$sweight)
    cross2 <- felm(fm, data = dat[subsample2, ],
                   weights = dat[subsample2, ]$sweight)
    across <- across + ((coef(cross1)[regressors] + coef(cross2)[regressors])/2)/s
  }
  # average cross over corrected
  acbc <- uc/(1 - q) - across*q/(1 - q)
  return(acbc)
}

#### 3. DGP Construction for each simulation
dgp_d <- function(data, mle) {
  # ar(21) structure
  res_sim <- matrix(0, nrow = dates, ncol = N) # placeholder
  u_it <- matrix(rnorm(N*dates, 0, sigma_ar), nrow = dates, ncol = N)
  res_sim[22:42, ] <- matrix(res, nrow = dates - 21, byrow = FALSE)[1:21, ] # initial obs
  for (t in 43:nrow(res_sim)) res_sim[t, ] <- arcoef %*% res_sim[(t - 1):(t - 21), ] + u_it[t, ]
  #for (i in 1:ncol(res_sim)) {
  #  for (t in 43:nrow(res_sim)) {
  #    res_sim[t, i] <- t(matrix(ar_coefs[, i])) %*% matrix(res_sim[(t - 1):(t - 21), i]) + u_it[t, i]
  #  }
  #}
  #for (t in 43:nrow(res_sim)) res_sim[t, ] <- arcoef %*% res_sim[(t - 1):(t - 21), ] + u_it[t, ]
  
  # ma(21) structure
  #res_sim2 <- matrix(0, nrow = dates, ncol = N) # placeholder
  #u_it <- matrix(rnorm(N*dates, 0, 1), nrow = dates, ncol = N)
  #res_sim2[22:42, ] <- matrix(res, nrow = dates - 21, byrow = FALSE)[1:21, ] # initial obs
  #for (t in 43:nrow(res_sim2)) res_sim2[t, ] <- arcoef %*% res_sim[(t - 1):(t - 21), ] + u_it[t, ]
  
  # synthetic logdd (an information variable)
  logdd_sim <- matrix(0, nrow = dates, ncol = N) # placeholder
  logdd_sim[1:21, ] <- matrix(data[, "logdd"], nrow = dates, byrow = FALSE)[1:21, ] # initial obs
  
  # synthetic dependent variable (dlogdd)
  y_sim <- matrix(0, nrow = dates, ncol = N) # placeholder
  y_sim[1:21, ] <- matrix(data[, "dlogdd"], nrow = dates, byrow = FALSE)[1:21, ] # initial obs
  u_iid <- matrix(rnorm(N*dates, 0, 1), nrow = dates, ncol = N)
  
  # Both logdd and dlogdd will be updated in the for loop
  for (t in 22:dates) {
  #  y_sim[t, ] <- unchanged[t - 21, ] + coef0["dlogdd21"]*y_sim[t - 21, ] + rnorm(N, mean =0, sd = sigma_ar2)
     y_sim[t, ] <- unchanged[t - 21, ] + coef0["dlogdd21"]*y_sim[t - 21, ] + coef0["logdd21"]*logdd_sim[t - 21, ] + res_sim[t, ]
    # update logdd_sim
     logdd_sim[t, ] <- y_sim[t, ] + logdd_sim[t - 7, ]
  }
  # construct the synthetic panel data
  y_sim <- matrix(y_sim, ncol = 1)
  logdd_sim <- matrix(logdd_sim, ncol = 1)
  #logdd_sim <- matrix(logdd_sim, ncol = 1)
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
  lagterms <- sapply(vars, function(x) out <- lag(data[, x], 21))
  colnames(lagterms) <- regressors
  dat_s <- cbind(data, lagterms)
  environment(fm) <- environment()
  d_fit_synthetic <- lm(fm, data = dat_s, weights = dat_s$sweight)
  
  # uncorrected
  coef_fe <- coef(d_fit_synthetic)[regressors]
  
  # crossover jackknife
  q <- 0.5
   unit <- as.double(dat_s$state)
   time <- as.double(dat_s$week)
    t_cutoff <- c(floor(quantile(time, q)), ceiling(quantile(time, 1 - q)))
    subsample1 <- (unit <= median(unit) & time <= t_cutoff[1]) | (unit >= median(unit) & time >= t_cutoff[2])
    subsample2 <- (unit <= median(unit) & time >= t_cutoff[2]) | (unit >= median(unit) & time <= t_cutoff[1])
    cross1 <- lm(fm, data = dat_s[subsample1, ], weights = dat_s[subsample1, ]$sweight)
    cross2 <- lm(fm, data = dat_s[subsample2, ], weights = dat_s[subsample2, ]$sweight)
    cbc <- coef_fe/(1 - q) - 0.5*(coef(cross1)[regressors] + coef(cross2)[regressors])*q/(1 - q)
  
  #cbc <- hybrid_mbcovp(10, 4539, unit, time, dat_s, coef_fe, fm, q)
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
                 ran.gen = data_wb, mle = 0, fm = fm, parallel = "multicore", 
                 ncpus = ncores, R = btimes)
  result <- structure(vapply(result$t, as.double, numeric(1)), dim = dim(result$t))
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


time_cal <- 50 # number of simulations
core <- 10
seed_bse <- 13 # seed for bootstrap standard error estimation
set.seed(88, kind = "L'Ecuyer-CMRG") # seed for calibration
#parallel::mc.reset.stream()
bd <- boot(data = df, statistic = bootstat, sim = "parametric", mle = 0, 
           formula = formula_piy, ncores = core, btimes = 10, bseed = seed_bse,
           ran.gen = dgp_d, R = time_cal, ncpus = core)

est_f2 <- bd$t[, 1:8]
est_c2 <- bd$t[, 9:16]
bse_f2 <- bd$t[, 17:24]
bse_c2 <- bd$t[, 25:32]

table_simulation <- function(est, bse, est0) {
  tab <- matrix(0, nrow = 5, ncol = length(est0))
  rownames(tab) <- c('Bias', 'Std Dev', 'RMSE', 'p.95 BSE Coverage', 'BSE/SD')
  colnames(tab) <- names(est0)
  tab[1, ] <- 100*(apply(est, 2, mean)/est0 - 1)
  tab[2, ] <- 100*(apply(t(est)/est0, 1, sd))
  tab[3, ] <- 100*sqrt((apply((t(est)/est0 - 1)^2, 1, mean)))
  tab[4, ] <- apply((est + qnorm(.05/2)*bse <= est0) &
                      (est + qnorm(1 - .05/2)*bse >= est0), 2, mean)
  tab[5, ] <- apply(bse, 2, mean)/apply(est, 2, sd)
  return(tab)
}

f2 <- table_simulation(est_f2, bse_f2, coef0)
cb2 <- table_simulation(est_c2, bse_c2, coef0)

