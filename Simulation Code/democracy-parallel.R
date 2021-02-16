# 2/13/2021 Shuowen Chen
# This script implements calibrated exercise using democracy data
# with panel weighted and nonparametric bootstrap
library(foreign)
library(plm)
library(sandwich)
library(boot)
library(speedglm)
library(lfe)
###  read-in Democracy data
#filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
#setwd(filepath)

######### Bootstrap Bias Correction Functions
# 1. analytical bias correction
abc <- function(data, form, lags, N, w = NULL) {
  # Changing environment of the formula, otherwise lm cannot access the weights
  environment(form) <- environment()
  fit <- lm(form, data, x = TRUE, na.action = na.omit, weights = w)
  res <- fit$residuals
  jac <- solve(t(fit$x) %*% fit$x / length(res))[2:6,2:6]
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



# 2. Synthetic Panel bias correction
syntheticpanel <- function(dat, fm, w = NULL) {
  environment(fm) <- environment()
  time <- as.double(dat$year)
  # relabel the ids for the second half of the panel
  newid <- paste0(dat[time >= median(time), ]$id, "20")
  # unfactor state so as to change 
  dat$id <- as.character(dat$id)
  dat[time >= median(time), ]$id <- newid
  # re-factor
  dat$id <- as.factor(dat$id)
  # Now estimate using the synthetic panel
  fit_syn <- lm(fm, data = dat, weights = w)
  coef_syn <- coef(fit_syn)[2:6]
  return(coef_syn)
}

# 3. crossover jackknife panel allowing for overlap
crossover <- function(unit, time, dat, uc, uc_lr, 
                      fm, q = 0.5, w = TRUE) {
  environment(fm) <- environment()
  t_cutoff <- c(floor(quantile(time, q)), 
                ceiling(quantile(time, 1 - q)))
  subsample1 <- (unit <= median(unit) & time <= t_cutoff[1]) | 
    (unit >= median(unit) & time >= t_cutoff[2])
  subsample2 <- (unit <= median(unit) & time >= t_cutoff[2]) | 
    (unit >= median(unit) & time <= t_cutoff[1])
  if (isTRUE(w)) {
    cross1 <- lm(fm, data = dat[subsample1, ], weights = dat[subsample1, ]$weights)
    cross2 <- lm(fm, data = dat[subsample2, ], weights = dat[subsample2, ]$weights)
  } else {
    cross1 <- lm(fm, data = dat[subsample1, ])
    cross2 <- lm(fm, data = dat[subsample2, ])
  }
  # debiased estimate
  bc <- uc/(1 - q) - 0.5*(coef(cross1)[2:6] + coef(cross2)[2:6])*q/(1 - q)
  lr1 <- coef(cross1)[2]/(1 - sum(coef(cross1)[3:6]))
  lr2 <- coef(cross2)[2]/(1 - sum(coef(cross2)[3:6]))
  lr_bc <- uc_lr/(1 - q) - 0.5*(lr1 + lr2)*q/(1 - q) 
  return(list(bc, lr_bc))
}

# 4. crossover with multiple sample splitting
multiplecrossover <- function(s, n, unit, time, dat, 
                              uc, uc_lr, fm, q = 0.5, w = TRUE) {
  environment(fm) <- environment()
  # placeholders for bias correction terms
  across <- 0 * uc
  across_lr <- 0
  t_cutoff <- c(floor(quantile(time, q)), 
                ceiling(quantile(time, 1 - q)))
  for (k in 1:s) {
    sample1 <- sample(n, ceiling(n/2), replace = FALSE)
    subsample1 <- (unit %in% sample1 & time <= t_cutoff[1]) | 
      (!(unit %in% sample1) & time >= t_cutoff[2])
    subsample2 <- (unit %in% sample1 & time >= t_cutoff[2]) | 
      (!(unit %in% sample1) & time <= t_cutoff[1])
    if (isTRUE(w)) {
      cross1 <- lm(fm, data = dat[subsample1, ], weights = dat[subsample1, ]$weights)
      cross2 <- lm(fm, data = dat[subsample2, ], weights = dat[subsample2, ]$weights)
    } else {
      cross1 <- lm(fm, data = dat[subsample1, ])
      cross2 <- lm(fm, data = dat[subsample2, ])
    }
    across <- across + ((coef(cross1)[2:6] + coef(cross2)[2:6])/2)/s
    lr1 <- coef(cross1)[2]/(1 - sum(coef(cross1)[3:6]))
    lr2 <- coef(cross2)[2]/(1 - sum(coef(cross2)[3:6]))
    across_lr <- across_lr + ((lr1 + lr2)/2)/s
  }
  # average cross over corrected
  acbc <- uc/(1 - q) - across*q/(1 - q)
  acbc_lr <- uc_lr/(1 - q) - across_lr*q/(1 - q)
  return(list(acbc, acbc_lr))
}

# Generate bootstrap data for standard errors 
# Nonparametric panel bootstrap
data_rg <- function(data, mle) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  ids <- kronecker(sample.int(N, N, replace = TRUE), rep(1, T))
  data_b <- data[(ids - 1)*T + rep(c(1:T), N), ]
  data_b$id <- kronecker(c(1:N), rep(1, T))
  data_b$year <- rep(c(1:T), N)
  data_b <- data.frame(data_b)
  data_b <- pdata.frame(data_b, index = c("id", "year"))
  return(data_b)
}

# Panel weighted bootstrap
data_wb <- function(data, mle) {
  # number of states
  n <- length(unique(data$id))
  t <- length(unique(data$year))
  # Exponential weights
  multipliers <- rexp(n)
  # For each state, weight is the same for all dates
  weights <- rep(multipliers/sum(multipliers), each = t)
  # Add to the data and treat it as sampling weight
  data$weights <- weights
  return(data)
}


estimator <- function(data, form) {
  # set the environment of formula to pass the weights
  environment(form) <- environment()
  fit_b <- lm(form, data, weights = data$weights)
  coefs_fe <- coef(fit_b)[2:6]
  lr_fe <- coefs_fe[1]/(1 - sum(coefs_fe[2:5]))
  
  # Various methods for bias correction
  year <- as.double(data$year)
  id <- as.double(data$id)
  
  # Split-sample bias correction
  d1 <- data[(year <= median(year)), ]
  d2 <- data[(year >= median(year)), ]
  fe_fit1 <- lm(form, data = d1, weights = d1$weights)
  fe_fit2 <- lm(form, data = d2, weights = d2$weights)
  coefs_sbc <- 19*coefs_fe/9 - 10*(coef(fe_fit1)[2:6] + coef(fe_fit2)[2:6])/18 
  lr_fe1 <- coef(fe_fit1)[2]/(1 - sum(coef(fe_fit1)[3:6]))
  lr_fe2 <- coef(fe_fit2)[2]/(1 - sum(coef(fe_fit2)[3:6]))
  lr_sbc <- 19*lr_fe/9 - 10*(lr_fe1 + lr_fe2)/18
  
  # Analytical bias correction
  bias_l4 <- abc(data, form, lags = 4, N = length(levels(data$id)),
                 w = data$weights)
  coefs_abc <- coefs_fe - bias_l4
  jac_lr <- c(1, rep(lr_fe, 4))/(1 - sum(coefs_fe[2:5]))
  lr_abc <- lr_fe - crossprod(jac_lr, bias_l4)
  
  # Crossover without overlapping
  cbc <- crossover(id, year, data, coefs_fe, lr_fe, form, q = 10/19)
  coefs_cbc <- cbc[[1]]
  lr_cbc <- cbc[[2]]
  
  # CROSS-OVER JACKKNIFE with 5/19 overlapping
  jbc <- crossover(id, year, data, coefs_fe, lr_fe, form, q = 12/19)
  coefs_jbc <- jbc[[1]]
  lr_jbc <- jbc[[2]]
  
  # Synthetic Panel BC
  syn <- syntheticpanel(data, form, data$weights)
  syn_bc <- 2*coefs_fe - syn
  lr_synbc <- 2*lr_fe - syn[1]/(1 - sum(syn[2:5]))
  
  # multiple sample splitting to crossover jackknife
  # without overlapping
  cbc5 <- multiplecrossover(s = 2, n = length(levels(data$id)), id, year,
                            data, coefs_fe, lr_fe, form, q = 10/19)
  coefs_cbc5 <- cbc5[[1]]
  lr_cbc5 <- cbc5[[2]]
  # with overlapping
  jbc5 <- multiplecrossover(s = 2, n = length(levels(data$id)), id, year,
                            data, coefs_fe, lr_fe, form, q = 12/19)
  coefs_jbc5 <- jbc5[[1]]
  lr_jbc5 <- jbc5[[2]]
  
  # compute analytical clustered standard errors
  fit <- felm(lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp | id + year | 0 | id,
              data, weights = data$weights)
  cse <- coef(summary(fit))[, 2]
  # Note: the following code yields the same clustered standard errors
  # se <- sqrt(diag(vcov(fit)))
  # vcov(fit2) enables us to compute se for lr_fe using delta method
  jac_lr <- c(1, rep(lr_fe, 4))/(1 - sum(coefs_fe[2:5]))
  cse_lr <- sqrt(t(jac_lr) %*% vcov(fit) %*% jac_lr)
  
  # Note: output of plm with weights doesn't support vcovXX function, maybe a future
  # version of the package will fix it
  #fit <- plm(lgdp ~ dem + lag(lgdp, 1:4) - 1, data, model = "within", 
  #           effect = "twoways", index = c("id", "year"), weights = data$weights)
  #HCV_coefs <- vcovHC(fit, cluster = 'group')
  #cse <- sqrt(diag(HCV_coefs)) 
  #jac_lr <- c(1, rep(lr_fe, 4))/(1 - sum(coefs_fe[2:5]))
  #cse_lr <- sqrt(t(jac_lr) %*% HCV_coefs %*% jac_lr)
  
  return(c(coefs_fe, coefs_sbc, coefs_abc, coefs_cbc, 
           coefs_jbc, syn_bc, coefs_cbc5, coefs_jbc5, 
           lr_fe, lr_sbc, lr_abc, lr_cbc, lr_jbc, 
           lr_synbc, lr_cbc5, lr_jbc5, cse, cse_lr))
} 

# A function that calls computes weighted bootstrap standard errors within in simulation
bse_wb <- function(data, fm, ncores, btimes, bseed) {
  # store seeds set in the global environment
  old <- .Random.seed
  # upon exiting the function environment, 
  # restore to the original seed
  on.exit({.Random.seed <<- old})
  # within this function environment, set
  # the new seed
  set.seed(bseed, kind = "L'Ecuyer-CMRG")
  result <- boot(data = data, statistic = estimator, sim = "parametric", 
                 ran.gen = data_wb, mle = 0, form = fm, 
                 parallel = "multicore", ncpus = ncores, R = btimes)
  result <- structure(vapply(result$t, as.double, numeric(1)), 
                      dim = dim(result$t))
  se <- apply(result, 2, function(x) {
    return((quantile(x, .75, na.rm = TRUE) - 
              quantile(x, .25, na.rm = TRUE))/(qnorm(.75) - qnorm(.25)))
  })
  return(se)
}

# A function that calls computes nonpar bootstrap standard errors within in simulation
bse_nb <- function(data, fm, ncores, btimes, bseed) {
  # store seeds set in the global environment
  old <- .Random.seed
  # upon exiting the function environment, 
  # restore to the original seed
  on.exit({.Random.seed <<- old})
  # within this function environment, set
  # the new seed
  set.seed(bseed, kind = "L'Ecuyer-CMRG")
  result <- boot(data = data, statistic = estimator, sim = "parametric", 
                 ran.gen = data_rg, mle = 0, form = fm, 
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
  se <- bse_wb(data, formula, ncores, btimes, bseed)
  return(c(bestimate, se))
}

bootstat_ng <- function(data, formula, ncores, btimes, bseed){
  bestimate <- estimator(data, formula)
  se <- bse_nb(data, formula, ncores, btimes, bseed)
  return(c(bestimate, se))
}

########## Calibration
data <- read.dta("democracy-balanced-l4.dta")
# claim to be panel
data <- pdata.frame(data, index = c("id", "year"))
# add lagged dependent variables
data$l1lgdp <- lag(data$lgdp, 1)
data$l2lgdp <- lag(data$lgdp, 2) 
data$l3lgdp <- lag(data$lgdp, 3) 
data$l4lgdp <- lag(data$lgdp, 4) 
# add in weights
data$weights <- 1
form <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id)
fe_fit <- lm(form, data, x = TRUE)
index <- fitted.values(fe_fit) - coef(fe_fit)[3]*fe_fit$x[, 3] - 
  coef(fe_fit)[4]*fe_fit$x[, 4] - coef(fe_fit)[5]*fe_fit$x[, 5] - 
  coef(fe_fit)[6]*fe_fit$x[, 6]
sigma <- sqrt(sum(fe_fit$residuals^2)/fe_fit$df.residual)
coefs0 <- fe_fit$coefficients
lr0 <- coefs0[2]/(1 - sum(coefs0[3:6]))
T <- length(levels(data$year)) - 4 # takes out 4 lags of obs
N <- length(levels(data$id))

# Generate simulated data
dgp <- function(data, mle) {
  y <- matrix(0, nrow = T + 4, ncol = N)
  # set the initial 4 obs for each unit as the original data
  for (t in 1:4) y[t, ] <- data[(t + c(0:(N - 1))*(T + 4)), 6]
  index <- matrix(index, nrow = T, ncol = N)
  # use estimated coefficient to construct data
  for (t in 5:(T + 4)) {
    y[t, ] <- index[t - 4, ] + coefs0["l1lgdp"]*y[t - 1, ] + 
      coefs0["l2lgdp"]*y[t - 2, ] + coefs0["l3lgdp"]*y[t - 3, ] + 
      coefs0["l4lgdp"]*y[t - 4, ] + rnorm(N, mean = 0, sd = sigma)
  }
  y <- matrix(y, ncol = 1)
  data_b <- data.frame(id = kronecker(sample(N), rep(1, T + 4)), 
                       year = kronecker(rep(1, N), c(1:(T + 4))), lgdp = y, 
                       dem = data[, "dem"])
  data_b <- pdata.frame(data_b, index = c("id", "year"))
  data_b$l1lgdp <- lag(data_b$lgdp, 1)
  data_b$l2lgdp <- lag(data_b$lgdp, 2) 
  data_b$l3lgdp <- lag(data_b$lgdp, 3) 
  data_b$l4lgdp <- lag(data_b$lgdp, 4)
  return(data_b)
}


##########
time_cal <- 500 # number of simulations
core <- 28
seed_bse <- 13 # seed for bootstrap standard error estimation
set.seed(88, kind = "L'Ecuyer-CMRG") # seed for calibration
parallel::mc.reset.stream()
bd <- boot(data = data, statistic = bootstat, sim = "parametric", mle = 0, 
           formula = form, ncores = core, btimes = 200, bseed = seed_bse,
           ran.gen = dgp, R = time_cal, ncpus = core)

set.seed(88, kind = "L'Ecuyer-CMRG") # seed for calibration
parallel::mc.reset.stream()
bd2 <- boot(data = data, statistic = bootstat_ng, sim = "parametric", mle = 0, 
             formula = form, ncores = core, btimes = 200, bseed = seed_bse,
             ran.gen = dgp, R = time_cal, ncpus = core)

######### results tabulation
########## Printing results
est_dem <- bd$t[, c(1, 6, 11, 16, 21, 26, 31, 36)]
est_l1 <- bd$t[, (1 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_l2 <- bd$t[, (2 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_l3 <- bd$t[, (3 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_l4 <- bd$t[, (4 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_lr <- bd$t[, c(41:48)]
est_cse_coef <- bd$t[, c(49:53)]
ese_cse_lr <- bd$t[, 54]
bse_dem <- bd$t[, (54 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l1 <- bd$t[, (54 + 1 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l2 <- bd$t[, (54 + 2 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l3 <- bd$t[, (54 + 3 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l4 <- bd$t[, (54 + 4 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_lr <- bd$t[, (54 + c(41:48))]

est_dem2 <- bd2$t[, c(1, 6, 11, 16, 21, 26, 31, 36)]
est_l12 <- bd2$t[, (1 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_l22 <- bd2$t[, (2 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_l32 <- bd2$t[, (3 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_l42 <- bd2$t[, (4 + c(1, 6, 11, 16, 21, 26, 31, 36))]
est_lr2 <- bd2$t[, c(41:48)]
est_cse_coef2 <- bd2$t[, c(49:53)]
ese_cse_lr2 <- bd2$t[, 54]
bse_dem2 <- bd2$t[, (54 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l12 <- bd2$t[, (54 + 1 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l22 <- bd2$t[, (54 + 2 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l32 <- bd2$t[, (54 + 3 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_l42 <- bd2$t[, (54 + 4 + c(1, 6, 11, 16, 21, 26, 31, 36))]
bse_lr2 <- bd2$t[, (54 + c(41:48))]
# We are interested in the properties of the following
# estimators: democracy, 4 lags of output, long run effect
table_simulation <- function(est, bse, ase, est0,
                             app = c("dem", "lfp")) {
  app <- match.arg(app)
  tab <- matrix(0, nrow = 8, ncol = 9)
  colnames(tab) <- c('Bias', 'Std Dev', 'RMSE', 'BSE/SD',
                     'ASE/SD', 'p.95 (BSE)', 'p.95 (ASE)',
                     'Length (BSE)', 'Length (ASE)')
  if (app == "dem") {
    rownames(tab) <- c('FE', 'SBC', 'ABC', 'CBC', 'CBC-5/19',
                       'Syn', 'CBC50', 'CBC50-5/19')
  } else {
    rownames(tab) <- c('FE', 'SBC', 'ABC', 'CBC', 'CBC-1/3',
                       'Syn', 'CBC50', 'CBC50-1/3')
  }
  tab[, 1] <- 100*(apply(est, 2, mean)/est0 - 1)
  tab[, 2] <- 100*(apply(est/est0, 2, sd))
  tab[, 3] <- 100*sqrt((apply((est/est0 - 1)^2, 2, mean)))
  tab[, 4] <- apply(bse, 2, mean)/apply(est, 2, sd)
  tab[, 5] <- mean(ase)/apply(est, 2, sd)
  tab[, 6] <- apply((est + qnorm(.05/2)*bse <= est0) &
                      (est + qnorm(1 - .05/2)*bse >= est0), 2, mean)
  tab[, 7] <- apply((est + qnorm(.05/2)*ase <= est0) &
                      (est + qnorm(1 - .05/2)*ase >= est0), 2, mean)
  tab[, 8] <- 2*qnorm(1 - .05/2)*apply(bse, 2, mean)/abs(est0)
  tab[, 9] <- 2*qnorm(1 - .05/2)*mean(ase)/abs(est0)
  return(tab)
}

dem_w <- table_simulation(est_dem, bse_dem, est_cse_coef[, 1], coefs0[2])
l1lgdp_w <- table_simulation(est_l1, bse_l1, est_cse_coef[, 2], coefs0[3])
l2lgdp_w <- table_simulation(est_l2, bse_l2, est_cse_coef[, 3], coefs0[4])
l3lgdp_w <- table_simulation(est_l3, bse_l3, est_cse_coef[, 4], coefs0[5])
l4lgdp_w <- table_simulation(est_l4, bse_l4, est_cse_coef[, 5], coefs0[6])
lr_w <- table_simulation(est_lr, bse_lr, ese_cse_lr, lr0)

dem_n <- table_simulation(est_dem2, bse_dem2, est_cse_coef2[, 1], coefs0[2])
l1lgdp_n <- table_simulation(est_l12, bse_l12, est_cse_coef2[, 2], coefs0[3])
l2lgdp_n <- table_simulation(est_l22, bse_l22, est_cse_coef2[, 3], coefs0[4])
l3lgdp_n <- table_simulation(est_l32, bse_l32, est_cse_coef2[, 4], coefs0[5])
l4lgdp_n <- table_simulation(est_l42, bse_l42, est_cse_coef2[, 5], coefs0[6])
lr_n <- table_simulation(est_lr2, bse_lr2, ese_cse_lr2, lr0)


save.image(file = "democracy_weighted_nonpar.RData")
