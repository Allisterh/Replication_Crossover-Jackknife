# 9/3/2020 Shuowen Chen
# This script implements calibrated exercise using democracy data
library(foreign)
library(plm)
library(sandwich)
library(boot)
library(speedglm)
###  read-in Democracy data
#filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
#setwd(filepath)

######### Bootstrap Bias Correction Functions
# 1. analytical bias correction
abc <- function(data, form, lags, N) {
  fit <- lm(form, data, x = TRUE, na.action = na.omit)
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

# 2. a function that randomly splits s times along the cross-section
# and pool to get mega half panels for crossover jackknife b
# Note: the dat has to be a panel (pdata.frame)
mega_data <- function(s, dat) {
  n <- length(levels(dat$id))
  unit <- as.double(dat$id)
  time <- as.double(dat$year)
  mega1 <- list()
  mega2 <- list()
  for (i in 1:s) {
    newid <- (i - 1)*n + unit
    index <- sample(n, ceiling(n/2), replace = FALSE)
    sub1 <- (unit %in% index & time <= median(time)) |
      (!(unit %in% index) & time >= median(time))
    sub2 <- (unit %in% index & time >= median(time)) | 
      (!(unit %in% index) & time <= median(time))
    mega1[[i]] <- dat[sub1, ]
    mega2[[i]] <- dat[sub2, ]
    mega1[[i]]$id <- newid[sub1]
    mega2[[i]]$id <- newid[sub2]
  }
  mega1 <- do.call(rbind, mega1)
  mega2 <- do.call(rbind, mega2)
  return(list(mega1, mega2))
}

# 3. Synthetic Panel bias correction
syntheticpanel <- function(dat, fm) {
  time <- as.double(dat$year)
  # relabel the ids
  newid <- paste0(dat[time >= median(time), ]$id, "20")
  # unfactor state so as to change 
  dat$id <- as.character(dat$id)
  dat[time >= median(time), ]$id <- newid
  # re-factor
  dat$id <- as.factor(dat$id)
  # Now estimate using the synthetic panel
  fit_syn <- lm(fm, data = dat)
  coef_syn <- coef(fit_syn)[2:6]
  return(coef_syn)
}

# 4. crossover jackknife panel allowing for overlap
crossover <- function(unit, time, dat, uc, uc_lr, 
                      fm, q = 0.5) {
  t_cutoff <- c(floor(quantile(time, q)), 
                ceiling(quantile(time, 1 - q)))
  subsample1 <- (unit <= median(unit) & time <= t_cutoff[1]) | 
    (unit >= median(unit) & time >= t_cutoff[2])
  subsample2 <- (unit <= median(unit) & time >= t_cutoff[2]) | 
    (unit >= median(unit) & time <= t_cutoff[1])
  cross1 <- lm(fm, data = dat[subsample1, ])
  cross2 <- lm(fm, data = dat[subsample2, ])
  # debiased estimate
  bc <- uc/(1 - q) - 0.5*(coef(cross1)[2:6] + coef(cross2)[2:6])*q/(1 - q)
  lr1 <- coef(cross1)[2]/(1 - sum(coef(cross1)[3:6]))
  lr2 <- coef(cross2)[2]/(1 - sum(coef(cross2)[3:6]))
  lr_bc <- uc_lr/(1 - q) - 0.5*(lr1 + lr2)*q/(1 - q) 
  return(list(bc, lr_bc))
}

# 5. crossover with multiple sample splitting
multiplecrossover <- function(s, n, unit, time, dat, 
                              uc, uc_lr, fm, q = 0.5) {
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
    cross1 <- lm(fm, data = dat[subsample1, ])
    cross2 <- lm(fm, data = dat[subsample2, ])
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
data.rg <- function(data, mle) {
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

# Estimators to be computed in each calibration and bootstrap estimation
estimator2 <- function(data, form) {
  # Fixed Effects
  fe_fit <- lm(form, data)
  coefs_fe <- coef(fe_fit)[2:6]
  lr_fe <- coefs_fe[1]/(1 - sum(coefs_fe[2:5]))

  # Various methods for bias correction
  year <- as.double(data$year)
  id <- as.double(data$id)
  # Split-sample bias correction
  fe_fit1 <- lm(form, data = data[(year <= median(year)), ])
  fe_fit2 <- lm(form, data = data[(year >= median(year)), ])
  coefs_sbc <- 19*coefs_fe/9 - 10*(coef(fe_fit1)[2:6] + coef(fe_fit2)[2:6])/18 
  lr_fe1 <- coef(fe_fit1)[2]/(1 - sum(coef(fe_fit1)[3:6]))
  lr_fe2 <- coef(fe_fit2)[2]/(1 - sum(coef(fe_fit2)[3:6]))
  lr_sbc <- 19*lr_fe/9 - 10*(lr_fe1 + lr_fe2)/18
  
  # Analytical bias correction
  bias_l4 <- abc(data, form, lags = 4, N = length(levels(data$id)))
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
  
  # Mega Split Crossover BC
  mega <- mega_data(s = 5, dat = data)
  fe_fit1 <- speedlm(form, mega[[1]])
  fe_fit2 <- speedlm(form, mega[[2]])
  coefs_mbc <- 19*coefs_fe/9 - 10*(coef(fe_fit1)[2:6] + coef(fe_fit2)[2:6])/18
  lr_fe1 <- coef(fe_fit1)[2]/(1 - sum(coef(fe_fit1)[3:6]))
  lr_fe2 <- coef(fe_fit2)[2]/(1 - sum(coef(fe_fit2)[3:6]))
  lr_mbc <- 19*lr_fe/9 - 10*(lr_fe1 + lr_fe2)/18
  
  # Synthetic Panel BC
  syn <- syntheticpanel(dat = data, fm = form)
  syn_bc <- 2*coefs_fe - syn
  lr_synbc <- 2*lr_fe - syn[1]/(1 - sum(syn[2:5]))
  
  # multiple sample splitting to crossover jackknife
  cbc5 <- multiplecrossover(s = 50, n = length(levels(data$id)), id, year,
                            data, coefs_fe, lr_fe, form, q = 10/19)
  coefs_cbc5 <- cbc5[[1]]
  lr_cbc5 <- cbc5[[2]]
  jbc5 <- multiplecrossover(s = 50, n = length(levels(data$id)), id, year,
                            data, coefs_fe, lr_fe, form, q = 12/19)
  coefs_jbc5 <- jbc5[[1]]
  lr_jbc5 <- jbc5[[2]]
  
  # compute analytical clustered standard errors
  fit <- plm(lgdp ~ dem + lag(lgdp, 1:4) - 1, data, model = "within", 
             effect = "twoways", index = c("id","year"))
  HCV_coefs <- vcovHC(fit, cluster = 'group')
  cse <- sqrt(diag(HCV_coefs)) 
  jac_lr <- c(1, rep(lr_fe, 4))/(1 - sum(coefs_fe[2:5]))
  cse_lr <- sqrt(t(jac_lr) %*% HCV_coefs %*% jac_lr)
  
  return(c(coefs_fe, coefs_sbc, coefs_abc, coefs_cbc, 
           coefs_jbc, coefs_mbc, syn_bc, coefs_cbc5,
           coefs_jbc5, lr_fe, lr_sbc, lr_abc, lr_cbc, 
           lr_jbc, lr_mbc, lr_synbc, lr_cbc5, lr_jbc5,
           cse, cse_lr))
}

# A function that calls computes bootstrap standard errors within in simulation
bse <- function(data, fm, ncores, btimes, bseed) {
  # store seeds set in the global environment
  old <- .Random.seed
  # upon exiting the function environment, 
  # restore to the original seed
  on.exit({.Random.seed <<- old})
  # within this function environment, set
  # the new seed
  set.seed(bseed, kind = "L'Ecuyer-CMRG")
  result <- boot(data = data, statistic = estimator2, sim = "parametric", 
                 ran.gen = data.rg, mle = 0, form = fm, 
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
  bestimate <- estimator2(data, formula)
  se <- bse(data, formula, ncores, btimes, bseed)
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
form <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + 
  factor(year) + factor(id)
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
data_calibration <- function(data, mle) {
  y  <- matrix(0, nrow = T, ncol = N)
  l1y  <- matrix(0, nrow = T, ncol = N)
  l2y  <- matrix(0, nrow = T, ncol = N)
  l3y  <- matrix(0, nrow = T, ncol = N)
  l4y  <- matrix(0, nrow = T, ncol = N)
  y[1, ] <- index[1 + c(0:(N - 1))*T] + coefs0["l1lgdp"]*l1y[1, ] + 
    coefs0["l2lgdp"]*l2y[1, ] + coefs0["l3lgdp"]*l3y[1, ] + 
    coefs0["l4lgdp"]*l4y[1, ] + rnorm(N, mean = 0, sd = sigma)
  for (t in 2:T) {
    l4y[t, ] <- l3y[t - 1, ]
    l3y[t, ] <- l2y[t - 1, ]
    l2y[t, ] <- l1y[t - 1, ]
    l1y[t, ] <- y[t - 1, ]
    y[t, ] <- index[t + c(0:(N - 1))*T] + coefs0["l1lgdp"]*l1y[t, ] + 
      coefs0["l2lgdp"]*l2y[t, ] + coefs0["l3lgdp"]*l3y[t, ] + 
      coefs0["l4lgdp"]*l4y[t, ] + rnorm(N, mean = 0, sd = sigma)
  }
  ys <- matrix(y, ncol = 1)
  l1ys <- matrix(l1y, ncol = 1)
  l2ys <- matrix(l2y, ncol = 1)
  l3ys <- matrix(l3y, ncol = 1)
  l4ys <- matrix(l4y, ncol = 1)
  data_b <- data.frame(id = kronecker(c(1:N), rep(1,T)), year = kronecker(rep(1, N), c(1:T)), 
                       lgdp = ys, l1lgdp = l1ys, l2lgdp = l2ys, l3lgdp = l3ys, l4lgdp = l4ys, 
                       dem = fe_fit$x[, "dem"])
  data_b <- pdata.frame(data_b, index = c("id", "year"))
  return(data_b)
}

data_calireshuffle <- function(data, mle) {
  y  <- matrix(0, nrow = T, ncol = N)
  l1y  <- matrix(0, nrow = T, ncol = N)
  l2y  <- matrix(0, nrow = T, ncol = N)
  l3y  <- matrix(0, nrow = T, ncol = N)
  l4y  <- matrix(0, nrow = T, ncol = N)
  y[1, ] <- index[1 + c(0:(N - 1))*T] + coefs0["l1lgdp"]*l1y[1, ] + 
    coefs0["l2lgdp"]*l2y[1, ] + coefs0["l3lgdp"]*l3y[1, ] + 
    coefs0["l4lgdp"]*l4y[1, ] + rnorm(N, mean = 0, sd = sigma)
  for (t in 2:T) {
    l4y[t, ] <- l3y[t - 1, ]
    l3y[t, ] <- l2y[t - 1, ]
    l2y[t, ] <- l1y[t - 1, ]
    l1y[t, ] <- y[t - 1, ]
    y[t, ] <- index[t + c(0:(N - 1))*T] + coefs0["l1lgdp"]*l1y[t, ] + 
      coefs0["l2lgdp"]*l2y[t, ] + coefs0["l3lgdp"]*l3y[t, ] + 
      coefs0["l4lgdp"]*l4y[t, ] + rnorm(N, mean = 0, sd = sigma)
  }
  ys <- matrix(y, ncol = 1)
  l1ys <- matrix(l1y, ncol = 1)
  l2ys <- matrix(l2y, ncol = 1)
  l3ys <- matrix(l3y, ncol = 1)
  l4ys <- matrix(l4y, ncol = 1)
  data_b <- data.frame(id = kronecker(sample(N), rep(1, T)), 
                       year = kronecker(rep(1, N), c(1:T)), lgdp = ys, 
                       l1lgdp = l1ys, l2lgdp = l2ys, l3lgdp = l3ys, 
                       l4lgdp = l4ys, dem = fe_fit$x[, "dem"])
  data_b <- pdata.frame(data_b, index = c("id", "year"))
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
            ran.gen = data_calireshuffle, R = time_cal, ncpus = core)

########## Printing results
est_dem <- bd$t[, c(1, 6, 11, 16, 21, 26, 31, 36, 41)]
est_l1 <- bd$t[, (1 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
est_l2 <- bd$t[, (2 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
est_l3 <- bd$t[, (3 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
est_l4 <- bd$t[, (4 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
est_lr <- bd$t[, c(46:54)]
est_cse_coef <- bd$t[, c(55:59)]
ese_cse_lr <- bd$t[, 60]
bse_dem <- bd$t[, (60 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
bse_l1 <- bd$t[, (60 + 1 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
bse_l2 <- bd$t[, (60 + 2 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
bse_l3 <- bd$t[, (60 + 3 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
bse_l4 <- bd$t[, (60 + 4 + c(1, 6, 11, 16, 21, 26, 31, 36, 41))]
bse_lr <- bd$t[, (60 + c(46:54))]
# We are interested in the properties of the following
# estimators: democracy, 4 lags of output, long run effect
table_simulation <- function(est, bse, ase, est0, 
                             app = c("dem", "lfp")) {
  app <- match.arg(app)
  tab <- matrix(0, nrow = 9, ncol = 9)
  colnames(tab) <- c('Bias', 'Std Dev', 'RMSE', 'BSE/SD', 
                     'ASE/SD', 'p.95 (BSE)', 'p.95 (ASE)', 
                     'Length (BSE)', 'Length (ASE)')
  if (app == "dem") {
    rownames(tab) <- c('FE', 'SBC', 'ABC', 'CBC', 'CBC-5/19', 
                       'MBC', 'Syn', 'CBC50', 'CBC50-5/19')
  } else {
    rownames(tab) <- c('FE', 'SBC', 'ABC', 'CBC', 'CBC-1/3', 
                       'MBC', 'Syn', 'CBC50', 'CBC50-1/3')
  }
  tab[, 1] <- 100*(apply(est, 2, mean)/est0 - 1)
  tab[, 2] <- 100*(apply(est/est0, 2, sd))
  tab[, 3] <- 100*sqrt((apply((est/est0 - 1)^2, 2, mean)))
  tab[, 4] <- apply((bse/apply(est, 2, sd)), 2, mean)
  tab[, 5] <- mean(ase)/apply(est, 2, sd)
  tab[, 6] <- apply((est + qnorm(.05/2)*bse <= est0) & 
                      (est + qnorm(1 - .05/2)*bse >= est0), 2, mean)
  tab[, 7] <- apply((est + qnorm(.05/2)*ase <= est0) & 
                      (est + qnorm(1 - .05/2)*ase >= est0), 2, mean)
  tab[, 8] <- 2*qnorm(1 - .05/2)*apply(bse, 2, mean)/abs(est0)
  tab[, 9] <- 2*qnorm(1 - .05/2)*mean(ase)/abs(est0)
  return(tab)
}

dem_co <- table_simulation(est_dem, bse_dem, est_cse_coef[, 1], coefs0[2])
l1lgdp_co <- table_simulation(est_l1, bse_l1, est_cse_coef[, 2], coefs0[3])
l2lgdp_co <- table_simulation(est_l2, bse_l2, est_cse_coef[, 3], coefs0[4])
l3lgdp_co <- table_simulation(est_l3, bse_l3, est_cse_coef[, 4], coefs0[5])
l4lgdp_co <- table_simulation(est_l4, bse_l4, est_cse_coef[, 5], coefs0[6])
lr_co <- table_simulation(est_lr, bse_lr, ese_cse_lr, lr0)

######## Save the workspace
#save.image(file = "democracymbc20.RData")

