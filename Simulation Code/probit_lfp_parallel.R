# 12/7/2020 Shuowen Chen
# This script implements calibrated exercise on LFP

### set working directory
filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
setwd(filepath)

library(foreign)
library(xtable)
library(alpaca)
library(boot)

############ Functions for Bias Corrections
# A function to implement split sample bias correction
sbc <- function(unit, time, dat, fm, uc) {
  # split along cross section
  deltahat11 <- feglm(fm, dat[(unit <= median(unit)), ], 
                      family = binomial(link = "probit"))
  betahat11 <- coef(deltahat11)
  deltahat12 <- feglm(fm, dat[(unit >= median(unit)), ], 
                      family = binomial(link = "probit"))
  betahat12 <- coef(deltahat12)
  # split along time series
  deltahat21 <- feglm(fm, dat[(time <= median(time)), ], 
                      family = binomial(link = "probit"))
  betahat21 <- coef(deltahat21)
  deltahat22 <- feglm(fm, dat[(time >= median(time)), ], 
                      family = binomial(link = "probit"))
  betahat22 <- coef(deltahat22)
  # bias correction based on Fernandez-Val and Weidner (2016)
  betahatsbc <- 3*uc - (betahat11 + betahat12 + betahat21 + betahat22)/2    
  return(betahatsbc)
}

# a function to implement multiple sample split sbc
# instead of partition along cross section by median, randomly split in half
sbc_multiple <- function(s, n, unit, time, fm, dat, uc) {
  # placeholder for bias correction term
  across <- 0*uc
  # split along the time series by median
  hat11 <- feglm(fm, dat[(time <= median(time)), ], 
                 family = binomial(link = "probit"))
  hat12 <- feglm(fm, dat[(time >= median(time)), ], 
                 family = binomial(link = "probit"))
  betahat11 <- coef(hat11)
  betahat12 <- coef(hat12)
  for (k in 1:s) {
    # randomly split the sample along the cross-section by half
    sample1 <- sample(n, ceiling(n/2), replace = FALSE)
    hat21 <- feglm(fm, dat[(unit %in% sample1), ], 
                   family = binomial(link = "probit"))
    hat22 <- feglm(fm, dat[!(unit %in% sample1), ], 
                   family = binomial(link = "probit"))
    betahat21 <- coef(hat21)
    betahat22 <- coef(hat22)
    across <- across + ((betahat11 + betahat12 + betahat21 + betahat22)/2)/s
  }
  # bias corrected estimate
  bc <- 3*uc - across
  return(bc)
}

# a function to create mega half panel
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

# A function to implement synthetic panel bias correction
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
  fit_syn <- feglm(fm, data = dat, 
                   family = binomial(link = "probit"))
  coef_syn <- coef(fit_syn)
  return(coef_syn)
}

# A function to implement crossover jackknife
# q denotes fraction of time used in subpanel
# q <- 0.5 means that there is no overlap, Larger than 0.5 
# features overlap in the time periods used: fraction is 2*q - 1
# For example, if fraction of overlap is 1/3, then q = 2/3.
crossover <- function(unit, time, dat, uc, fm, q = 0.5) {
  t_cutoff <- c(floor(quantile(time, q)), 
                ceiling(quantile(time, 1 - q)))
  subsample1 <- (unit <= median(unit) & time <= t_cutoff[1]) | 
    (unit >= median(unit) & time >= t_cutoff[2])
  subsample2 <- (unit <= median(unit) & time >= t_cutoff[2]) | 
    (unit >= median(unit) & time <= t_cutoff[1])
  cross1 <- feglm(fm, data = dat[subsample1, ], 
                  family = binomial(link = "probit"))
  cross2 <- feglm(fm, data = dat[subsample2, ], 
                  family = binomial(link = "probit"))
  # debiased estimate
  bc <- uc/(1 - q) - 0.5*(coef(cross1) + coef(cross2))*q/(1 - q)
  return(bc)
}

# Apply multiple sample splitting (along cross-section) to crossover
hybrid_mbcovp <- function(s, n, unit, time, dat, uc, fm, q = 0.5) {
  across <- 0 * uc
  t_cutoff <- c(floor(quantile(time, q)), 
                ceiling(quantile(time, 1 - q)))
  for (k in 1:s) {
    sample1 <- sample(n, ceiling(n/2), replace = FALSE)
    subsample1 <- (unit %in% sample1 & time <= t_cutoff[1]) | 
      (!(unit %in% sample1) & time >= t_cutoff[2])
    subsample2 <- (unit %in% sample1 & time >= t_cutoff[2]) | 
      (!(unit %in% sample1) & time <= t_cutoff[1])
    cross1 <- feglm(fm, data = dat[subsample1, ], 
                    family = binomial(link = "probit"))
    cross2 <- feglm(fm, data = dat[subsample2, ],
                    family = binomial(link = "probit"))
    across <- across + ((coef(cross1) + coef(cross2))/2)/s
  }
  # average cross over corrected
  acbc <- uc/(1 - q) - across*q/(1 - q)
  return(acbc)
}
############ Functions for Bootstrap
# statistics to be computed in each bootstrap draw 
estimators <- function(data, form) {
  # Fixed Effects
  deltahatfete <- feglm(form, data, 
                        family = binomial(link = "probit"))
  betahatfete <- coef(deltahatfete)
  # Various Bias Correction Methods
  time <- as.double(data$year)
  unit <- as.double(data$id)
  # 1. Split-sample bias correction
  betahatsbc <- sbc(unit, time, data, form, betahatfete)
  # 2. Split-sample with Multiple Splitting bias correction
  #betahatsbc10 <- sbc_multiple(s = 5, n = length(levels(data$id)), 
  #                             unit, time, form, data, betahatfete)
  # 3. Analytical bias correction
  deltahatbc <- biasCorr(deltahatfete, L = 2)
  betahatabc <- coef(deltahatbc)
  # 4. Crossover without overlapping;
  betahatcbc <- crossover(unit, time, data, betahatfete, form)
  # 5. Crossover with 1/3 overlapping;
  betahatjbc <- crossover(unit, time, data, betahatfete, form, q = 2/3)
  # 6. Crossover with multiple sample splitting
  betahatcbc5 <- hybrid_mbcovp(s = 50, n = length(levels(data$id)), unit,
                               time, data, betahatfete, form)
  betahatjbc5 <- hybrid_mbcovp(s = 50, n = length(levels(data$id)), unit,
                               time, data, betahatfete, form, q = 2/3)
  # 7. Mega sample split cross over
  mega <- mega_data(s = 2, dat = data)
  deltahat1 <- feglm(form, data = mega[[1]], 
                     family = binomial(link = "probit"))
  deltahat2 <- feglm(form, data = mega[[2]], 
                     family = binomial(link = "probit"))
  betahat1 <- coef(deltahat1)
  betahat2 <- coef(deltahat2)
  betahatmbc <- 2*betahatfete - (betahat1 + betahat2)/2
  # 8. Synthetic panel bias correction
  syn <- syntheticpanel(dat = data, fm = form)
  syn_bc <- 2*betahatfete - syn
  
  # Analytical standard errors for the uncorrected FE
  se_betahatfete <- summary(deltahatfete)[1]$cm[, 2]
  
  return(c(betahatfete, betahatsbc, betahatabc, 
           betahatcbc, betahatjbc, betahatcbc5, 
           betahatjbc5, betahatmbc, syn_bc, 
           se_betahatfete))
}

# Function to generate bootstrap data sets
# We use nonparametric panel bootstrap
data.rg <- function(data, mle) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  ids <- kronecker(sample.int(N, N, replace = TRUE), rep(1, T))
  data.b <- data[(ids - 1)*T + rep(c(1:T),N), ]
  data.b$id <- kronecker(c(1:N), rep(1, T))
  data.b$year <- rep(c(1:T), N)
  data.b <- data.frame(data.b)
  return(data.b)
}

# A function to compute bootstrap standard errors
bse <- function(data, form.fe, ncores, btimes, bseed) {
  # store seeds set in the global environment
  old <- .Random.seed
  # upon exiting the function environment, 
  # restore to the original seed
  on.exit({.Random.seed <<- old})
  # within this function environment, set
  # the new seed for bootstrap standard errors
  set.seed(bseed, kind = "L'Ecuyer-CMRG")
  result <- boot(data = data, statistic = estimators, sim = "parametric", 
                 ran.gen = data.rg, mle = 0, form = form.fe,
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
bootstat <- function(data, form.fe, ncores, btimes, bseed){
  bestimate <- estimators(data, form.fe)
  se <- bse(data, form.fe, ncores, btimes, bseed)
  return(c(bestimate, se))
}

################# Start of Calibration Exercise
lfpdata <- read.dta("LFPmovers.dta")
lfpdata <- pdata.frame(lfpdata, index = c("id", "year"))
# VALUES OF THE PARAMETERS;
N0 <- length(unique(lfpdata$id))
T0 <- length(unique(lfpdata$year))
P0 <- 7 # number of nonfe regressors in the specification
P <- 7
B <- 5 # number of bootstrap
N <- N0 # number of individuals
T <- T0 # number of time periods

# Regression specifications
spec_fe <- lfp ~ laglfp + kids0_2 + kids3_5 + kids6_17 + 
  loghusbandincome + age + age2 | id + year

# Check number of missing LFP values
nomissing <- !is.na(lfpdata$lfp)
pop <- length(lfpdata$lfp[nomissing])

# TRIMMING PARAMETERS FOR ANALYTICAL CORRECTIONS
L <- 2

# NUMBER OF DRAWS FOR RANDOM LEAVE-ONE-OUT JACKKNIFE;
draws <- 10

# CALIBRATION OF PARAMETERS
deltahatfete <- feglm(spec_fe, data = lfpdata, 
                      family = binomial(link = "probit"))
betahatfete <- deltahatfete$coefficients
# Compute theta*y_{i,t-1}+beta*X_{i,t}+\alpha_{i}+\gamma_{t}
index <- predict(deltahatfete)
index0 <- index - betahatfete[1] * lfpdata$laglfp
index1 <- index + betahatfete[1] * (1 - lfpdata$laglfp)

# Estimated coefficients from the data
beta0 <- betahatfete[c(1:P)]

# Interpretation: beta*X_{i,t}+\alpha_{i}+\gamma_{t}
true_index0 <- matrix(index0, nrow = T, byrow = FALSE)
# Interpretation: theta-theta*y_{i,t-1}+beta*X_{i,t}+\alpha_{i}+\gamma_{t}
true_index1 <- matrix(index1, nrow = T, byrow = FALSE)

rm(deltahatfete, index, index0, index1)

################# Functions to Simulate Data for Calibration
# Simulate data with original ID
data_calibration <- function(data, mle) {
  lfp_s <- matrix(0, nrow = T0, ncol = N0)
  lfp0 <- data$laglfp[data$year == 1] #initial value
  for (t in 1:T0) {
    # Interpretation: lfp0*betahatfete[1] + X'*betahatfete[2:end] > rnorm(N0),
    # Determines the lfp in the next period
    lfp_s[t, ] <- (lfp0*true_index1[t, ] + (1 - lfp0)*true_index0[t, ] > rnorm(N0))
    lfp0 <- lfp_s[t, ] # update so as to simulate the next period lfp
  }
  # simulate lagged dependent variable
  laglfp_s <- matrix(rbind(data$laglfp[data$year == 1], lfp_s[-T0, ]), 
                     ncol = 1, byrow = FALSE)
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp = lfp_s, laglfp = laglfp_s, kids0_2 = data$kids0_2, 
                       kids3_5 = data$kids3_5, kids6_17 = data$kids6_17, 
                       loghusbandincome = data$loghusbandincome, age = data$age, 
                       age2 = data$age2, id = data$id, year = data$year)
  # eliminate obs with all 0's or all 1's
  D <- kronecker(c(1:N0), matrix(1, T0, 1))
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T0, 1)) 
  insample <- ((td > 0) & (td < T0))
  data_s <- data_s[insample, ]
  return(data_s)
}

# Simulate data with IDs reshuffled
data_calireshuffle <- function(data, mle) {
  lfp_s <- matrix(0, nrow = T0, ncol = N0)
  lfp0 <- data$laglfp[data$year == 1] #initial value
  for (t in 1:T0) {
    # Interpretation: lfp0*betahatfete[1] + X'*betahatfete[2:end] > rnorm(N0),
    # Determines the lfp in the next period
    lfp_s[t, ] <- (lfp0*true_index1[t, ] + (1 - lfp0)*true_index0[t, ] > rnorm(N0))
    lfp0 <- lfp_s[t, ] # update so as to simulate the next period lfp
  }
  # simulate lagged dependent variable
  laglfp_s <- matrix(rbind(data$laglfp[data$year == 1], lfp_s[-T0, ]), 
                     ncol = 1, byrow = FALSE)
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp = lfp_s, laglfp = laglfp_s, kids0_2 = data$kids0_2, 
                       kids3_5 = data$kids3_5, kids6_17 = data$kids6_17, 
                       loghusbandincome = data$loghusbandincome, age = data$age, 
                       age2 = data$age2, id = kronecker(sample(N0), rep(1, T0)), 
                       year = data$year)
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T0, 1)) 
  insample <- ((td > 0) & (td < T0))
  # sort the data by id
  data_s <- data_s[order(data_s$id), ]
  data_s <- data_s[insample, ]
  return(data_s)
}

# For actual simulation set R and btimes to 500 and 200
core <- 28
seed_bse <- 13 # seed for bootstrap standard error estimation
set.seed(88, kind = "L'Ecuyer-CMRG") # seed for calibration
parallel::mc.reset.stream()
bda <- boot(data = lfpdata, statistic = bootstat, sim = "parametric", 
            mle = 0, form.fe = spec_fe, ncores = core, bseed = seed_bse,
            ran.gen = data_calibration, btimes = 200, R = 500, ncpus = core)

# a function that produces the diagnostic statistics
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

est_laglfp <- bda$t[, c(1, 8, 15, 22, 29, 36, 43, 50, 57)] 
est_kids02 <- bda$t[, 1 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
est_kids35 <- bda$t[, 2 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
est_kids617 <- bda$t[, 3 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
est_lhi <- bda$t[, 4 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
est_age <- bda$t[, 5 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
est_age2 <- bda$t[, 6 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
est_ase <- bda$t[, c(64:70)]

bse_laglfp <- bda$t[, 70 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
bse_kids02 <- bda$t[, 70 + 1 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
bse_kids35 <- bda$t[, 70 + 2 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
bse_kids617 <- bda$t[, 70 + 3 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
bse_lhi <- bda$t[, 70 + 4 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
bse_age <- bda$t[, 70 + 5 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]
bse_age2 <- bda$t[, 70 + 6 + c(1, 8, 15, 22, 29, 36, 43, 50, 57)]

laglfp_co <- table_simulation(est_laglfp, bse_laglfp, 
                              est_ase[, 1], betahatfete[1], app = "lfp")
kids02_co <- table_simulation(est_kids02, bse_kids02, 
                              est_ase[, 2], betahatfete[2], app = "lfp")
kids35_co <- table_simulation(est_kids35, bse_kids35, 
                              est_ase[, 3], betahatfete[3], app = "lfp")
kids617_co <- table_simulation(est_kids617, bse_kids617,
                               est_ase[, 4], betahatfete[4], app = "lfp")
lhi_co <- table_simulation(est_lhi, bse_lhi, est_ase[, 5], 
                           betahatfete[5], app = "lfp")
age_co <- table_simulation(est_age, bse_age, est_ase[, 6], 
                           betahatfete[6], app = "lfp")
age2_co <- table_simulation(est_age2, bse_age2, est_ase[, 7], 
                            betahatfete[7], app = "lfp")

