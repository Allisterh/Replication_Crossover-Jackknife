# 7/12/2020 Shuowen Chen
# The script contains code that implement debiased fixed effects and panel bootstrap 
# standard errors for dataset in Chernozhukov, Kasahara and Schrimpf (2020)
# Serves as application for Chen, Chernozhukov, Fernandez-Val, Kasahara and Schripmf (2020)
# Acknowledgements: functions extend those written by Paul Schrimpf. 

library(lfe)
library(speedglm)
library(boot)
# load data (sdf)
load("~/Desktop/Research/Panel-Cross-over-jackknife/application/covidsince0307.RData")

# add weight to sdf
sdf$sweight <- 1
# create week factor to the data
sdf$week <- as.factor(strftime(sdf$date, format = "%V"))

#### Functions ####
# Create formula 
# The function extends the createfmla to accommodate FE specifications
# yvar: dependent variable for bpiy and piy
# xvars: regressors. 
#  (a) For bpiy, regressors are behavior, policy and information
#  (b) For pib, regressors are policy and information, dependent variable is 1 
#      of the 4 behavior variables
#  (c) For piy, regressors are policy and information variables
# iv: This is to accommodate the syntax of felm. "0" means no iv. 
# cluster: cluster at state level (treatment). 
# Note: interacted fixed effect hasn't been tested. 
createfmla_fe <- function(yvar, xvars, iv = "0", cluster = "state",
                          fe = c("state + month", "state * month", 
                                 "week + state", "date + state")) {
  fe <- match.arg(fe)
  rhs <- paste(xvars, collapse = " + ")
  if (fe == "state * month") {
    # we will use speedlm, which is faster
    rhs <- paste(rhs, "+", fe, collapse = " ")
  } else {
    rhs <- paste(rhs," | ", fe, " | ", iv, " | ", paste(cluster, collapse = " + ")) 
  }
  return(list(as.formula(paste(yvar, "~", rhs, sep = " ")), rhs)) 
}

# a function to conduct multiple split (along the cross section) bias correction
# to be called in the policyreg function
# s: number of splits
# n: number of cross section units
# dat: data in use
# uc: uncorrected estimators
# fm: regression formula
multiple_split <- function(s, n, unit, time, dat, uc, fm) {
  across <- 0 * uc
  for (k in 1:s) {
    sample1 <- sample(n, ceiling(n/2), replace = FALSE)
    subsample1 <- (unit %in% sample1 & time <= median(time)) | 
      (!(unit %in% sample1) & time >= median(time))
    subsample2 <- (unit %in% sample1 & time >= median(time)) | 
      (!(unit %in% sample1) & time <= median(time))
    cross1 <- felm(fm, data = dat[subsample1, ], 
                   weights = dat[subsample1, ]$sweight)
    cross2 <- felm(fm, data = dat[subsample2, ],
                   weights = dat[subsample2, ]$sweight)
    across <- across + ((coef(cross1) + coef(cross2))/2)/s
  }
  # average cross over corrected
  acbc <- 2 * uc - across
  return(acbc)
}

# Cross-over jackknife with overlapping in the tiem dimension
# q: fraction of time periods used in each subpabel,
#    default is 0.5. Larger than 0.5 features overlap in the 
#    time periods used: fraction is 2*q - 1
crossover <- function(unit, time, dat, uc, fm, q = 0.5) {
  t_cutoff <- c(floor(quantile(time, q)), 
                ceiling(quantile(time, 1 - q)))
  subsample1 <- (unit <= median(unit) & time <= t_cutoff[1]) | 
    (unit >= median(unit) & time >= t_cutoff[2])
  subsample2 <- (unit <= median(unit) & time >= t_cutoff[2]) | 
    (unit >= median(unit) & time <= t_cutoff[1])
  cross1 <- felm(fm, data = dat[subsample1, ], 
                 weights = dat[subsample1, ]$sweight)
  cross2 <- felm(fm, data = dat[subsample2, ], 
                 weights = dat[subsample2, ]$sweight)
  # debiased estimate
  bc <- uc/(1 - q) - 0.5*(coef(cross1) + coef(cross2))*q/(1 - q)
  return(bc)
}

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
    across <- across + ((coef(cross1) + coef(cross2))/2)/s
  }
  # average cross over corrected
  acbc <- uc/(1 - q) - across*q/(1 - q)
  return(acbc)
}

##### DFE estimation based on cross-over jackknife
# Note: The output of this function facilitates calculation of bootstrap se
reg_fe <- function(df, # data
                   yvar, # dependent variable (1 of 4 bevahiors if pib)
                   pols, # policies
                   bvars, # behavior (NULL if not behavior)
                   x, # controls 
                   iv, # iv for felm, if not needed just input "0" 
                   l = 14,
                   frac, # fraction of time periods used in subpanel for crossover
                   num_split, # number of sample split for multiple splitting bc
                   fixedeffect) {
  if (l == 0) {
    # This is for pib
    p <- pols
    b <- bvars
  } else {
    # This is for pbiy and piy
    p <- sprintf("lag(%s, %d)", pols, l)
    b <- sprintf("lag(%s, %d)", bvars, l)
  }
  # create regression formula
  xvars <- c(p, b, x)
  fmla <- createfmla_fe(yvar, xvars, fe = fixedeffect, iv = iv)
  
  # Uncorrected Estimator
  whole <- felm(fmla[[1]], data = df, weights = df$sweight)  
  coef_nobc <- coef(whole)
  
  # Various Bias Correction Estimators
  unit <- as.double(df$state)
  time <- as.double(df$date)
  # 1. Crossover Jackknife without overlap
  coef_cbc <- crossover(unit, time, df, coef_nobc, fmla[[1]], q = 0.5)
  # 2. Multiple split bias correction
  coef_acbc <- multiple_split(s = num_split, length(levels(df$state)), unit, 
                              time, df, coef_nobc, fmla[[1]])
  # 3. Cross-over with overlap in time dimension
  coef_overlap <- crossover(unit, time, df, coef_nobc, fmla[[1]], q = frac)
  
  # Sum of policy coefficients (uncorrected and corrected)
  peff <- c(sum(coef_nobc[p]), sum(coef_cbc[p]), 
            sum(coef_acbc[p]), sum(coef_overlap[p]))
  
  # Weighted sum of behavioral coefficients (uncorrected and corrected)
  if (!is.null(bvars)) {
    # weight is from April 1st to 10th
    w <- colMeans(subset(df, df$date >= as.Date("2020-04-01") &
                           df$date <= as.Date("2020-04-10"))[, bvars])
    beff <- c(sum(coef_nobc[b]*w), sum(coef_cbc[b]*w),
              sum(coef_acbc[b]*w), sum(coef_overlap[b]*w))
  } else {
    # for piy and pib no such metric
    beff <- NULL
  }
  # outputs: noncorrected, corrected, sum of policy and behavioral coefficients
  return(list(nobc = coef_nobc, cbc = coef_cbc, acbc = coef_acbc, 
              ovp = coef_overlap, sumpolicy = peff, sumbehavior = beff))
}


# The following function extends the function mainregressions
# NOTE: there is one change in the code
# the original code put xlist in the controls, but xlist is specified as an empty list
# in the program, so I'm not sure what it stands for. I think it's meant for creating a list
# so as to loop over choices of specifications. I deleted it in the argument. I checked the 
# results with and without xlist in the control for fe = "0" in the original code 
# and they are identical. 

mainregressions_fe <- function(df, # data
                               yvar, # dependent variable (Y)
                               pols, # policy variables
                               bvars, # behavior variables (B)
                               infovars, # information structures
                               tvars, # key confounder from SIR
                               ivlist = "0",
                               L = 14, 
                               frac = 0.6,
                               num_split = 5,
                               fixed) {
  # This is considering both pols and interacted of pmask with months
  plist <- list(pols, c("pmask.april","pmask.may", pols[-1]))
  # for pols only
  # plist <- list(pols)
  
  ijs <- expand.grid(1:length(plist), 1:length(infovars))
  
  # The function loops over plist, infovars to run different specifications
  pbiy <- apply(ijs, 1, function(ij) {
    reg_fe(df, yvar, plist[[ij[1]]], bvars, 
           c(sprintf("lag(%s, %d)", infovars[[ij[2]]], L), tvars), 
           ivlist, l = L, frac, num_split, fixedeffect = fixed)
  })

  piy <- apply(ijs, 1, function(ij) {
    reg_fe(df, yvar, plist[[ij[1]]], NULL, 
           c(sprintf("lag(%s, %d)", infovars[[ij[2]]], L), tvars), 
           ivlist, l = L, frac, num_split, fixedeffect = fixed)
  })
  
  # for pib, note the four behavioral variables each can be a 
  # dependent variable, so loop over bvars, plist and infovars
  ijs <- expand.grid(1:length(bvars), 1:length(plist))
  pib <- list()
  if (!is.null(infovars)) {
    infonum <- length(infovars) 
  } else {
    infonum <- 1
  }
  for (k in 1:infonum) {
    pib[[k]] <- apply(ijs, 1, function(ij) {
      reg_fe(df, bvars[ij[1]], plist[[ij[2]]], NULL,
             c(infovars[[k]]), ivlist, l = 0, frac, num_split, 
             fixedeffect = fixed)
    })
  }
  # return output as three lists
  return(list(pib = pib, pbiy = pbiy, piy = piy))
}

## Auxiliary function: produce data for nonparametric panel bootstrap
# Note: if data gets updated please change the date in the data_b$date command
data_rg <- function(data, mle) {
  n <- length(unique(data$state))
  t <- length(unique(data$date))
  # swap state index
  ids <- kronecker(sample.int(n, n, replace = TRUE), rep(1, t))
  data_b <- data[(ids - 1)*t + rep(c(1:t), n), ]
  data_b$state <- kronecker(c(1:n), rep(1, T))
  data_b$date <- rep(seq(as.Date("2020-03-07"), as.Date("2020-06-03"), by = "day"), n)
  data_b <- data.frame(data_b)
  #data_b <- pdata.frame(data_b, index = c("state", "date"))
  return(data_b)
}

# Function for weighted bootstrap
# More robust to close to singular designs because it resamples all the states
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

                                      #### Empirical Analysis ####
#### 1. Defining lhs and rhs variables
# dependent variable for pbiy and piy
yvar <- "dlogdc" # cases
yvar2 <- "dlogdd" # deaths
# policy variables
pols <- c("pmaskbus", "pk12", "pshelter", "pmovie", "prestaurant", "pnonessential") 
# behavioral variables
bvars <- c("workplaces", "retail", "grocery", "transit") 
# confounder from the SIR
tvars <- "dlogtests"
# information structure I and III 
# Note: to run only one structure, say I, just specify 
# infovars <- list(c("dlogdc", "logdc"))
infovars <- list(c("dlogdc", "logdc"),
                 c("dlogdc", "logdc", "dlogdc.national", "logdc.national")) # cases

infovarsd <- list(c("dlogdd", "logdd"),
                  c("dlogdd", "logdd", "dlogdd.national", "logdd.national")) # deaths

#### 2. Bootstrap Standard Errors 
# Compute estimates in each bootstrap
bootstat_case <- function(data) {
  # change the specification here
  tot <- mainregressions_fe(data, yvar, pols, bvars, infovars, tvars, fixed = "state + month")
  return(unlist(tot))
}

bootstat_deathcase <- function(data) {
  tot <- mainregressions_fe(data, yvar2, pols, bvars, infovars, NULL, L = 21, fixed = "state + month")
  return(unlist(tot))
}

bootstat_deathdeath <- function(data) {
  tot <- mainregressions_fe(data, yvar2, pols, bvars, infovarsd, NULL, L = 21, fixed = "state + month")
  return(unlist(tot))
}

# Call boot command to conduct bootstrap
set.seed(88) # seed for replication
num_boot <- 500 # number of bootstraps
ncores <- 1 # number of cpus (speed)

# nonpar will throw warning wrt multiple splitting bias correction
case_nonpar <- boot(data = sdf, statistic = bootstat_case, sim = "parametric", 
                    ran.gen = data_rg, mle = 0, parallel = "multicore", 
                    ncpus = ncores, R = num_boot)

case_wb <- boot(data = sdf, statistic = bootstat_case, sim = "parametric", 
                ran.gen = data_wb, mle = 0, parallel = "multicore", 
                ncpus = ncores, R = num_boot)

deathdeath_nonpar <- boot(data = sdf, statistic = bootstat_deathdeath, 
                          sim = "parametric", ran.gen = data_rg, mle = 0, 
                          parallel = "multicore", ncpus = ncores, R = num_boot)

deathdeath_wb <- boot(data = sdf, statistic = bootstat_deathdeath, 
                      sim = "parametric", ran.gen = data_wb, mle = 0, 
                      parallel = "multicore", ncpus = ncores, R = num_boot)

# a function that produces input for table output
bsemat <- function(bootlist) {
  # obtain the matrix of boostrap estimators
  result <- structure(vapply(bootlist$t, as.double, numeric(1)), dim = dim(bootlist$t))
  # compute the boostrap standard errors
  bse <- apply(result, 2, function(x) {
    # Normal scaled IQR
    return((quantile(x, 0.75, na.rm = TRUE) - 
              quantile(x, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25)))
  })
  # combine estimators and bse
  ref <- rbind(bootlist$t0, bse)
  # call starnum function to obtain level of siginficance
  starp <- c(0.1, 0.05, 0.01)
  num_star <- matrix(0, nrow = 1, ncol = length(bse))
  for (i in 1:length(bse)) {
    num_star[1, i] <- sum(pnorm(-abs(ref[1, i]), sd = ref[2, i]) < starp/2)
  }
  ref <- rbind(ref, num_star)
  return(ref)
}

