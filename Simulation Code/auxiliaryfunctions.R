# 12/23/2020 Shuowen Chen
# This script stores some auxiliary functions for the crossover 
# jackknife bias correction project

###### Create regression formula 
# Inputs:
# yvar:    dependent variable for bpiy and piy
# xvars:   regressors. 
#          (a) For bpiy, regressors are behavior, policy and information
#          (b) For pib, regressors are policy and information, dependent variable
#              is 1 of the 4 behavior variables
#          (c) For piy, regressors are policy and information variables
# iv:      This is to accommodate the syntax of felm. "0" means no iv. 
# cluster: cluster at state level (treatment). 
# Outputs:
# fm1:     formula for felm 
# fm2:     formula for lm (for analytical bias correction) 
createfmla_fe <- function(yvar, xvars, iv = "0", cluster = "state",
                          fe = c("state + month", "state + week")) {
  fe <- match.arg(fe)
  rhs <- paste(xvars, collapse = " + ")
  rhs_felm <- paste(rhs," | ", fe, " | ", iv, " | ", paste(cluster, collapse = " + ")) 
  rhs_lm <- paste(rhs, "+", fe, collapse = " ")
  # formula for felm
  fm1 <- as.formula(paste(yvar, "~", rhs_felm, sep = " "))
  # formula for lm (for the use of analytical bias correction)
  fm2 <- as.formula(paste(yvar, "~", rhs_lm, sep = " "))
  return(list(fm1, fm2)) 
}

############# Different Bias Correction Functions
# Analytical bias correction for dynamic linear panel
# Inputs:
# fit:       output from the lm command
# trim:      trimming parameter for bias correction
# N:         number of cross section identity
# rhs:       regressors for bias correction
# unc:       uncorrecteed fixed effect estimators
#
# Output:
# corrected: analytically bias corrected estimators
abc <- function(data, fit, trim, n, rhs, unc) {
  res <- fit$residuals
  jac <- solve(t(fit$x) %*% fit$x / length(res))[rhs, rhs]
  indexes <- c(1:length(res))
  bscore <- rep(0, length(rhs))
  T <- length(res)/n
  for (i in 1:trim) {
    indexes <- indexes[-c(1 + c(0:(n - 1))*T)]
    lindexes <- indexes - i
    bscore <- bscore + t(fit$x[indexes, rhs]) %*% res[lindexes]/length(indexes)
  }
  bias <- -jac %*% bscore
  bias <- as.vector(bias/T)
  corrected <- unc - bias
  return(corrected)
}

# 2. Cross-over jackknife with option for overlapping in the time dimension
# Inputs:
# unit: number of cross section units
# time: number of time series units
# dat:  data
# uc:   uncorrected fixed effects estimators
# fm:   regression formula (felm)
# q:    fraction of time periods used in each subpabel. Default is 0.5, and there
#       is no overlap. Larger than 0.5 features overlap in the time periods used: 
#       fraction is 2*q - 1
# Outputs:
# bc:   bias corrected estimators
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

# 3. Multiple split (along the cross section) bias correction
# Inputs: 
# s: number of splits
# n: number of cross section units
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

# A variation that tends to stablize the multiple split
mega_split <- function(s, n, unit, time, dat, uc, fm) {
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
    # reassign the ids
    mega1[[i]]$state <- newid[sub1]
    mega2[[i]]$state <- newid[sub2]
  }
  mega1 <- do.call(rbind, mega1)
  mega2 <- do.call(rbind, mega2)
  cross1 <- felm(fm, data = mega1, weights = mega1$sweight)
  cross2 <- felm(fm, data = mega2, weights = mega2$sweight)
  cbc <- 2*uc - 0.5*(coef(cross1) + coef(cross2))
  return(cbc)
}


# a function that combines mbc and crossover with overlap 
# when q - 0.5, same as multiple_split
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

# The only difference from hybrid_mbcovp: instead of aggregating
# by average, take the median
hybrid_median <- function(s, n, unit, time, dat, uc, fm, q) {
  cross_mat <- matrix(0, nrow = length(uc), ncol = s)
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
    cross_mat[, k] <- (coef(cross1) + coef(cross2))/2
  }
  # median as correction
  mcross <- apply(cross_mat, 2, median)
  # average cross over corrected
  acbc <- uc/(1 - q) - mcross*q/(1 - q)
  return(acbc)
}

# syntheticpanel: relabel the id of states in the second 
# half of the time periods. This creates a synthetic panel
# with 2n states, T periods and T/2 observations per state
syntheticpanelbc <- function(dat, uc, fm) {
  time <- as.double(dat$date)
  # relabel the ids
  newid <- paste0(dat[time >= median(time), ]$state, "2")
  # unfactor state so as to change 
  dat$state <- as.character(dat$state)
  dat[time >= median(time), ]$state <- newid
  # re-factor
  dat$state <- as.factor(dat$state)
  # Now estimate using the synthetic panel
  fit_syn <- felm(fm, data = dat, weights = dat$sweight)
  syn_bc <- 2*uc - coef(fit_syn)
  return(syn_bc)
}

# a function that computes bootstrap bias correction
bbc <- function(uc, bdraws) {
  corrected <- 2*uc - apply(bdraws, 2, mean)
  bse <- apply(bdraws, 2, function(x) {
    # Normal scaled IQR
    return((quantile(x, 0.75, na.rm = TRUE) - 
              quantile(x, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25)))
  })
  ref <- rbind(corrected, bse)
  starp <- c(0.1, 0.05, 0.01)
  num_star <- matrix(0, nrow = 1, ncol = length(bse))
  for (i in 1:length(bse)) {
    num_star[1, i] <- sum(pnorm(-abs(ref[1, i]), sd = ref[2, i]) < starp/2)
  }
  return(rbind(ref, num_star))
}

##### DFE estimation based on cross-over jackknife
# Inputs:
# 1. df:        data
# 2. yvar:      dependent variable. Note: for pib regression it 
#               is 1 of the 4 behavior variables
# 3. plos:      policies
# 4. bvars:     behavior variables (NULL if runs piy regression)
# 5. x:         controls
# 6. l:         lags (14 for cases and 21 for deaths)
# 7. frac:      fraction of time periods used in subpanel for crossover
# 8. num_split: number of sample split for multiple splitting bias correction
# 9. purpose:   specify whether the regression is for calibration or estimation.
#               If "cali", then only report uncorrected estimator for DGP 
#               construction. If "est", will report other bias corrected estimators.
reg_fe <- function(df, yvar, pols, bvars, x, l = 14, frac, num_split, 
                   fixedeffect, purpose = c("cali", "est")) {
  purpose <- match.arg(purpose)
  if (l == 0) {
    # pib regression
    p <- pols
    b <- bvars
  } else {
    # pbiy or piy regression
    p <- sprintf("lag(%s, %d)", pols, l)
    b <- sprintf("lag(%s, %d)", bvars, l)
  }
  # create regression formula
  xvars <- c(p, b, x)
  fmla <- createfmla_fe(yvar, xvars, fe = fixedeffect, iv = "0")
  # Uncorrected Estimator
  whole <- felm(fmla[[1]], data = df, weights = df$sweight)  
  coef_nobc <- coef(whole)
  if (purpose == "cali") {
    peff <- sum(coef_nobc[p])
    # Weighted sum of behavioral coefficients
    if (!is.null(bvars)) {
      # weight is from April 1st to 10th
      w <- colMeans(subset(df, df$date >= as.Date("2020-04-01") &
                             df$date <= as.Date("2020-04-10"))[, bvars])
      beff <- sum(coef_nobc[b]*w)
    } else {
      # for piy and pib no such metric
      beff <- NULL
    }
    return(list(nobc = coef_nobc, sumpolicy = peff, sumbehavior = beff))
  } else {
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
}


# The main regression command
# Inputs:
#  1. df:        data
#  2. yvar:      dependent variable
#  3. pols:      policy variable
#  4. bvars:     behavior variables
#  5. infovars:  information variables
#  6. tvars:     key confounder from SIR (when yvar is case)
#  7. L:         lags
#  8. frac:      fraction of time periods used in subpanel for crossover
#  9. num_split: number of sample split for multiple splitting bias correction
#  10. spec:     specify which regression to run
#  11. purpose:  specify whether the regression is for calibration or estimation.
#                If "cali", then only report uncorrected estimator for DGP 
#                construction. If "est", will report other bias corrected estimators.
mainregressions <- function(df, yvar, pols, bvars, infovars, tvars, 
                            L = 14, frac = 0.6, num_split = 5, fixed,
                            spec = c("pbiy", "piy", "pib"), 
                            purpose = c("cali", "est")) {
  spec <- match.arg(spec)
  purpose <- match.arg(purpose)
  # This is considering both pols and interacted of pmask with months
  # plist <- list(pols, c("pmask.april","pmask.may", pols[-1]))
  # for pols only
  plist <- list(pols)
  
  # The function loops over plist, infovars to run different specifications
  ijs <- expand.grid(1:length(plist), 1:length(infovars))
  if (spec == "pbiy") {
    result <- apply(ijs, 1, function(ij) {
      reg_fe(df, yvar, plist[[ij[1]]], bvars, 
             c(sprintf("lag(%s, %d)", infovars[[ij[2]]], L), tvars), 
             l = L, frac, num_split, fixedeffect = fixed, purpose)
    })
  } else if (spec == "piy") {
    result <- apply(ijs, 1, function(ij) {
      reg_fe(df, yvar, plist[[ij[1]]], NULL, 
             c(sprintf("lag(%s, %d)", infovars[[ij[2]]], L), tvars), 
             l = L, frac, num_split, fixedeffect = fixed, purpose)
    })
  } else {
    # for pib, note the four behavioral variables each can be a 
    # dependent variable, so loop over bvars, plist and infovars
    ijs <- expand.grid(1:length(bvars), 1:length(plist))
    result <- list()
    if (!is.null(infovars)) {
      infonum <- length(infovars) 
    } else {
      infonum <- 1
    }
    for (k in 1:infonum) {
      result[[k]] <- apply(ijs, 1, function(ij) {
        reg_fe(df, bvars[ij[1]], plist[[ij[2]]], NULL, c(infovars[[k]]), 
               l = 0, frac, num_split, fixedeffect = fixed, purpose)
      })
    }
  }
  # return output 
  return(result)
}

############ Bootstrap Data Functions
# Produce data for nonparametric panel bootstrap
# Note: if data gets updated please change the date in the data_b$date command
data_rg <- function(data, mle) {
  n <- length(unique(data$state))
  t <- length(unique(data$date))
  # swap state index
  ids <- kronecker(sample.int(n, n, replace = TRUE), rep(1, t))
  data_b <- data[(ids - 1)*t + rep(c(1:t), n), ]
  data_b$state <- kronecker(c(1:n), rep(1, T))
  data_b$date <- rep(seq(as.Date("2020-03-07"), as.Date("2020-06-03"), 
                         by = "day"), n)
  data_b <- data.frame(data_b)
  return(data_b)
}

# Produce data for panel weighted bootstrap
# More robust to close to singular designs because it resamples 
# all the states
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

##### Table Production
# An intermediate output function
# Inputs:
# bootlist: output of boot command
# Outputs:
# A table containing:
# (1) Estimate using original data
# (2) Bootstrap Standard Errors
# (3) Significance levels (number of stars)
bsemat <- function(bootlist) {
  # obtain the matrix of boostrap estimators
  result <- structure(vapply(bootlist$t, as.double, numeric(1)), 
                      dim = dim(bootlist$t))
  # compute the boostrap standard errors
  bse <- apply(result, 2, function(x) {
    # Normal scaled IQR
    return((quantile(x, 0.75, na.rm = TRUE) - 
              quantile(x, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25)))
  })
  # combine estimators and bse
  ref <- rbind(bootlist$t0, bse)
  starp <- c(0.1, 0.05, 0.01)
  num_star <- matrix(0, nrow = 1, ncol = length(bse))
  for (i in 1:length(bse)) {
    num_star[1, i] <- sum(pnorm(-abs(ref[1, i]), sd = ref[2, i]) < starp/2)
  }
  ref <- rbind(ref, num_star)
  return(ref)
}

# A function to fill in entries for tables
fill2 <- function(mat, covnames, type, numest) {
  numcov <- length(covnames)
  if (type == "pib") {
    # First four columns: bvar1 
    tab <- matrix(NA, ncol = 4*numest, nrow = 3*numcov)
    # fill in contents
    for (i in 1:numest) {
      tab[1:(3*numcov - 3), i] <- matrix(mat[, (1 + (i - 1)*(numcov - 1)):(i*(numcov - 1))], 
                                         ncol = 1)
    }
    for (i in (numest + 1):(2*numest)) {
      tab[1:(3*numcov - 3), i] <- matrix(mat[, (1 + numest*numcov + (i - numest - 1)*(numcov - 1)):(numest*numcov + (i - numest)*(numcov - 1))], 
                                         ncol = 1)
    }
    for (i in (2*numest + 1):(3*numest)) {
      tab[1:(3*numcov - 3), i] <- matrix(mat[, (1 + 2*numest*numcov + (i - 2*numest - 1)*(numcov - 1)):(2*numest*numcov + (i - 2*numest)*(numcov - 1))],
                                         ncol = 1)
    }
    for (i in (3*numest + 1):(4*numest)) {
      tab[1:(3*numcov - 3), i] <- matrix(mat[, (1 + 3*numest*numcov + (i - 3*numest - 1)*(numcov - 1)):(3*numest*numcov + (i - 3*numest)*(numcov - 1))], 
                                         ncol = 1)
    }
    # fill in sum of policies
    for (i in 1:4) {
      tab[(3*numcov - 2):(3*numcov), 
          (numest*(i - 1) + 1):(numest*i)] <- matrix(mat[, (numest*(i*numcov - 1) + 1):(numest*i*numcov)],
                                                     nrow = 3, ncol = numest)
    }
  } else if (type == "pbiy") {
    tab <- matrix(NA, ncol = numest, nrow = 3*numcov)
    for (i in 1:ncol(tab)) {
      tab[1:(3*numcov - 6), i] <- matrix(mat[, (1 + (i - 1)*(numcov - 2)):(i*(numcov - 2))],
                                         ncol = 1)
      # fill in sumpol and wsumb
      tab[(3*numcov - 5):(3*numcov - 3), i] <- matrix(mat[, numest*(numcov - 2) + i])
      tab[(3*numcov - 2):(3*numcov), i] <- matrix(mat[, numest*(numcov - 1) + i])
    }
  } else {
    tab <- matrix(NA, ncol = numest, nrow = 3*numcov)
    for (i in 1:ncol(tab)) {
      tab[1:(3*numcov - 3), i] <- matrix(mat[, (1 + (i - 1)*(numcov - 1)):(i*(numcov - 1))],
                                         ncol = 1)
      # fill in sumpol
      tab[(3*numcov - 2):(3*numcov), i] <- matrix(mat[, numest*(numcov - 1) + i])
    }
  }
  # assign colnames and rownames
  #colnames(tab) <- rep(c("nobc","bc"), ncol(tab)/numest)
  rownames(tab) <- rep("bse", nrow(tab))
  for (i in (1:length(covnames))) {
    rownames(tab)[1 + 3*(i - 1)] <- covnames[i]
    rownames(tab)[3*i] <- "stars"
  }
  # round up to 3 decimal
  tab <- round(tab, digits = 3)
  for (i in 1:nrow(tab)) {
    # adding parentheses to bse
    if (rownames(tab)[i] == "bse") tab[i, ] <- paste0("(", format(unlist(tab[i, ])), ")")
    # change to stars
    if (rownames(tab)[i] == "stars") tab[i, ] <- sapply(1:length(tab[i, ]), function(x) {
      return(strrep("*", as.numeric(tab[i, ][x])))
    })
  }
  for (i in 1:nrow(tab)) {
    if (i %% 3 == 1) {
      tab[i, ] <- sapply(1:length(tab[i, ]), function(x) {
        if (tab[(i + 2), ][x] != "") {
          out <- paste0(tab[i, ][x], tab[(i + 2), ][x])
        } else {
          out <- tab[i, ][x]
        }
        return(out)
      })
    }
  }
  # finally remove the stars row
  # note: this uses Hmisc package
  tab <- tab[which(rownames(tab) %nin% "stars"), ]
  return(tab)
}

# Final function to produce tables
# Inputs:
# ref: output of the bsemat function
# type: whether to prduce results for pib, pbiy or piy regressions
# pols: policy variables
# bvars: behavioral variables
# tvars: SIR variable
# infovars: information variables (not a list, infovars[[1]])
# lag: lag of regressors (for printing names)
# num_est: number of estimators (uncorrected, corrected, etc.)
hybridtab <- function(ref, type = c("pib", "pbiy", "piy"), 
                      pols, bvars, tvars, infovars, lag,
                      num_est) {
  type <- match.arg(type)
  # rows for pbi
  xpib1 <- c(pols, infovars[[1]], "sumpol")
  xpib3 <- c(pols, infovars[[2]], "sumpol")
  # Here 4 denotes number of bvars
  pib_length <- c(length(xpib1)*4*num_est, length(xpib3)*4*num_est)
  # Lagged covariates name for pbiy and piy
  lag_pol <- sprintf("lag(%s, %d)", pols, lag)
  lag_b <- sprintf("lag(%s, %d)", bvars, lag)
  lag_info1 <- sprintf("lag(%s, %d)", infovars[[1]], lag)
  lag_info3 <- sprintf("lag(%s, %d)", infovars[[2]], lag)
  # rows for pbiy
  xpbiy1 <- c(lag_pol, lag_b, lag_info1, tvars, "sumpol", "wsumb")
  xpbiy3 <- c(lag_pol, lag_b, lag_info3, tvars, "sumpol", "wsumb")
  pbiy_length <- c(length(xpbiy1)*num_est, length(xpbiy3)*num_est)
  # rows for piy
  xpiy1 <- c(lag_pol, lag_info1, tvars, "sumpol")
  xpiy3 <- c(lag_pol, lag_info3, tvars, "sumpol")
  piy_length <- c(length(xpiy1)*num_est, length(xpiy3)*num_est)
  if (type == "pib") {
    # subset
    pib1 <- ref[, 1:pib_length[1]]
    pib3 <- ref[, (pib_length[1] + 1):sum(pib_length)]
    # call function fill2 to actually produce tables
    tab1 <- fill2(pib1, xpib1, type = type, num_est)
    tab3 <- fill2(pib3, xpib3, type = type, num_est)
    output <- list(info1 = tab1, info3 = tab3)
    #output <- list(info1 = pib1, info3 = pib3)
  } else if (type == "pbiy") {
    pbiy1 <- ref[,  (sum(pib_length) + 1):(sum(pib_length) + pbiy_length[1])]
    pbiy3 <- ref[, (sum(pib_length) + 
                      pbiy_length[1] + 1):(sum(pib_length) + sum(pbiy_length))] 
    tab1 <- fill2(pbiy1, xpbiy1, type = type, num_est)
    tab3 <- fill2(pbiy3, xpbiy3, type = type, num_est)
    output <- list(info1 = tab1, info3 = tab3)
  } else {
    piy1 <- ref[, (sum(pib_length) + sum(pbiy_length) + 
                     1):(sum(pib_length) + sum(pbiy_length) + piy_length[1])]
    piy3 <- ref[, (sum(pib_length) + sum(pbiy_length) + 
                     piy_length[1] + 1):(sum(pib_length) + sum(pbiy_length) + sum(piy_length))]
    tab1 <- fill2(piy1, xpiy1, type = type, num_est)
    tab3 <- fill2(piy3, xpiy3, type = type, num_est)
    output <- list(info1 = tab1, info3 = tab3)
  }
  return(output)
}

