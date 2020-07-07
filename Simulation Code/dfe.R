# 7/3/2020 Shuowen Chen
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
#    time periods used. 
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
                               frac = 2/3,
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

# Regression
#case <- mainregressions_fe(sdf, yvar, pols, bvars, infovars, tvars, L = 14, fixed = "state + month")
#death_case <- mainregressions_fe(sdf, yvar2, pols, bvars, infovars, NULL, L = 21, fixed = "state + month")
#death_death <- mainregressions_fe(sdf, yvar2, pols, bvars, infovarsd, NULL, L = 21, fixed = "state + month")

# Week and state fixed effects
#case_week <- mainregressions_fe(sdf, yvar, pols, bvars, infovars, tvars, L = 14, fixed = "state + week")
#death_case_week <- mainregressions_fe(sdf, yvar2, pols, bvars, infovars, NULL, L = 21, fixed = "state + week")
#death_death_week <- mainregressions_fe(sdf, yvar2, pols, bvars, infovarsd, NULL, L = 21, fixed = "state + week")

# Consider date and state fixed effect. Drop information structure as there is no policy variation
#case_date <- mainregressions_fe(sdf, yvar, pols, bvars, NULL, tvars, L = 14, fixed = "state + date")
#death_date <- mainregressions_fe(sdf, yvar2, pols, bvars, NULL, NULL, L = 21, fixed = "state + date")


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
num_boot <- 5 # number of bootstraps
ncores <- 1 # number of cpus (speed)

# nonpar will throw warning wrt multiple splitting bias correction
case_nonpar <- boot(data = sdf, statistic = bootstat_case, sim = "parametric", 
                    ran.gen = data_rg, mle = 0, parallel = "multicore", 
                    ncpus = ncores, R = num_boot)

case_wb <- boot(data = sdf, statistic = bootstat_case, sim = "parametric", 
                ran.gen = data_wb, mle = 0, parallel = "multicore", 
                ncpus = ncores, R = num_boot)

deathcase_nonpar <- boot(data = sdf, statistic = bootstat_deathcase, 
                         sim = "parametric", ran.gen = data_rg, mle = 0, 
                         parallel = "multicore", ncpus = ncores, R = num_boot)

deathcase_wb <- boot(data = sdf, statistic = bootstat_deathcase, 
                     sim = "parametric", ran.gen = data_wb, mle = 0, 
                     parallel = "multicore", ncpus = ncores, R = num_boot)

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

#### Produce tables 

# auxiliary function for table productions
# 1. Filling tables
fill <- function(mat, covnames, type) {
  numcov <- length(covnames)
  if (type == "pib") {
    # First four columns: bvar1 
    tab <- matrix(NA, ncol = 16, nrow = 3*numcov)
    # fill in contents
    for (i in (1:16)) {
      tab[1:(3*numcov - 3), i] <- matrix(mat[, (1 + (i - 1)*(numcov - 1)):(i*(numcov - 1))], 
                                         ncol = 1)
    }
    for (i in (1:4)) {
      tab[(3*numcov - 2):(3*numcov), (4*i - 3):(4*i)] <- matrix(mat[, (4*i*numcov - 3):(4*i*numcov)],
                                                                nrow = 3, ncol = 4)
    }
  } else if (type == "pbiy") {
    tab <- matrix(NA, ncol = 4, nrow = 3*numcov)
    for (i in (1:4)) {
      tab[1:(3*numcov - 6), i] <- matrix(mat[, (1 + (i - 1)*(numcov - 2)):(i*(numcov - 2))],
                                         ncol = 1)
      # fill in sumpol and wsumb
      tab[(3*numcov - 5):(3*numcov - 3), i] <- matrix(mat[, 4*(numcov - 2) + i])
      tab[(3*numcov - 2):(3*numcov), i] <- matrix(mat[, 4*(numcov - 1) + i])
    }
  } else {
    tab <- matrix(NA, ncol = 4, nrow = 3*numcov)
    for (i in (1:4)) {
      tab[1:(3*numcov - 3), i] <- matrix(mat[, (1 + (i - 1)*(numcov - 1)):(i*(numcov - 1))],
                                         ncol = 1)
      # fill in sumpol
      tab[(3*numcov - 2):(3*numcov), i] <- matrix(mat[, 4*(numcov - 1) + i])
    }
  }
  # assign colnames and rownames
  colnames(tab) <- rep(c("nobc","bc","mbc","ovp"), ncol(tab)/4)
  rownames(tab) <- rep("bse", nrow(tab))
  for (i in (1:length(covnames))) {
    rownames(tab)[1 + 3*(i - 1)] <- covnames[i]
    rownames(tab)[3*i] <- "stars"
  }
  tab <- round(tab, digits = 3)
  # adding parentheses to bse
  for (i in 1:nrow(tab)) {
    if (rownames(tab)[i] == "bse") tab[i, ] <- paste0("(", format(unlist(tab[i, ])), ")")
  }
  return(tab)
}

pols_interact <- c("pmask.april","pmask.may", pols[-1])

# Depending on what regression type, produce a correspondent table
# lag: number of lags for covariates
tableproduction <- function(ref, type = c("pib", "pbiy", "piy"), pols, 
                            pols_interact, bvars, tvars, infovars, 
                            lag) {
  type <- match.arg(type)
  # rows for pib
  xpib1 <- c(pols, infovars[[1]], "sumpol")
  xpib3 <- c(pols, infovars[[2]], "sumpol")
  s1 <- length(xpib1)*16
  s2 <- length(xpib3)*16
  # Lagged covariates name for pbiy and piy
  lag_pol <- sprintf("lag(%s, %d)", pols, lag)
  lag_b <- sprintf("lag(%s, %d)", bvars, lag)
  lag_pol_int <- sprintf("lag(%s, %d)", pols_interact, lag)
  lag_info1 <- sprintf("lag(%s, %d)", infovars[[1]], lag)
  lag_info3 <- sprintf("lag(%s, %d)", infovars[[2]], lag)
  # rows for pbiy
  xpbiy1 <- c(lag_pol, lag_b, lag_info1, tvars, "sumpol", "wsumb")
  xpbiy3 <- c(lag_pol, lag_b, lag_info3, tvars, "sumpol", "wsumb")
  xpbiy1_int <- c(lag_pol_int, lag_b, lag_info1, tvars, "sumpol", "wsumb")
  xpbiy3_int <- c(lag_pol_int, lag_b, lag_info3, tvars, "sumpol", "wsumb")
  s3 <- length(xpbiy1)*4 # number for pbiy case info1
  s4 <- length(xpbiy1_int)*4
  s5 <- length(xpbiy3)*4
  s6 <- length(xpbiy3_int)*4
  # rows for piy
  xpiy1 <- c(lag_pol, lag_info1, tvars, "sumpol")
  xpiy3 <- c(lag_pol, lag_info3, tvars, "sumpol")
  xpiy1_int <- c(lag_pol_int, lag_info1, tvars, "sumpol")
  xpiy3_int <- c(lag_pol_int, lag_info3, tvars, "sumpol")
  piy_length <- c(length(xpiy1)*4, length(xpiy1_int)*4,
                  length(xpiy3)*4, length(xpiy3_int)*4)
  if (type == "pib") {
    # subset from ref
    pib_info1 <- ref[, 1:s1]
    pib_info3 <- ref[, (s1 + 161):(s1 + 160 + s2)]
    # return two tables: one for info1, another for info3
    # cols: four behaviors as dependent, each has nobc, bc, mbc and ovp
    # rows: 3 * num of pols and info (estimate, bse and stars)
    tab_info1 <- fill(pib_info1, xpib1, type = type)
    tab_info3 <- fill(pib_info3, xpib3, type = type)
    output <- list(pib_info1 = tab_info1, pib_info3 = tab_info3)
  } else if (type == "pbiy") {
    # subset from ref
    pbiy_info1 <- ref[, (s1 + s2 + 353):(s1 + s2 + 352 + s3)]
    pbiy_info1_int <- ref[, (s1 + s2 + 353 + s3):(s1 + s2 + 352 + s3 + s4)]
    pbiy_info3 <- ref[, (s1 + s2 + 353 + s3 + s4):(s1 + s2 + 352 + s3 + s4 + s5)]
    pbiy_info3_int <- ref[, (s1 + s2 + 353 + s3 + s4 + s5):(s1 + s2 + 352 + s3 + s4 + s5 + s6)]
    tab_info1 <- fill(pbiy_info1, xpbiy1, type = type)
    tab_info1_int <- fill(pbiy_info1_int, xpbiy1_int, type = type)
    tab_info3 <- fill(pbiy_info3, xpbiy3, type = type)
    tab_info3_int <- fill(pbiy_info3_int, xpbiy3_int, type = type)
    output <- list(pbiy_info1 = tab_info1, pbiy_info1_int = tab_info1_int,
                   pbiy_info3 = tab_info3, pbiy_info3_int = tab_info3_int)
  } else {
    # subset from ref
    initial <- s1 + s2 + 353 + s3 + s4 + s5 + s6
    piy_info1 <- ref[, initial:(initial + piy_length[1] - 1)]
    piy_info1_int <- ref[, (initial + piy_length[1]):(initial + sum(piy_length[1:2]) - 1)]
    piy_info3 <- ref[, (initial + sum(piy_length[1:2])):((initial + sum(piy_length[1:3]) - 1))]
    piy_info3_int <- ref[, (initial + sum(piy_length[1:3])):((initial + sum(piy_length) - 1))]
    tab_info1 <- fill(piy_info1, xpiy1, type = type)
    tab_info1_int <- fill(piy_info1_int, xpiy1_int, type = type)
    tab_info3 <- fill(piy_info3, xpiy3, type = type)
    tab_info3_int <- fill(piy_info3_int, xpiy3_int, type = type)
    output <- list(piy_info1 = tab_info1, piy_info1_int = tab_info1_int, 
                   piy_info3 = tab_info3, piy_info3_int = tab_info3_int)
  }
  return(output)
}

# Call bsemat to get inputs from the boot results
mat_case_nonpar <- bsemat(case_nonpar)
mat_case_wb <- bsemat(case_wb)
mat_deathcase_nonpar <- bsemat(deathcase_nonpar)
mat_deathcase_wb <- bsemat(deathcase_wb)
mat_deathdeath_nonpar <- bsemat(deathdeath_nonpar)
mat_deathdeath_wb <- bsemat(deathdeath_wb)

# Call producetable function to get table
case_nonpar_pib <- tableproduction(mat_case_nonpar, type = "pib", pols, pols_interact,
                                   bvars, tvars, infovars, 14)
death_nonpar_pib <- tableproduction(mat_deathcase_wb, type = "pib", pols, pols_interact, 
                               bvars, NULL, infovars, 21)

