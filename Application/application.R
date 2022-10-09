# 6/25/2022
# application for the crossover jackknife
# Shuowen Chen
#########################################################
#########################################################
library(plm)
library(readxl)
library(tidyverse)
library(data.table)
library(dplyr)
library(zoo)
library(boot)
library(lfe)
library(stargazer)
#########################################################
#########################################################
rootdir <- setwd("/Users/shuowenchen/downloads/covid-schools-master")
countyvars <- c("PopulationEstimate2018", "PopulationDensityperSqMile2010",
                "Smokers_Percentage", "Social.Association.Rate",
                "DiabetesPercentage", "MedianAge2010", "dem_to_rep_ratio",
                "HPSAShortage")

source(paste(rootdir,"rmd/generatetables.R",sep="/"))

# Loading stored data
load("/Users/shuowenchen/Downloads/covid-schools-master/data/data_preped.Rdata")

df <- as.data.frame(df)
df <- subset(df, df$date>=as.Date("2020-02-01"))
# define the subset data we use
startdate  <- "2020-04-01"
enddate  <- "2020-12-02"
sdf <- subset(df, df$date>=as.Date(startdate))
sdf <- subset(sdf, sdf$date<=as.Date(enddate))
sdf <- subset(sdf, sdf$fips!=0)
# define some new variables
sdf$pgather <- sdf$pgather50
sdf <- as.data.frame(sdf)
sdf$pschoolfull[sdf$portion.Unknown>0.5] <- NA
sdf$pschoolhybrid[sdf$portion.Unknown>0.5] <- NA
sdf$pschoolremote[sdf$portion.Unknown>0.5] <- NA
sdf$staffmask.No[sdf$staffmask.Unknown>0.5] <- NA
sdf$staff_dum <- sdf$staffmask.No
q <- 0
sdf$staff_dum[sdf$staffmask.No>q] <- 1
sdf$staff_dum[sdf$staffmask.No<=q] <- 0
sdf$staff_dum[is.na(sdf$staffmask.No)] <- NA
sdf$school_staff_dum <- sdf$staff_dum*sdf$school
sdf$full_staff_dum <- sdf$staff_dum*sdf$pschoolfull
sdf$hybrid_staff_dum <- sdf$staff_dum*sdf$pschoolhybrid

#########################################################
#########################################################
# summarize the number of counties in each state
fips_level <- as.character(unique(sdf$fips))
temp <- sapply(fips_level, function(x) head(strsplit(x, split="")[[1]], -3))
sum_cs <- as.data.frame(table(as.numeric(unlist(lapply(temp, function(x) paste0(x, collapse=""))))))

#########################################################
#########################################################
# Summary statistics
# This corresponds to the SI table in CKS
pols  <-  c("college", "school", "pschoolfull", "pschoolhybrid",
            "pschoolremote", "staff_dum", "pmask", "pgather50", "pshelter")
bvars <- c("fullwork", "partwork", "home", "bar", "restaurant", "gym", "church")
xvars <- c("dlogtests", "pop")
yvars <- c("dlogdc", "dlogdd", "logdc", "logdd", "d3logdd")
extravars <- c("school_staff_dum", "full_staff_dum", "hybrid_staff_dum")
sdf$pop <- sdf$PopulationEstimate2018/1000000
sdf$week <- week(sdf$date)
sdf2 <- subset(sdf, sdf$date>=as.Date("2020-04-15"))
sdf_temp <- sdf2[, c("date", "week", yvars, pols, bvars, xvars, extravars)]
sdf2 <- sdf2[,c(yvars,pols,bvars,xvars)]
tbl <- stargazer(sdf2)

#########################################################
#########################################################
# now convert daily data to weekly data by a 7-day moving average
# there are 3142 county, each will have 32 observations
# this is because we only include weeks that have all 7 days, thus ruling out week 16 and 49
# in our setting, starting date is 4/21 while ending date is 11/30
sdf_temp <- subset(sdf_temp, sdf_temp$week > 16)
sdf_temp <- subset(sdf_temp, sdf_temp$week < 49)
sdf_week <- matrix(0, nrow = 3142*32, ncol = 27)
for (i in 1:3142) {
  temp <- sdf_temp[(1+224*(i-1)):(224*i), ]
  sdf_week[(1+32*(i-1)):(32*i),] <- as.matrix(aggregate(temp[,-(1:2)], list(temp$week), FUN=mean))
}
sdf_week <- as.data.frame(sdf_week)
names(sdf_week) <- c("week", names(sdf_temp)[-(1:2)])
# add an id
sdf_week$fips <- rep(unique(sdf$fips),each = 32)
sdf_week <- pdata.frame(sdf_week, index = c("fips", "week"))
# report the summary statistics
tbl2 <- stargazer(sdf_week[,2:23])

# the original data have a lot of missing values, now do some preprocessing to
# obtain a balanced panel
## probably need to drop the counties with missing dependent variables
# missing values positions
detect_missing <- function(dat, x) {
  pos <- names(which(is.na(dat[, x])))
  final <- str_split(pos, "-")
  missing <- unique(unlist(lapply(final, `[[`, 1)))
  return(missing)
}

# for the case
tvars <- "dlogtests"
pols1  <-  c("pmask","pgather50","pshelter","college","school")
pols2  <-  c("pmask","pgather50","pshelter","college","school","school_staff_dum")
pols3  <-  c("pmask","pgather50","pshelter","college","pschoolfull","pschoolhybrid","pschoolremote")
pols4  <-  c("pmask","pgather50","pshelter","college","pschoolfull","pschoolhybrid","pschoolremote",
             "full_staff_dum","hybrid_staff_dum")
yvars <- "dlogdc"
vars_all <- as.matrix(c(tvars, pols4, yvars, "logdc", "school", "school_staff_dum"))
lala <- apply(vars_all, 1, detect_missing, dat = sdf_week)
county_missing <-unique(unlist(lala))
experiment <- base::setdiff(levels(sdf_week$fips), county_missing)
sdf_new <- sdf_week[sdf_week[, "fips"] %in% experiment, c(vars_all, "week")]
# this step is necessary as the droplevels does not work for this dataset
# need to redefine the levels so as to drop unused levels
sdf_new$fips <- as.factor(rep(experiment, each = 32))
# in total, 1901 counties

#######
# this part of the code explicitly adds the lags
sdf_new$l2logdc <- plm::lag(sdf_new$logdc, 2)
sdf_new$l3logdc <- plm::lag(sdf_new$logdc, 3)
sdf_new$l4logdc <- plm::lag(sdf_new$logdc, 4)
sdf_new$l2pmask <- plm::lag(sdf_new$pmask, 2)
sdf_new$l2pgather50 <- plm::lag(sdf_new$pgather50, 2)
sdf_new$l2pshelter <- plm::lag(sdf_new$pshelter, 2)
sdf_new$l2college <- plm::lag(sdf_new$college, 2)
sdf_new$l2school <- plm::lag(sdf_new$school, 2)
sdf_new$l2school_staff_dum <- plm::lag(sdf_new$school_staff_dum, 2)

#########################################################
#########################################################
# main regression function

### Using speedlm as a first check, replaced by felm as the latter obtains the same
# result but is faster
####
form <- dlogdc ~ l2pmask + l2pgather50 + l2pshelter + l2college + l2school + dlogtests + l2logdc + l3logdc + l4logdc + factor(fips) + factor(week)
robust <- speedlm(form, sdf_new)

f2 <- dlogdc ~ l2pmask + l2pgather50 + l2pshelter + l2college + l2school +
  dlogtests + l2logdc + l3logdc + l4logdc | fips + week | 0 | fips

f_pos2 <- dlogdc ~ l2pmask + l2pgather50 + l2pshelter + l2college + l2school + l2school_staff_dum +
  dlogtests + l2logdc + l3logdc + l4logdc | fips + week | 0 | fips

fe_uc <- felm(f_pos2, sdf_new, keepX = TRUE)

# now do crossover jackknife
multiple_split <- function(s, dat, uc, fm) {
  n <-length(unique(dat$fips))
  unit <- as.double(dat$fips)
  time <- as.double(dat$week)
  across <- 0 * uc
  # store number of splitting
  k <- 0
  while (k < s) {
    sample1 <- sample(n, floor(n/2), replace = FALSE)
    subsample1 <- (unit %in% sample1 & time < median(time)) |
      (!(unit %in% sample1) & time >= median(time))
    subsample2 <- (unit %in% sample1 & time > median(time)) |
      (!(unit %in% sample1) & time <= median(time))
    cross1 <- felm(fm, data = dat[subsample1, ])
    cross2 <- felm(fm, data = dat[subsample2, ])
    across <- across + ((coef(cross1) + coef(cross2))/2)/s
    k <- k + 1
  }
  # average cross over corrected
  acbc <- 2 * uc - across
  return(list(acbc = acbc, uc = uc, cross1 = dat[subsample1, ], cross2 = dat[subsample2, ]))
}

testing <- multiple_split(1, sdf_new, coef(fe_uc), f_pos2)
testing

# inputs:
# fm: regression formula
# dat: full panel data
# weights: regression weights
# s: number of splits
piy_2 <- function(fm, dat, weights = NULL, s) {
  piy <- felm(fm, dat, weights = weights, keepX = TRUE)
  est_fe <- coef(piy) # fixed effect estimators

  id <- as.double(dat$fips)
  time <- as.double(dat$week)
  n <- length(unique(sdf_new$fips)) # number of cross section units
  cbc_multiple <- multiple_split(s, n, id, time, dat, est_fe, fm)
  out <- list(fe = est_fe, bc = cbc_multiple, X = piy$X)
  return(out)
}

#piy_analysis <- function(yvar, dat, infovar, policyvar, exo_test,
#                         cluster = "fips", fe = "fips + week",
#                         lags = 2, debias = FALSE, weights = NULL) {
  # first, create formula
#  iv <- "0" # this is for felm formula
  # lags
#  if (lags != 0) policyvar <- sprintf("lag(%s, %d)", policyvar, lags)
#  xvars <- c(policyvar, exo_test, infovar)
#  rhs <- paste(xvars, collapse = " + ")
#  rhs_felm <- paste(rhs, " | ", fe, " | ", iv, " | ", paste(cluster, collapse = " + "))
#  fm <- as.formula(paste(yvar, "~", rhs_felm, sep = " "))
  # second, invoke felm to obtain the results
#  piy <- felm(fm, dat, weights = weights, keepX = TRUE)
#  if (debias == FALSE) {
#    return(piy)
#  } else {
#    id <- as.double(dat$fips)
#    time <- as.double(dat$week)
#    sub1 <- (id <= median(id) & time <= median(time)) |
#      (id >= median(id) & time >= median(time))
#    sub2 <- (id <= median(id) & time >= median(time)) |
#      (id >= median(id) & time <= median(time))
#    cross1 <- felm(fm, dat[sub1, ], weights = dat[sub1, ]$weights)
#    cross2 <- felm(fm, dat[sub2, ], weights = dat[sub2, ]$weights)
#    bc <- 2*piy$coefficients - 0.5*(cross1$coefficients + cross2$coefficients)
#    out <- list(fe = piy$coefficients, bc = bc, form = fm, X = piy$X)
#    return(out)
#  }
#}

#########################################################
#########################################################
# now run the regression using weekly data

### Case growth

tvars <- "dlogtests"
pols1  <-  c("pmask","pgather50","pshelter","college","school")
pols2  <-  c("pmask","pgather50","pshelter","college","school","school_staff_dum")
pols3  <-  c("pmask","pgather50","pshelter","college","pschoolfull","pschoolhybrid","pschoolremote")
pols4  <-  c("pmask","pgather50","pshelter","college","pschoolfull","pschoolhybrid","pschoolremote",
             "full_staff_dum","hybrid_staff_dum")


yvars <- "dlogdc"


info <- c("lag(logdc, 2)", "lag(logdc, 3)", "lag(logdc, 4)")


results_pols1 <- piy_analysis(yvars, info, pols1, tvars,
                              dat = sdf_new, debias = TRUE)



results_pols2 <- piy_analysis(yvars, info, pols2, tvars,
                              dat = sdf_week, debias = TRUE)
results_pols3 <- piy_analysis(yvars, info, pols3, tvars,
                              dat = sdf_week, debias = TRUE)
results_pols4 <- piy_analysis(yvars, info, pols4, tvars,
                              dat = sdf_week, debias = TRUE)

### Death growth
info_d <- c("lag(logdd, 5)", "lag(logdd, 6)", "lag(logdd, 7)")
yvars_d <- "d3logdd"

results_d_pols1 <- piy_analysis(yvars_d, info_d, pols1, exo_test = NULL,
                                dat = sdf_week, lags = 5, debias = TRUE)
results_d_pols2 <- piy_analysis(yvars_d, info_d, pols2, exo_test = NULL,
                                dat = sdf_week, lags = 5, debias = TRUE)
results_d_pols3 <- piy_analysis(yvars_d, info_d, pols3, exo_test = NULL,
                                dat = sdf_week, lags = 5, debias = TRUE)
results_d_pols4 <- piy_analysis(yvars_d, info_d, pols4, exo_test = NULL,
                                dat = sdf_week, lags = 5, debias = TRUE)

#########################################################
# Generate bootstrap data for standard errors

# Nonparametric panel bootstrap
data_rg <- function(dat, mle) {
  n <- length(unique(dat$fips))
  T <- length(unique(dat$week))
  ids <- kronecker(sample.int(n, n, replace = TRUE), rep(1, T))
  dat_b <- dat[(ids - 1)*T + rep(c(1:T), n), ]
  dat_b$fips <- kronecker(c(1:n), rep(1, T))
  dat_b$week <- rep(c(1:T), n)
  dat_b <- data.frame(dat_b)
  dat_b <- pdata.frame(dat_b, index = c("fips", "week"))
  return(dat_b)
}

# Panel weighted bootstrap
data_wb <- function(dat, mle) {
  # number of states
  n <- length(unique(dat$fips))
  t <- length(unique(dat$week))
  # Exponential weights
  multipliers <- rexp(n)
  # For each state, weight is the same for all dates
  weights <- rep(multipliers/sum(multipliers), each = t)
  # Add to the data and treat it as sampling weight
  dat$weights <- weights
  return(dat)
}

# estimator in each bootstrap iteration
estimator <- function(dat_sim) {
  # essentially, use the bootstrap data and piy_analysis function
  case_pols1 <- piy_analysis(yvars, info, pols1, tvars, dat = dat_sim,
                             debias = TRUE, weights = dat_sim$weights)
  case_pols2 <- piy_analysis(yvars, info, pols2, tvars, dat = dat_sim,
                             debias = TRUE, weights = dat_sim$weights)
  case_pols3 <- piy_analysis(yvars, info, pols3, tvars, dat = dat_sim,
                             debias = TRUE, weights = dat_sim$weights)
  case_pols4 <- piy_analysis(yvars, info, pols4, tvars, dat = dat_sim,
                             debias = TRUE, weights = dat_sim$weights)
  death_pols1 <- piy_analysis(yvars_d, info_d, pols1, exo_test = NULL, dat = dat_sim,
                              lags = 5, debias = TRUE, weights = dat_sim$weights)
  death_pols2 <- piy_analysis(yvars_d, info_d, pols1, exo_test = NULL, dat = dat_sim,
                              lags = 5, debias = TRUE, weights = dat_sim$weights)
  death_pols3 <- piy_analysis(yvars_d, info_d, pols1, exo_test = NULL, dat = dat_sim,
                              lags = 5, debias = TRUE, weights = dat_sim$weights)
  death_pols4 <- piy_analysis(yvars_d, info_d, pols1, exo_test = NULL, dat = dat_sim,
                              lags = 5, debias = TRUE, weights = dat_sim$weights)
  output <- c(case_pols1$fe, case_pols1$bc, case_pols2$fe, case_pols2$bc,
              case_pols3$fe, case_pols3$bc, case_pols4$fe, case_pols4$bc,
              death_pols1$fe, death_pols1$bc, death_pols2$fe, death_pols2$bc,
              death_pols3$fe, death_pols3$bc, death_pols4$fe, death_pols4$bc)
  return(output)
}

# pass to the cluster
ncores <- 28
btimes <- 500
set.seed(88, kind = "L'Ecuyer-CMRG")

result <- boot(data = sdf_week, statistic = estimator, sim = "parametric",
               ran.gen = data_wb, mle = 0, parallel = "multicore",
               ncpus = ncores, R = btimes)
result_wb <- structure(vapply(result$t, as.double, numeric(1)), dim = dim(result$t))
se_wb <- apply(result_wb, 2, function(x) {
  return((quantile(x, .75, na.rm = TRUE) - quantile(x, .25, na.rm = TRUE))/(qnorm(.75) - qnorm(.25)))
})

set.seed(88, kind = "L'Ecuyer-CMRG")

result2 <- boot(data = sdf_week, statistic = estimator, sim = "parametric",
                ran.gen = data_rg, mle = 0, parallel = "multicore",
                ncpus = ncores, R = btimes)
result_rg <- structure(vapply(result2$t, as.double, numeric(1)), dim = dim(result$t))
se_rg <- apply(result_rg, 2, function(x) {
  return((quantile(x, .75, na.rm = TRUE) - quantile(x, .25, na.rm = TRUE))/(qnorm(.75) - qnorm(.25)))
})


tab <- matrix(0, nrow = 164, ncol = 2)
colnames(tab) <- c("est", "se")
estimate <- c(results_pols1$fe, results_pols1$bc, results_pols2$fe, results_pols2$bc,
              results_pols3$fe, results_pols3$bc, results_pols4$fe, results_pols4$bc,
              results_d_pols1$fe, results_d_pols1$bc, results_d_pols2$fe, results_d_pols2$bc,
              results_d_pols3$fe, results_d_pols3$bc, results_d_pols4$fe, results_d_pols4$bc)
tab[, 1] <- estimate
rownames(tab) <- c(rownames(results_pols1$fe), rownames(results_pols1$bc),
                   rownames(results_pols2$fe), rownames(results_pols2$bc),
                   rownames(results_pols3$fe), rownames(results_pols3$bc),
                   rownames(results_pols4$fe), rownames(results_pols4$bc),
                   rownames(results_d_pols1$fe), rownames(results_d_pols1$bc),
                   rownames(results_d_pols2$fe), rownames(results_d_pols2$bc),
                   rownames(results_d_pols3$fe), rownames(results_d_pols3$bc),
                   rownames(results_d_pols4$fe), rownames(results_d_pols4$bc))
tab[, 2] <- se_wb

# now let's do the calibration
piy_cali <- function(yvar, dat, infovar, policyvar, exo_test,
                     cluster = "fips", fe = "fips + week",
                     lags = 2, weights = NULL) {
  # first, create formula
  iv <- "0" # this is for felm formula
  # lags
  if (lags != 0) policyvar <- sprintf("lag(%s, %d)", policyvar, lags)
  xvars <- c(policyvar, exo_test, infovar)
  rhs <- paste(xvars, collapse = " + ")
  rhs_felm <- paste(rhs, " | ", fe, " | ", iv, " | ", paste(cluster, collapse = " + "))
  fm <- as.formula(paste(yvar, "~", rhs_felm, sep = " "))
  # invoke felm to run the main regression
  piy <- felm(fm, dat, weights = weights, keepX = TRUE)
  # sum of TWFE and residuals
  index <- fitted.values(piy) - piy$X %*% as.matrix(coef(piy))
  # obtain the SE of the residuals
  sigma <- sqrt(sum(experiment$residuals^2)/experiment$df.residual)
  # obtain the coefficients
  coef0 <- coef(piy)
  results <- list(index = index, sigma = sigma, coef0 = coef0, X = piy$X)
  return(results)
}

# using the results from the piy_cali, generate the data for the case regression
data_generation <- function(results, dat0) {
  index <- results$index
  sigma <- results$sigma
  coef0 <- results$coef0
  n <- 3142 # number of counties
  T <- 32 # number of weeks
  y <- matrix(0, nrow = T, ncol = n)
  # simulate the dependent variable
  # retain the first four lags
  for (t in 1:4) y[t, ] <- dat0[(t + c(0:(n - 1))*T), "dlogdc"]
  # starting from t = 5, use the regression specification to generate data
  for (t in 5:T) {
    y[t, ] <- index[]
  }

}




