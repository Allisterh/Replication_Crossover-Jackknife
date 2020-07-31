# 7/12/2020 Shuowen Chen
# The script contains code that implement debiased fixed effects and panel bootstrap 
# standard errors for dataset in Chernozhukov, Kasahara and Schrimpf (2020)
# Serves as application for Chen, Chernozhukov, Fernandez-Val, Kasahara and Schripmf (2020)
# Acknowledgements: functions extend those written by Paul Schrimpf. 

library(lfe)
library(speedglm)
library(boot)
pathSC <- "/Users/shuowenchen/desktop/research/Panel-Cross-over-jackknife/application"
setwd(pathSC)
source('Code/auxiliaryfunctions.R')
# load data (sdf)
load("covidsince0307.RData")

# add weight to sdf
sdf$sweight <- 1
# create week factor to the data
sdf$week <- as.factor(strftime(sdf$date, format = "%V"))

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
  # plist <- list(pols, c("pmask.april","pmask.may", pols[-1]))
  # for pols only
   plist <- list(pols)
  
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
# cases as information
infovars <- list(c("dlogdc", "logdc"),
                 c("dlogdc", "logdc", "dlogdc.national", "logdc.national")) 
# deaths as information
infovarsd <- list(c("dlogdd", "logdd"),
                  c("dlogdd", "logdd", "dlogdd.national", "logdd.national")) 

test <- mainregressions_fe(sdf, yvar2, pols, bvars, infovarsd, NULL, L = 21, fixed = "state + month")
test2 <- mainregressions_fe(sdf, yvar2, pols, bvars, infovarsd, NULL, L = 21, fixed = "state + month")

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
num_boot <- 500 # number of bootstraps
ncores <- 1 # number of cpus (speed)

set.seed(88) 
# nonpar will throw warning wrt multiple splitting bias correction
case_nonpar <- boot(data = sdf, statistic = bootstat_case, sim = "parametric", 
                    ran.gen = data_rg, mle = 0, parallel = "multicore", 
                    ncpus = ncores, R = num_boot)
set.seed(88)
case_wb <- boot(data = sdf, statistic = bootstat_case, sim = "parametric", 
                ran.gen = data_wb, mle = 0, parallel = "multicore", 
                ncpus = ncores, R = num_boot)
set.seed(88)
deathdeath_nonpar <- boot(data = sdf, statistic = bootstat_deathdeath, 
                          sim = "parametric", ran.gen = data_rg, mle = 0, 
                          parallel = "multicore", ncpus = ncores, R = num_boot)
set.seed(88) 
deathdeath_wb <- boot(data = sdf, statistic = bootstat_deathdeath, 
                      sim = "parametric", ran.gen = data_wb, mle = 0, 
                      parallel = "multicore", ncpus = ncores, R = num_boot)

