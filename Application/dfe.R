# 11/10/2020 Shuowen Chen
# The script contains code that implement debiased fixed effects and panel bootstrap 
# standard errors for application in the crossover jackknife project

# load libraries, auxiliary functions and data
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

#### Empirical Analysis ####
# 1. Define lhs and rhs variables
# dependent variable for pbiy and piy
yvar <- "dlogdc" # case growth
yvar2 <- "dlogdd" # death growth

# policy variables
pols <- c("pmaskbus", "pk12", "pshelter", "pmovie", "prestaurant", "pnonessential") 

# behavioral variables
bvars <- c("workplaces", "retail", "grocery", "transit") 

# confounder from the SIR: change of test rates
tvars <- "dlogtests"

# information variables
# cases growth
infovars <- list(c("dlogdc", "logdc"),
                 c("dlogdc", "logdc", "dlogdc.national", "logdc.national")) 
# deaths growth
infovarsd <- list(c("dlogdd", "logdd"),
                  c("dlogdd", "logdd", "dlogdd.national", "logdd.national")) 

# run regression
results_case <- mainregressions_fe(sdf, yvar, pols, bvars, infovarsd, tvars, L = 14, fixed = "state + week")
results_death <- mainregressions_fe(sdf, yvar2, pols, bvars, infovarsd, NULL, L = 21, fixed = "state + week")

#### 2. Bootstrap Standard Errors 
# Compute estimates in each bootstrap
bootstat_case <- function(data) {
  # change the specification here
  tot <- mainregressions_fe(data, yvar, pols, bvars, infovars, tvars, fixed = "state + week")
  return(unlist(tot))
}

bootstat_deathdeath <- function(data) {
  tot <- mainregressions_fe(data, yvar2, pols, bvars, infovarsd, NULL, L = 21, fixed = "state + week")
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

