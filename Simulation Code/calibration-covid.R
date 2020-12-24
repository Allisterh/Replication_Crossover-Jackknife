# 12/23/2020 Shuowen Chen
# This script conducts calibrated exercise using Covid data

library(lfe)
library(speedglm)
library(boot)

filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
setwd(filepath)
# load functions for regression, bias correction and bootstrap
source('auxiliaryfunctions.R')
# load data (sdf)
load("covidsince0307.RData")

# add weight to sdf
sdf$sweight <- 1
# create week factor to the data
sdf$week <- as.factor(strftime(sdf$date, format = "%V"))

#### Empirical Analysis ####
#### 1. Defining lhs and rhs variables
# dependent variable for pbiy and piy
yvar <- "dlogdc" # cases growth
yvar2 <- "dlogdd" # deaths growth
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

# regression specification
casereg <- mainregressions(sdf, yvar, pols, bvars, infovars, tvars, L = 14, 
                           fixed = "state + week", spec = "piy", purpose = "cali")

deathreg <- mainregressions(sdf, yvar2, pols, bvars, infovarsd, NULL, L = 21, 
                            fixed = "state + week", spec = "piy", purpose = "cali")




