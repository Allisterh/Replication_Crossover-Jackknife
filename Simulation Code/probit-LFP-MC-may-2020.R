####################################################################
# MONTE CARLO CALIBRATED TO LFP:  DYNAMIC PROBIT MODEL
####################################################################
#  14.382 L10 MIT.  V. Chernozhukov and I. Fernandez-Val
#  and Shuowen Chen
# Data source: Ivan Fernandez-Val, "Fixed effects estimation of structural parameters and marginal effects in 
# panel probit models," Journal of Econometrics, Elsevier, vol. 150(1), pages 71-85, May.
  
# The variables in the data set include:
#
# id = "woman identifier"
# year = "year"
# lfp = "= 1 if woman participates in labor force"
# laglfp = "lag of lfp"
# kids0_2 = "number of children of ages <= 2 years"
# kids3_5 = "number of children of ages > 2 years and < 6 years"
# kids6_17 = "number of children of ages > 5 years and < 18 years"
# loghusbandincome = "log of husband income ($1000 of 1995)"
# age = "age"
# age2 = "age squared"
####################################################################
  
  
#rm(list = ls());

### set working directory
filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
setwd(filepath)

# Load functions to estimate probit fast and perform bias corrections;

library(foreign)
library(speedglm)
library(xtable)
library(boot)
library(doParallel)
library(doRNG)

source('probit-functions.R')
# Contains the following three functions: 
# (1) bias_dynamic_parameter, (2) bias_parameter and (3) spglm

# Reading the data;
data  <- read.dta("LFPmovers.dta")

# VALUES OF THE PARAMETERS;
N0 <- length(unique(data$id))
T0 <- length(unique(data$year))
P0 <- 7 # number of nonfe regressors in the specification
P <- 7
B <- 2 # number of bootstrap
N <- N0 # number of individuals
T <- T0 # number of time periods

# Regression specifications
spec_fe <- lfp ~ laglfp + kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age+ age2 + factor(id) + factor(year)
# specification for simulated data (without intercept)
# do not include intercept because do not have it in simulation in the first place
spec_fe_s <- lfp_s ~ laglfp_s + kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age + age2 + factor(id) + factor(year) - 1

# Check number of missing LFP values
nomissing <- !is.na(data$lfp)
pop <- length(data$lfp[nomissing])

# TRIMMING PARAMETERS FOR ANALYTICAL CORRECTIONS
L2 <- 1
L3 <- 0
L4 <- 0

# CALIBRATION OF PARAMETERS
deltahatfete <- spglm(spec_fe, data = data, family = binomial(link = "probit"), subset = nomissing)
betahatfete <- deltahatfete$coef[2:(P0 + 1)] # subset the coefficients of nonfe regressors

# predicted values within the link function
index <- deltahatfete$linear.predictors 
# X'*theta
index0 <- index - betahatfete[1] * data$laglfp[nomissing] 
# X'*theta + coef of laglfp
index1 <- index + betahatfete[1] * (1 - data$laglfp[nomissing]) 

# Estimated coefficients from the data
beta0 <- betahatfete[c(1:P)]

true_index0 <- matrix(index0, nrow = T, byrow = FALSE)
true_index1 <- matrix(index1, nrow = T, byrow = FALSE)

rm(deltahatfete, index, index0, index1)

# VARIABLE INITIALIZIN

betasfete <- matrix(0, nrow = B, ncol = P)
betasabc <- matrix(0, nrow = B, ncol = P)
betasjbc <- matrix(0, nrow = B, ncol = P)
betascbc <- matrix(0, nrow = B, ncol = P)
se_betasfete <- matrix(0, nrow = B, ncol = P)

# SIMULATION
set.seed(88)
D <- kronecker(c(1:N0), matrix(1, T0, 1))

lfpstore <- matrix(0, nrow = N0*T0, ncol = B)
# loop over bootstrap simulations
for (i in 1:B) {
  # CREATE OUTCOME VARIABLE
  lfp_s <- matrix(0, nrow = T0, ncol = N0)
  lfp0 <- data$laglfp[data$year == 1] #initial value
  for (t in 1:T0) {
    # Interpretation: lfp0*betahatfete[1] + X'*betahatfete[2:end] > rnorm(N0),
    # Determines the lfp in the next period
    lfp_s[t, ] <- (lfp0 * true_index1[t, ] + (1 - lfp0) * true_index0[t, ] > rnorm(N0))
    lfp0 <- lfp_s[t, ] # update so as to simulate the next period lfp
  }
  # simulate lagged dependent variable
  laglfp_s <- matrix(rbind(data$laglfp[data$year == 1], lfp_s[-T0, ]), ncol = 1, byrow = FALSE)
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  lfpstore[, i] <- lfp_s
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp_s = lfp_s, laglfp_s = laglfp_s, kids0_2 = data$kids0_2, kids3_5 = data$kids3_5, 
                       kids6_17 = data$kids6_17, loghusbandincome = data$loghusbandincome, age = data$age, 
                       age2 = data$age2, id = data$id, year = data$year)
  # ELIMINATING INDIVIDUALS WITH ALL 1's OR 0's
  # tapply(lfp_s, D, sum): counts LFP for each individual over the spell of the years (0-9)
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T0, 1)) 
  insample <- ((td > 0) & (td < T0))
  N <- sum(insample)/T0 # population of this simulated data
  # SELECT SUBSAMPLE;
  innomissing <- nomissing[insample]
  
  ### Estimation using Simulated Data ###
  # FIXED EFFECTS on the simulated data;
  deltahatfete <- spglm(spec_fe_s, data = data_s[insample, ], family = binomial(link = "probit"), x = TRUE)
  betahatfete <- deltahatfete$coef[1:P]
  
  # ANALYTICAL BIAS CORRECTION;
  betahatabc <- t(betahatfete - bias_dynamic_parameter(N, T0, deltahatfete$linear.predictors, 
                                                        deltahatfete$x[, 1:P], lfp_s[insample], 
                                                        L2 = L2, L3 = L3, L4 = L4))               

  # SPLIT-PANEL JACKKNIFE IN BOTH DIMENSIONS;
  
  # LEAVING HALF PANEL OUT, ALL TIME PERIODS
  # i = 1, ..., N/2, j = 1, ..., T
  index <- (data_s$id <= median(data_s$id[nomissing & insample]))
  deltahat11 <- spglm(spec_fe_s, data = data_s[index, ], family = binomial(link = "probit"))
  betahat11 <- deltahat11$coef[1:P]
  # i = N/2 + 1, ..., N, t = 1, ..., T
  index <- (data_s$id >= median(data_s$id[nomissing & insample]))  
  deltahat12 <- spglm(spec_fe_s, data = data_s[index, ], family = binomial(link = "probit"))
  betahat12 <- deltahat12$coef[1:P]
  
  # LEAVING HALF PANEL OUT, ALL INDIVIDUALS
  # i = 1, ..., N, j = 1, ..., T/2
  index <- (data_s$year <= median(data_s$year[nomissing & insample]))  
  deltahat21 <- spglm(spec_fe_s, data = data_s[index, ], family = binomial(link = "probit"))
  betahat21 <- deltahat21$coef[1:P]
  # i = 1, ..., N, j = T/2 + 1, ..., T
  index <- (data_s$year >= median(data_s$year[nomissing & insample]))  
  deltahat22 <- spglm(spec_fe_s, data = data_s[index, ], family = binomial(link = "probit"))
  betahat22 <- deltahat22$coef[1:P]
  # SBC (Fernandez-Val and Weidner, 2016)
  betahatjbc <- 3*betahatfete - (betahat11 + betahat12 + betahat21 + betahat22)/2    
  
# CROSS-OVER JACKKNIFE;
  # FIRST CROSS: i = 1, ..., N/2, j = 1, ..., T/2; i = N/2 + 1, ..., N, j = T/2 + 1, ..., T
  index <- (data_s$id <= median(data_s$id[nomissing & insample]) & data_s$year <= median(data_s$year[nomissing & insample])) | 
    (data_s$id >= median(data_s$id[nomissing & insample]) & data_s$year >= median(data_s$year[nomissing & insample]))
  deltahat1 <- spglm(spec_fe_s, data = data_s[index, ], family = binomial(link = "probit"))
  betahat1 <- deltahat1$coef[1:P]
  # SECOND CROSS: i = 1, ..., N/2,  j = T/2 + 1, ..., T; i = N/2 + 1, ..., N, j = 1, ..., T/2
  index <- (data_s$id <= median(data_s$id[nomissing & insample]) & data_s$year >= median(data_s$year[nomissing & insample])) | 
    (data_s$id >= median(data_s$id[nomissing & insample]) & data_s$year <= median(data_s$year[nomissing & insample]))
  deltahat2 <- spglm(spec_fe_s,  data = data_s[index, ], family = binomial(link = "probit"))
  betahat2 <- deltahat2$coef[1:P]
  # CBC
  betahatcbc <- 2*betahatfete - (betahat1 + betahat2)/2    
  # STANDARD ERRORS;
  se_betahatfete <- summary(deltahatfete)$coeff[1:P, 2]
  
  # SAVE ITERATION RESULTS
  betasfete[i, ] <- betahatfete
  betasabc[i, ] <- betahatabc
  betasjbc[i, ] <- betahatjbc
  betascbc[i, ] <- betahatcbc
  se_betasfete[i, ] <- se_betahatfete

  # REPORT RESULTS OF EACH ITERATION FOR CONTROL; 
  if (i > 1) {    
  print(i)  
  print(rbind(beta0, apply(betasfete[1:i,], 2, mean), apply(betasabc[1:i,], 2, mean), 
              apply(betasjbc[1:i,], 2, mean), apply(betascbc[1:i,], 2, mean)))
  print(rbind(apply(betasfete[1:i,], 2, sd), apply(betasabc[1:i, ], 2, sd), 
              apply(betasjbc[1:i, ], 2, sd), apply(betascbc[1:i, ], 2, sd), 
              apply(se_betasfete[1:i, ], 2, mean)))
  }
}

# CHECK CONVERGENCE;
conver <- abs(betasjbc[, 1]) < 100 & abs(betasjbc[, 2]) < 100 & abs(betascbc[, 1]) < 100 & abs(betascbc[, 2]) < 100
print(B - sum(conver))

for (ind in 1:P) {
  table.beta <- matrix(0, ncol = 5, nrow = 4, dimnames = list(c('MLE-FETE','Analytical','Jackknife','Cross-Over'), 
                                                               c('Bias', 'Std. Dev.', 'RMSE', 'SE/SD', 'p:.95')))
  
  table.beta[1,1] <- 100*mean(betasfete[conver, ind]/beta0[ind] - 1)
  table.beta[1,2]	<- 100*sd(betasfete[conver, ind]/beta0[ind])
  table.beta[1,3]	<- 100*sqrt(mean((betasfete[conver, ind]/beta0[ind] - 1)^2))
  table.beta[1,4] <- (mean((se_betasfete[conver, ind]/sd(betasfete[conver, ind]))))
  table.beta[1,5] <- (mean((betasfete[conver, ind] + qnorm(.05/2) * se_betasfete[conver, ind] <= beta0[ind]) &
                             (betasfete[conver, ind] + qnorm(1 - .05/2) * se_betasfete[conver, ind] >= beta0[ind])))
  
  table.beta[2,1] <- 100*mean(betasabc[conver, ind]/beta0[ind] - 1)
  table.beta[2,2]	<- 100*sd(betasabc[conver, ind]/beta0[ind])
  table.beta[2,3]	<- 100*sqrt(mean((betasabc[conver, ind]/beta0[ind] - 1)^2))
  table.beta[2,4] <- (mean((se_betasfete[conver, ind]/sd(betasabc[conver, ind]))))
  table.beta[2,5] <- (mean((betasabc[conver, ind] + qnorm(.05/2) * se_betasfete[conver, ind] <= beta0[ind]) & 
                             (betasabc[conver, ind] + qnorm(1 - .05/2) * se_betasfete[conver, ind] >= beta0[ind])))
  
  table.beta[3,1] <- 100*mean(betasjbc[conver, ind]/beta0[ind] - 1)
  table.beta[3,2] <- 100*sd(betasjbc[conver, ind]/beta0[ind])
  table.beta[3,3]	<- 100*sqrt(mean((betasjbc[conver, ind]/beta0[ind] - 1)^2))
  table.beta[3,4] <- (mean((se_betasfete[conver, ind]/sd(betasjbc[conver, ind]))))
  table.beta[3,5] <- (mean((betasjbc[conver, ind] + qnorm(.05/2) * se_betasfete[conver, ind] <= beta0[ind]) & 
                             (betasjbc[conver, ind] + qnorm(1 - .05/2) * se_betasfete[conver, ind] >= beta0[ind])))
  
  table.beta[4,1] <- 100*mean(betascbc[conver, ind]/beta0[ind] - 1)
  table.beta[4,2] <- 100*sd(betascbc[conver, ind]/beta0[ind])
  table.beta[4,3]	<- 100*sqrt(mean((betascbc[conver, ind]/beta0[ind] - 1)^2))
  table.beta[4,4] <- (mean((se_betasfete[conver, ind]/sd(betascbc[conver, ind]))))
  table.beta[4,5] <- (mean((betascbc[conver, ind] + qnorm(.05/2) * se_betasfete[conver, ind] <= beta0[ind]) & 
                             (betascbc[conver, ind] + qnorm(1 - .05/2) * se_betasfete[conver, ind] >= beta0[ind])))
  
  print(xtable(table.beta, digits = 2))
}


