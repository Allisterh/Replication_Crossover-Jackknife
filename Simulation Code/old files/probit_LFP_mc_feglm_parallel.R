####################################################################
# MONTE CARLO CALIBRATED TO LFP:  DYNAMIC PROBIT MODEL USING FEGLM
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
###################################################################

### set working directory
filepath <- "/Users/shuowenchen/desktop/Research/Panel-Cross-over-jackknife/Crossover-Jackknife/Simulation Code"
setwd(filepath)

# Load functions to estimate probit fast and perform bias corrections;

library(foreign)
library(xtable)
library(alpaca)
library(boot)

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
B <- 5 # number of bootstrap
N <- N0 # number of individuals
T <- T0 # number of time periods

# Regression specifications
spec_fe <- lfp ~ laglfp + kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age+ age2 | id + year
# specification for simulated data (without intercept)
# do not include intercept because do not have it in simulation in the first place
spec_fe_s <- lfp ~ laglfp + kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age + age2 | id + year
# Check number of missing LFP values
nomissing <- !is.na(data$lfp)
pop <- length(data$lfp[nomissing])

# TRIMMING PARAMETERS FOR ANALYTICAL CORRECTIONS
L2 <- 1
L3 <- 0
L4 <- 0

# CALIBRATION OF PARAMETERS
deltahatfete <- feglm(spec_fe, data = data, family = binomial(link = "probit"))
betahatfete <- deltahatfete$coefficients
index <- predict(deltahatfete)
index0 <- index - betahatfete[1] * data$laglfp
index1 <- index + betahatfete[1] * (1 - data$laglfp)

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

# (1) Data Simulation
data_gen <- function(data, mle) {
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
  # All other covariates are taken from the original data
  data_s <- data.frame(year = data$year, lfp = lfp_s, laglfp = laglfp_s, kids0_2 = data$kids0_2, 
                       kids3_5 = data$kids3_5, kids6_17 = data$kids6_17, 
                       loghusbandincome = data$loghusbandincome, age = data$age, age2 = data$age2, id = data$id)
  return(data_s)
}
# (2) Bootstrap statistics function
boot_stat <- function(data, spec_fe_s) {
  # ELIMINATING INDIVIDUALS WITH ALL 1's OR 0's
  D <- kronecker(c(1:N0), matrix(1, T0, 1))
  # tapply(lfp_s, D, sum): counts LFP for each individual over the spell of the years (0-9)
  td <- kronecker(tapply(data$lfp, D, sum),  matrix(1, T0, 1)) 
  insample <- ((td > 0) & (td < T0))
  N <- sum(insample)/T0 # population of this simulated data
  # SELECT SUBSAMPLE;
  innomissing <- nomissing[insample]
  ### Estimation using Simulated Data ###
  
  # FIXED EFFECTS on the simulated data;
  deltahatfete <- feglm(spec_fe_s, data = data[insample, ], family = binomial(link = "probit"))
  betahatfete <- deltahatfete$coefficients
  
  # ANALYTICAL BIAS CORRECTION;
  betahatabc <- t(betahatfete - bias_dynamic_parameter(N, T0, predict(deltahatfete), 
                                                       as.matrix(deltahatfete$data[, 2:(P + 1)]), 
                                                       data$lfp[insample], L2 = L2, L3 = L3, L4 = L4))               
  # SPLIT-PANEL JACKKNIFE IN BOTH DIMENSIONS;
  # LEAVING HALF PANEL OUT, ALL TIME PERIODS
  # i = 1, ..., N/2, j = 1, ..., T
  index <- (data$id <= median(data$id[nomissing & insample]))
  deltahat11 <- feglm(spec_fe_s, data = data[index, ], family = binomial(link = "probit"))
  betahat11 <- deltahat11$coefficients
  # i = N/2 + 1, ..., N, t = 1, ..., T
  index <- (data$id >= median(data$id[nomissing & insample]))  
  deltahat12 <- feglm(spec_fe_s, data = data[index, ], family = binomial(link = "probit"))
  betahat12 <- deltahat12$coefficients
  
  # LEAVING HALF PANEL OUT, ALL INDIVIDUALS
  # i = 1, ..., N, j = 1, ..., T/2
  index <- (data$year <= median(data$year[nomissing & insample]))  
  deltahat21 <- feglm(spec_fe_s, data = data[index, ], family = binomial(link = "probit"))
  betahat21 <- deltahat21$coefficients
  # i = 1, ..., N, j = T/2 + 1, ..., T
  index <- (data$year >= median(data$year[nomissing & insample]))  
  deltahat22 <- feglm(spec_fe_s, data = data[index, ], family = binomial(link = "probit"))
  betahat22 <- deltahat22$coefficients
  # SBC (Fernandez-Val and Weidner, 2016)
  betahatjbc <- 3*betahatfete - (betahat11 + betahat12 + betahat21 + betahat22)/2    
  
  # CROSS-OVER JACKKNIFE;
  # FIRST CROSS: i = 1, ..., N/2, j = 1, ..., T/2; i = N/2 + 1, ..., N, j = T/2 + 1, ..., T
  index <- (data$id <= median(data$id[nomissing & insample]) & data$year <= median(data$year[nomissing & insample])) | 
    (data$id >= median(data$id[nomissing & insample]) & data$year >= median(data$year[nomissing & insample]))
  deltahat1 <- feglm(spec_fe_s, data = data[index, ], family = binomial(link = "probit"))
  betahat1 <- deltahat1$coefficients
  # SECOND CROSS: i = 1, ..., N/2,  j = T/2 + 1, ..., T; i = N/2 + 1, ..., N, j = 1, ..., T/2
  index <- (data$id <= median(data$id[nomissing & insample]) & data$year >= median(data$year[nomissing & insample])) | 
    (data$id >= median(data$id[nomissing & insample]) & data$year <= median(data$year[nomissing & insample]))
  deltahat2 <- feglm(spec_fe_s, data = data[index, ], family = binomial(link = "probit"))
  betahat2 <- deltahat2$coefficients
  # CBC
  betahatcbc <- 2*betahatfete - (betahat1 + betahat2)/2    
  # STANDARD ERRORS;
  se_betahatfete <- summary(deltahatfete)[1]$cm[, 2]
  # return output
  return(c(betahatfete, betahatabc, betahatjbc, betahatcbc, se_betahatfete))
}
set.seed(88)
sim_fe <- boot(data = data, statistic = boot_stat, sim = "parametric", ran.gen = data_gen, mle = 0, 
               spec_fe_s = spec_fe_s, parallel = "multicore", ncpus = 1, R = B)

# Post Simulation Summary
sim_tot <- matrix(as.numeric(sim_fe$t), nrow = B, ncol = 35)
sim_nobc <- sim_tot[, 1:7]
sim_abc <- sim_tot[, 8:14]
sim_sbc <- sim_tot[, 15:21]
sim_cbc <- sim_tot[, 22:28]
sim_se <- sim_tot[, 29:35]

sumstats <- function(sim, real, se_sim, conv) {
  output <- matrix(0, nrow = dim(sim)[2], ncol = 5)
  colnames(output) <- c("bias", "std", "RMSE", "SE/SD", "p95")
  rownames(output) <- names(real)
  for (i in 1:dim(sim)[2]) {
    # bias
    output[i, 1] <- 100*(mean(sim[conv, i], na.rm = TRUE)/real[i] - 1)
    # standard deviation
    output[i, 2] <- 100*sd(sim[conv, i]/real[i], na.rm = TRUE)
    # rmse
    output[i, 3] <- 100*sqrt(mean((sim[conv, i]/real[i] - 1)^2, na.rm = TRUE))
    # se_sd
    output[i, 4] <- (mean((se_sim[conv, i]/sd(sim[conv, i], na.rm = TRUE)), na.rm = TRUE))
    # pvalue
    output[i, 5] <- (mean((sim[conv, i] + qnorm(.05/2) * se_sim[conv, i] <= real[i]) &
                            (sim[conv, i] + qnorm(1 - .05/2) * se_sim[conv, i] >= real[i]), na.rm = TRUE))
  }
  return(output)
}

# covergence 
convergence <- abs(sim_sbc[, 1]) < 100 & abs(sim_sbc[, 2]) < 100 & abs(sim_cbc[, 1]) < 100 & abs(sim_cbc[, 2]) < 100

sum_nobc <- sumstats(sim_nobc, beta0, sim_se, convergence)
sum_abc <- sumstats(sim_abc, beta0, sim_se, convergence)
sum_sbc <- sumstats(sim_sbc, beta0, sim_se, convergence)
sum_cbc <- sumstats(sim_cbc, beta0, sim_se, convergence)
