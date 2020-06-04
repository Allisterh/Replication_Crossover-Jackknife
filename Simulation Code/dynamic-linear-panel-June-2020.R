####################################################################
# MONTE CARLO CALIBRATED TO DYNAMIC LINEAR PANEL MODEL
####################################################################
# Cross-Over Jackknife of Nonlinear Panel
# Shuowen Chen, Victor Chernozhukov and Ivan Fernandez-Val
# Data source: Daron Acemoglu (MIT), N = 147 countries, T = 23 years (1987-2009)
# include data for four lags of dependent variable, balanced panel
#
# Description of the data: the sample selection and variable contruction follow ANRR
#
# The variables in the data set include:
#
# country_name    = Country name
# wbcode          = World Bank country code
# year            = Year 
# id              = Generated numeric country code
# dem             = Democracy measure by ANRR
# lgdp            = log of GDP per capita in 2000 USD from World Bank
####################################################################

library(foreign)
library(xtable)
library(plm)
library(gmm)

data2 <- read.dta("democracy-balanced-l4.dta")
data2 <- pdata.frame(data2, index = c("id", "year"))

################################
## Calibration by Fixed Effects
################################
# Regression Specification
form <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id)
# Create lagged dependent variable
data2$l1lgdp <- lag(data2$lgdp, 1)
data2$l2lgdp <- lag(data2$lgdp, 2) 
data2$l3lgdp <- lag(data2$lgdp, 3) 
data2$l4lgdp <- lag(data2$lgdp, 4) 

##### 1. FE using data at hand ####
fit <- lm(form, data2) 
coefs <- coef(fit)[2:6]
# Long run effects of treatment
lr <- coefs[1]/(1 - sum(coefs[2:5]))

##### 2. Compute X*beta_hat + time FE + ind FE (will be used in the data simulation) ####
# Drop the first four years at they contain NA due to lagging
speriod <- !((data2$year == 1987) | (data2$year == 1988) | (data2$year == 1989) | (data2$year == 1990))
sdata <- data2[speriod, ]
sdata$l1lgdp <- 0
sdata$l2lgdp <- 0
sdata$l3lgdp <- 0
sdata$l4lgdp <- 0
index0 <- predict(fit, sdata)
rm(sdata)

## Auxiliary Functions 
# (a) Generate data sets
data_rg <- function(data, fit_original = fit, index0 = index0, mle) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  coefs <- coef(fit_original)[3:6] # coefs on lagged dependent variable
  sigma <- summary(fit_original)$sigma 
  slgdp <- matrix(data$lgdp, nrow = T, ncol = N)
  index0 <- matrix(index0, nrow = T - 4, ncol = N)
  # Generate simulated dependent variable from t=5 
  for (t in 5:T) {
    slgdp[t, ] <- index0[t - 4, ] + coefs[1]*slgdp[t - 1, ] + 
      coefs[2]*slgdp[t - 2, ] + coefs[3]*slgdp[t - 3, ] + 
      coefs[4]*slgdp[t - 4, ] + sigma*rnorm(N)
  }
  data_b <- data
  # Replace with simulated dependent variable
  data_b$lgdp <- matrix(slgdp, ncol = 1)
  data_b <- data.frame(data_b)
  data_b <- pdata.frame(data_b, index = c("id", "year"))    
  return(data_b)
}
# (b) Analytical Bias correction
abc <- function(data, form, lags, N) {
  data$l1lgdp <- lag(data$lgdp, 1)
  data$l2lgdp <- lag(data$lgdp, 2) 
  data$l3lgdp <- lag(data$lgdp, 3) 
  data$l4lgdp <- lag(data$lgdp, 4) 
  fit <- lm(form, data, x = TRUE, na.action = na.omit)
  res <- fit$residuals
  jac <- solve(t(fit$x) %*% fit$x / length(res))[2:6, 2:6]
  indexes <- c(1:length(res))
  bscore <- rep(0, 5)
  T <- length(res)/N
  for (i in 1:lags) {
    indexes <- indexes[-c(1 + c(0:(N - 1))*T)]
    lindexes <- indexes - i
    bscore <- bscore + t(fit$x[indexes, 2:6]) %*% res[lindexes] / length(indexes)
  }
  bias <- -jac %*% bscore
  return(as.vector(bias/T))
}

# (c) Split Sample Jackknife
sbc <- function(data, form) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  # Full sample estimation
  whole <- plm(form, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Split across time series dimension
  t_lower <- plm(form, data, subset = (as.double(data$year) <= ceiling(T/2)), 
                 model = "within", effect = "twoways", index = c("id", "year"))
  t_upper <- plm(form, data, subset = (as.double(data$year) >= floor(1 + T/2)), 
                 model = "within", effect = "twoways", index = c("id", "year"))
  bc_t <- 0.5*(coef(t_lower) + coef(t_upper))
  # Split across individual dimension
  n_lower <- plm(form, data, subset = (as.double(data$id) <= ceiling(N/2)), 
                 model = "within", effect = "twoways", index = c("id", "year"))
  n_upper <- plm(form, data, subset = (as.double(data$id) >= floor(1 + N/2)), 
                 model = "within", effect = "twoways", index = c("id", "year"))
  bc_n <- 0.5*(coef(n_lower) + coef(n_upper))
  # SBC bias correction
  coef_sbc <- 3*coef(whole) - bc_t - bc_n
  return(coef_sbc)
}

# (d) Cross-over Bias Correction
cbc <- function(data, form) {
  N <- length(unique(data$id))
  T <- length(unique(data$year))
  # Full sample estimation
  whole <- plm(form, data, model = "within", effect = "twoways", index = c("id", "year"))
  # Two cross-over subsamples
  subsample1 <- as.double(data$id) <= 
  subsample2 
}
#### 3. Simulation ####
R <- 500

## FIXED EFFECTS APPROACH 
## statistics to be computed in each bootstrap draw ##
sim_fe <- function(data, form_fe, form_abc){
  # Fixed Effects
  fe_fit <- plm(form_fe, data, model = "within", effect = "twoways", index = c("id", "year"))
  coefs_fe <- coef(fe_fit)
  lr_fe <- coefs_fe[1]/(1 - sum(coefs_fe[2:5]))
  # jackknife bias correction
  fe_fit1 <- plm(form_fe, data, subset = (as.double(year) <= 14), model = "within", 
                 effect = "twoways", index = c("id", "year"))
  fe_fit2 <- plm(form_fe, data, subset = (as.double(year) >= 14), model = "within", 
                 effect = "twoways", index = c("id", "year"))
  coefs_jbc  <- 19*coef(fe_fit)/9 - 10*(coef(fe_fit1) + coef(fe_fit2))/18
  lr_fe1     <- coef(fe_fit1)[1]/(1 - sum(coef(fe_fit1)[2:5]))
  lr_fe2     <- coef(fe_fit2)[1]/(1 - sum(coef(fe_fit2)[2:5]))
  lr_jbc     <- 19*lr_fe/9 - 10*(lr_fe1 + lr_fe2)/18
  
  # Analytical bias correction
  bias_l4 <- abc(data, form_abc, lags = 4, N = length(levels(id)))
  coefs_abc4 <- coefs_fe - bias_l4
  jac_lr <- c(1,rep(lr_fe, 4))/(1 - sum(coefs_fe[2:5]))
  lr_abc4 <- lr_fe - crossprod(jac_lr, bias_l4)
  
  return(c(coefs_fe, coefs_jbc, coefs_abc4, lr_fe, lr_jbc, lr_abc4))
}

library(boot) # library to do bootstrap with paralell computing

set.seed(888)
form_fe <- lgdp ~ dem + lag(lgdp, 1:4) - 1
form_abc <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id)
result.sim.fe <- boot(data = data, statistic = sim_fe, sim = "parametric", ran.gen = data_rg, 
                      mle = 0, form_fe = form_fe, form_abc = form_abc,  parallel = "multicore", 
                      ncpus = 10, R = R)
# robust estimator of std deviation based on IQR
rsd <- function(x) {return((quantile(x, .75, na.rm = TRUE) - quantile(x, .25, na.rm = TRUE))/(qnorm(.75) - qnorm(.25)))} 

# Robust bootstrap std errors;
result      <- structure(vapply(result.boot.SE.fe$t, as.double, numeric(1)), dim=dim(result.boot.SE.fe$t)); # transforms "Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'\n" to NA
bse.fe      <- apply(result[,1:5], 2, rsd);
bse.jbc     <- apply(result[,6:10], 2, rsd);
bse.abc4    <- apply(result[,11:15], 2, rsd);
bse.lr.fe   <- rsd(result[,16]);
bse.lr.jbc  <- rsd(result[,17]);
bse.lr.abc4 <- rsd(result[,18]);

## ARELLANO-BOND APPROACH #######################################################################
########## statistics to be computed in each bootstrap draw #####################################
boot.SE.ab <- function(data, form.ab){
  
  # Arellano-Bond
  N         <- length(unique(data$id));
  ab.fit <- pgmm(form.ab, data, model = "onestep", effect = "twoways" );
  coefs.ab  <- coef(ab.fit); 
  lr.ab     <- coefs.ab[1]/(1-sum(coefs.ab[2:5]));
  
  # Split sample bias correction - 1 partition
  S2        <- 1;
  acoeff.ab <- 0*coef(ab.fit);
  alr.ab    <- 0;
  for (s in 1:S2) {
    sample1   <- sample(N, ceiling(N/2), replace = FALSE);
    ab.fit1 <- pgmm(form.ab, data[as.double(id) %in% sample1, ], model = "onestep", effect = "twoways" );
    ab.fit2 <- pgmm(form.ab, data[!(as.double(id) %in% sample1), ], model = "onestep", effect = "twoways" );
    lr.ab1     <- coef(ab.fit1)[1]/(1-sum(coef(ab.fit1)[2:5]));
    lr.ab2     <- coef(ab.fit2)[1]/(1-sum(coef(ab.fit2)[2:5]));
    acoeff.ab <- acoeff.ab + ((coef(ab.fit1) + coef(ab.fit2))/2)/S2;
    alr.ab    <- alr.ab + ((lr.ab1 + lr.ab2)/2)/S2;
  }
  coefs.ab.jbc  <- 2*coef(ab.fit) - acoeff.ab;
  lr.ab.jbc     <- 2*lr.ab - alr.ab;
  
  # Split sample bias correction - 5 partitions
  S2        <- 5;
  acoeff.ab <- 0*coef(ab.fit);
  alr.ab    <- 0;
  for (s in 1:S2) {
    sample1   <- sample(N, ceiling(N/2), replace = FALSE);
    ab.fit1 <- pgmm(form.ab, data[as.double(id) %in% sample1, ], model = "onestep", effect = "twoways" );
    ab.fit2 <- pgmm(form.ab, data[!(as.double(id) %in% sample1), ], model = "onestep", effect = "twoways" );
    lr.ab1     <- coef(ab.fit1)[1]/(1-sum(coef(ab.fit1)[2:5]));
    lr.ab2     <- coef(ab.fit2)[1]/(1-sum(coef(ab.fit2)[2:5]));
    acoeff.ab <- acoeff.ab + ((coef(ab.fit1) + coef(ab.fit2))/2)/S2;
    alr.ab    <- alr.ab + ((lr.ab1 + lr.ab2)/2)/S2;
  }
  coefs.ab.jbc5  <- 2*coef(ab.fit) - acoeff.ab;
  lr.ab.jbc5     <- 2*lr.ab - alr.ab;
  
  
  return(c(coefs.ab[1:5], coefs.ab.jbc[1:5], coefs.ab.jbc5[1:5], lr.ab, lr.ab.jbc, lr.ab.jbc5));
}

library(boot); # library to do bootstrap with paralell computing

set.seed(888);
result.boot.SE.ab <- boot(data = data, statistic=boot.SE.ab, sim = "parametric", ran.gen = data.rg, mle = 0, form.ab = form.ab, 
                          parallel="multicore", ncpus = 20, R=R);

# Robust bootstrap std errors;

result      <- structure(vapply(result.boot.SE.ab$t, as.double, numeric(1)), dim=dim(result.boot.SE.ab$t)); # transforms "Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'\n" to NA
bse.ab      <- apply(result[,1:5], 2, rsd);
bse.ab.jbc  <- apply(result[,6:10], 2, rsd);
bse.ab.jbc5 <- apply(result[,11:15], 2, rsd);
bse.lr.ab   <- rsd(result[,16]);
bse.lr.ab.jbc  <- rsd(result[,17]);
bse.lr.ab.jbc5  <- rsd(result[,18]);



######## Table of results;

options(digits=2);
table.all <- matrix(NA, nrow = 18, ncol = 6, dimnames = list(c("Democracy", "CSE", "BSE", "L1.log(gdp)",  "CSE1", "BSE1", "L2.log(gdp)",  "CSE2", "BSE2","L3.log(gdp)",  "CSE3", "BSE3", "L4.log(gdp)",  "CSE4", "BSE4", "LR-Democracy","CSE5","BSE5"), c("FE", "SBC", "ABC4", "AB", "SBC1", "SBC5" )));


table.all[c(1,4,7,10,13), 1] <- coefs.fe[1:5];
table.all[c(2,5,8,11,14), 1] <- cse.fe[1:5];
table.all[c(3,6,9,12,15), 1] <- bse.fe[1:5];

table.all[c(1,4,7,10,13), 2] <- coefs.jbc[1:5];
#table.all[c(2,5,8,11,14), 2] <- cse.fe[1:5];
table.all[c(3,6,9,12,15), 2] <- bse.jbc[1:5];

table.all[c(1,4,7,10,13), 3] <- coefs.abc4[1:5]
table.all[c(3,6,9,12,15), 3] <- bse.abc4[1:5]

table.all[c(1,4,7,10,13), 4] <- coefs.ab[1:5]
table.all[c(2,5,8,11,14), 4] <- cse.ab[1:5]
table.all[c(3,6,9,12,15), 4] <- bse.ab[1:5]

table.all[c(1,4,7,10,13), 5] <- coefs.ab.jbc[1:5]
table.all[c(3,6,9,12,15), 5] <- bse.ab.jbc[1:5]

table.all[c(1,4,7,10,13), 6] <- coefs.ab.jbc5[1:5]
table.all[c(3,6,9,12,15), 6] <- bse.ab.jbc5[1:5]


table.all[16, 1] <- lr.fe
table.all[17, 1] <- cse.lr.fe
table.all[18, 1] <- bse.lr.fe

table.all[16, 2] <- lr.jbc
#table.all[17, 2] <- cse.lr.fe
table.all[18, 2] <- bse.lr.jbc

table.all[16, 3] <- lr.abc4
table.all[18, 3] <- bse.lr.abc4

table.all[16, 4] <- lr.ab
table.all[18, 4] <- bse.lr.ab

table.all[16, 5] <- lr.ab.jbc
table.all[17, 5] <- cse.lr.ab
table.all[18, 5] <- bse.lr.ab.jbc

table.all[16, 6] <- lr.ab.jbc5
table.all[18, 6] <- bse.lr.ab.jbc5


table.all[1, ] <- 100 * table.all[1, ]
table.all[2, ] <- 100 * table.all[2, ]
table.all[3, ] <- 100 * table.all[3, ]
table.all[16, ] <- 100 * table.all[16, ]
table.all[17, ] <- 100 * table.all[17, ]
table.all[18, ] <- 100 * table.all[18, ]


xtable(table.all, digits=2)