# 7/29/2020 Shuowen Chen
# This script produces tables for abc with different trimming parameters and bse
library(lfe)
library(boot)
library(Hmisc)
# load data (sdf)
load("covidsince0307.RData")
# add weight to sdf
sdf$sweight <- 1
# create week factor to the data
sdf$week <- as.factor(strftime(sdf$date, format = "%V"))
sdf$grouppol <- (sdf$pshelter + sdf$pmovie + sdf$prestaurant)/3
sdf$grouppol2 <- (sdf$pshelter + sdf$pmovie + sdf$prestaurant + sdf$pnonessential)/4

createfmla_fe <- function(yvar, xvars, iv = "0", cluster = "state",
                          fe = c("state + month", "state*month")) {
  fe <- match.arg(fe)
  rhs <- paste(xvars, collapse = " + ")
  rhs_felm <- paste(rhs," | ", fe, " | ", iv, " | ", 
                    paste(cluster, collapse = " + ")) 
  rhs_lm <- paste(rhs, "+", fe, collapse = " ")
  # formula for felm
  fm1 <- as.formula(paste(yvar, "~", rhs_felm, sep = " "))
  # formula for lm (for the use of analytical bias correction)
  fm2 <- as.formula(paste(yvar, "~", rhs_lm, sep = " "))
  return(list(fm1, fm2)) 
}

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

reg <- function(df, # data
                yvar, # dependent variable (1 of 4 bevahiors if pib)
                pols, # policies
                bvars, # behavior (NULL if not behavior)
                infovars,
                x, # controls 
                iv, # iv for felm, if not needed just input "0" 
                l = 14,
                fixedeffect) {
  if (l == 0) {
    # This is for pib
    p <- pols
    b <- bvars
  } else {
    # This is for pbiy and piy
    p <- sprintf("lag(%s, %d)", pols, l)
    b <- sprintf("lag(%s, %d)", bvars, l)
    info <- sprintf("lag(%s, %d)", infovars, l)
    # this is for lm command because it cannot run if data 
    # columns are not actually created
    lagterm <- sapply(c(pols, bvars, infovars), 
                      function(x) out <- lag(sdf[, x], l))
    colnames(lagterm) <- c(p, b, info)
    df <- cbind(df, lagterm)
  }
  # create regression formula
  xvars <- c(p, b, x)
  fmla <- createfmla_fe(yvar, xvars, fe = fixedeffect, iv = iv)
  environment(fmla[[2]]) <- environment()
  # Uncorrected Estimator
  whole <- felm(fmla[[1]], data = df, weights = df$sweight)  
  coef_nobc <- coef(whole)
  whole_lm <- lm(fmla[[2]], data = df, x = TRUE, weights = df$sweight,
                 na.action = na.omit)
  
  abc7 <- abc(df, whole_lm, trim = 7, n = length(levels(df$state)),
              rhs = xvars, unc = coef_nobc)
  abc14 <- abc(df, whole_lm, trim = 14, n = length(levels(df$state)),
               rhs = xvars, unc = coef_nobc)
  abc21 <- abc(df, whole_lm, trim = 21, n = length(levels(df$state)),
               rhs = xvars, unc = coef_nobc)
  abc28 <- abc(df, whole_lm, trim = 28, n = length(levels(df$state)),
               rhs = xvars, unc = coef_nobc)
  abc35 <- abc(df, whole_lm, trim = 35, n = length(levels(df$state)),
               rhs = xvars, unc = coef_nobc)
  abc42 <- abc(df, whole_lm, trim = 42, n = length(levels(df$state)),
               rhs = xvars, unc = coef_nobc)
  abc49 <- abc(df, whole_lm, trim = 49, n = length(levels(df$state)),
               rhs = xvars, unc = coef_nobc)
  
  # Sum of policy coefficients (uncorrected and corrected)
  peff <- c(sum(coef_nobc[p]), sum(abc7[p]),
            sum(abc14[p]), sum(abc21[p]),
            sum(abc28[p]), sum(abc35[p]),
            sum(abc42[p]), sum(abc49[p]))
  
  # Weighted sum of behavioral coefficients (uncorrected and corrected)
  if (!is.null(bvars)) {
    # weight is from April 1st to 10th
    w <- colMeans(subset(df, df$date >= as.Date("2020-04-01") &
                           df$date <= as.Date("2020-04-10"))[, bvars])
    beff <- c(sum(coef_nobc[b]*w), sum(abc7[b]*w),
              sum(abc14[b]*w), sum(abc21[b]*w),
              sum(abc28[b]*w), sum(abc35[b]*w),
              sum(abc42[b]*w), sum(abc49[b]*w))
  } else {
    # for piy and pib no such metric
    beff <- NULL
  }
  # outputs: noncorrected, corrected, sum of policy and behavioral coefficients
  return(list(nobc = coef_nobc, abc7 = abc7, abc14 = abc14,
              abc21 = abc21, abc28 = abc28, abc35 = abc35,
              abc42 = abc42, abc49 = abc49, sumpolicy = peff, 
              sumbehavior = beff))
}

# The following function extends the function mainregressions
# NOTE: there is one change in the code
# the original code put xlist in the controls, but xlist is specified as an empty list
# in the program, so I'm not sure what it stands for. I think it's meant for creating a list
# so as to loop over choices of specifications. I deleted it in the argument. I checked the 
# results with and without xlist in the control for fe = "0" in the original code 
# and they are identical. 

mainregressions <- function(df, # data
                            yvar, # dependent variable (Y)
                            pols, # policy variables
                            bvars, # behavior variables (B)
                            infovars, # information structures
                            tvars, # key confounder from SIR
                            ivlist = "0",
                            L = 14, 
                            fixed = "state + month") {
  # This is considering both pols and interacted of pmask with months
  # plist <- list(pols, c("pmask.april","pmask.may", pols[-1]))
  # for pols only
  plist <- list(pols)
  
  ijs <- expand.grid(1:length(plist), 1:length(infovars))
  
  # The function loops over plist, infovars to run different specifications
  pbiy <- apply(ijs, 1, function(ij) {
    reg(df, yvar, plist[[ij[1]]], bvars, infovars[[ij[2]]],
        c(sprintf("lag(%s, %d)", infovars[[ij[2]]], L), tvars), 
        ivlist, l = L, fixedeffect = fixed)
  })
  
  piy <- apply(ijs, 1, function(ij) {
    reg(df, yvar, plist[[ij[1]]], NULL, infovars[[ij[2]]],
        c(sprintf("lag(%s, %d)", infovars[[ij[2]]], L), tvars), 
        ivlist, l = L, fixedeffect = fixed)
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
      reg(df, bvars[ij[1]], plist[[ij[2]]], NULL, 
          infovars[[k]], c(infovars[[k]]), ivlist, 
          l = 0, fixedeffect = fixed)
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
pols <- c("pmaskbus", "pk12", "grouppol", "pnonessential") 
pols2 <- c("pmaskbus", "pk12", "grouppol2") 
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

# statstics in each bootstrap stat
bstat <- function(data, mle) {
  tot <- mainregressions(data, yvar2, pols, bvars, infovarsd, NULL, L = 21)
  return(unlist(tot))
}

bstat2 <- function(data, mle) {
  tot <- mainregressions(data, yvar2, pols2, bvars, infovarsd, NULL, L = 21)
  return(unlist(tot))
}


ncores <- 1
set.seed(88, kind = "L'Ecuyer-CMRG")
#parallel::mc.reset.stream()
curious7 <- boot(data = sdf, statistic = bstat, 
                sim = "parametric", ran.gen = data_wb, mle = 0, 
                parallel = "multicore", ncpus = ncores, R = 500)
set.seed(88, kind = "L'Ecuyer-CMRG")
curious2 <- boot(data = sdf, statistic = bstat2, 
                 sim = "parametric", ran.gen = data_wb, mle = 0, 
                 parallel = "multicore", ncpus = ncores, R = 500)

# Produce Tables
bse <- bsemat(curious7)
bse2 <- bsemat(curious2)

pib_g1 <- hybridtab(bse, "pib", pols, bvars, NULL, infovarsd, 21, 8)
pib_g2 <- hybridtab(bse2, "pib", pols2, bvars, NULL, infovarsd, 21, 8)

pbiy_g1 <- hybridtab(bse, "pbiy", pols, bvars, NULL, infovarsd, 21, 8)
pbiy_g2 <- hybridtab(bse2, "pbiy", pols2, bvars, NULL, infovarsd, 21, 8)

piy_g1 <- hybridtab(bse, "piy", pols, bvars, NULL, infovarsd, 21, 8)
piy_g2 <- hybridtab(bse2, "piy", pols2, bvars, NULL, infovarsd, 21, 8)



