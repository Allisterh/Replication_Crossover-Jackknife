# 7/19/2020 Shuowen Chen
# This scripts implements analytical bias correction for two way fixed effects
library(lfe)
pathSC <- "/Users/shuowenchen/desktop/research/Panel-Cross-over-jackknife/application"
setwd(pathSC)
source('Code/auxiliaryfunctions.R')
# load data (sdf)
load("covidsince0307.RData")
# add weight to sdf
sdf$sweight <- 1
# create week factor to the data
sdf$week <- as.factor(strftime(sdf$date, format = "%V"))

reg <- function(df, # data
                yvar, # dependent variable (1 of 4 bevahiors if pib)
                pols, # policies
                bvars, # behavior (NULL if not behavior)
                infovars,
                x, # controls 
                iv, # iv for felm, if not needed just input "0" 
                l = 14,
                trim, # trimming parameter for abc
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
                      function(x) out <- lag(df[, x], l))
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
  
  coef_abc <- abc(df, whole_lm, trim, n = length(levels(df$state)),
                  rhs = xvars, unc = coef_nobc)
  
  # Sum of policy coefficients (uncorrected and corrected)
  peff <- c(sum(coef_nobc[p]), sum(coef_abc[p]))
  
  # Weighted sum of behavioral coefficients (uncorrected and corrected)
  if (!is.null(bvars)) {
    # weight is from April 1st to 10th
    w <- colMeans(subset(df, df$date >= as.Date("2020-04-01") &
                           df$date <= as.Date("2020-04-10"))[, bvars])
    beff <- c(sum(coef_nobc[b]*w), sum(coef_abc[b]*w))
  } else {
    # for piy and pib no such metric
    beff <- NULL
  }
  # outputs: noncorrected, corrected, sum of policy and behavioral coefficients
  return(list(nobc = coef_nobc, abc = coef_abc, 
              sumpolicy = peff, sumbehavior = beff))
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
                            trim, # trimming parameter for abc
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
           ivlist, l = L, trim, fixedeffect = fixed)
  })
  
  piy <- apply(ijs, 1, function(ij) {
    reg(df, yvar, plist[[ij[1]]], NULL, infovars[[ij[2]]],
           c(sprintf("lag(%s, %d)", infovars[[ij[2]]], L), tvars), 
           ivlist, l = L, trim, fixedeffect = fixed)
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
      reg(df, bvars[ij[1]], plist[[ij[2]]], NULL, infovars[[k]],
             c(infovars[[k]]), ivlist, l = 0, trim,
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
# Note: to run only one structure, say I, just specify 
# infovars <- list(c("dlogdc", "logdc"))
infovars <- list(c("dlogdc", "logdc"),
                 c("dlogdc", "logdc", "dlogdc.national", "logdc.national")) # cases

infovarsd <- list(c("dlogdd", "logdd"),
                  c("dlogdd", "logdd", "dlogdd.national", "logdd.national")) # deaths


lala <- lapply(1:7, function(x) {
  result <- mainregressions(sdf, yvar2, pols, bvars, infovarsd, NULL, 
                            L = 21, trim = 7*x)
  pib <- result$pib
  pbiy <- result$pbiy
  piy <- result$piy
  ijs <- expand.grid(1:length(infovarsd), 1:length(bvars))
  pibmat <- apply(ijs, 1, function(ij){
    out <- rbind(matrix(head(unlist(pib[[ij[[1]]]][[ij[[2]]]]), -2), ncol = 2),
                 matrix(tail(unlist(pib[[ij[[1]]]][[ij[[2]]]]), 2), ncol = 2))
  })
  pbiy1 <- rbind(matrix(head(unlist(pbiy[[1]]), -4), ncol = 2), 
                 unlist(pbiy[[1]][[3]]), unlist(pbiy[[1]][[4]]))
  
  pbiy3 <- rbind(matrix(head(unlist(pbiy[[2]]), -4), ncol = 2), 
                 unlist(pbiy[[2]][[3]]), unlist(pbiy[[2]][[4]]))
  
  piymat <- sapply(c(1, 2), function(x){
    out <- rbind(matrix(head(unlist(piy[[x]]), -2), ncol = 2), 
                 tail(unlist(piy[[x]]), 2))
  })
  return(list(pib = pibmat, pbiy1 = pbiy1, pbiy3 = pbiy3, piy = piymat))
})


piy_abc <- function(lala) {
  piy_1 <- cbind(lala[[1]]$piy[[1]], lala[[2]]$piy[[1]][, 2])
  piy_3 <- cbind(lala[[1]]$piy[[2]], lala[[2]]$piy[[2]][, 2])
  for (i in 3:7) {
    piy_1 <- cbind(piy_1, lala[[i]]$piy[[1]][, 2])
    piy_3 <- cbind(piy_3, lala[[i]]$piy[[2]][, 2])
  }
  rownames(piy_1) <- c(sprintf("lag(%s, %d)", c(pols, infovarsd[[1]]), 21), 
                       "sumpol")
  rownames(piy_3) <- c(sprintf("lag(%s, %d)", c(pols, infovarsd[[2]]), 21),
                       "sumpol")
  colnames(piy_1) <- c("nobc", "abc7", "abc14", "abc21", "abc28",
                       "abc35", "abc42", "abc49")
  colnames(piy_3) <- colnames(piy_1)
  return(list(piy1 = piy_1, piy3 = piy_3))
} 

pbiy_abc <- function(lala) {
  pbiy_1 <- cbind(lala[[1]]$pbiy1, lala[[2]]$pbiy1[, 2])
  pbiy_3 <- cbind(lala[[1]]$pbiy3, lala[[2]]$pbiy3[, 2])
  for (i in 3:7) {
    pbiy_1 <- cbind(pbiy_1, lala[[i]]$pbiy1[, 2])
    pbiy_3 <- cbind(pbiy_3, lala[[i]]$pbiy3[, 2])
  }
  rownames(pbiy_1) <- c(sprintf("lag(%s, %d)", c(pols, bvars, infovarsd[[1]]), 21), 
                        "sumpol", "sumb")
  rownames(pbiy_3) <- c(sprintf("lag(%s, %d)", c(pols, bvars, infovarsd[[2]]), 21),
                        "sumpol", "sumb")
  colnames(pbiy_1) <- c("nobc", "abc7", "abc14", "abc21", "abc28",
                        "abc35", "abc42", "abc49")
  colnames(pbiy_3) <- colnames(pbiy_1)
  return(list(pbiy1 = pbiy_1, pbiy3 = pbiy_3))
}

pib_final <- sapply(1:4, function(x){
  out1 <- cbind(lala[[1]]$pib[[(2*x - 1)]], lala[[2]]$pib[[(2*x - 1)]][, 2])
  out3 <- cbind(lala[[1]]$pib[[(2*x)]], lala[[2]]$pib[[(2*x)]][, 2])
  for (i in 3:7) {
    out1 <- cbind(out1, lala[[i]]$pib[[(2*x - 1)]][, 2])
    out3 <- cbind(out3, lala[[i]]$pib[[(2*x)]][, 2])
  }
  colnames(out1) <- c("nobc", "abc7", "abc14", "abc21", "abc28",
                      "abc35", "abc42", "abc49")
  colnames(out3) <- colnames(out1)
  rownames(out1) <- c(pols, infovarsd[[1]], "sumpol")
  rownames(out3) <- c(pols, infovarsd[[2]], "sumpol")
  return(list(out1, out3))
})
# naming for ease of tabulation
names(pib_final) <- paste(c("info1", "info3"), rep(bvars, each = 2), sep = "")
piy_final <- piy_abc(lala)
pbiy_final <- pbiy_abc(lala)

# use xtable to tabulate results in latex









