####################################################################;
# FUNCTIONS FOR PANEL PROBIT WITH INDIVIDUAL AND TIME EFFECTS ;
#####################################################################;


bias_parameter <- function(N, T, index, x, nomissing)
  {
  ws <- dnorm(index)^2 / (pnorm(index) * pnorm(-index))
  ws[is.na(ws)] <- 0
  FE <- kronecker(diag(1, N), matrix(1, nrow = T))
  TE <- kronecker(matrix(1, nrow = N), diag(1, T))
  FE <- FE[nomissing, ]
  TE <- TE[nomissing, ] 
  resx <- lm(x ~ cbind(FE[,-1], TE[,-1]), weights = ws)$res;
  B <- (1/2) * apply((t(index * ws * resx) %*% FE) / kronecker(ws %*% FE, matrix(1, nrow = ncol(x))), 1, mean)
  D <- (1/2) * apply((t(index * ws * resx) %*% TE) / kronecker(ws %*% TE, matrix(1, nrow = ncol(x))), 1, mean)
  W <- (t(ws * resx) %*% resx)/nrow(x)
  return(solve(W) %*% (B/T + D/N))
 }



bias_dynamic_parameter <- function(N, T, index, x, y, L2 = 0, L3 = 0, L4 = 0) {
 ws     <- dnorm(index)^2 / (pnorm(index) * pnorm(-index))
 FE     <- kronecker(diag(1, N), matrix(1, nrow = T))
 TE     <- kronecker(matrix(1, nrow = N), diag(1, T))
 resx   <- lm(x ~ cbind(FE[, -1], TE[, -1]), weights = ws)$res
 psi    <- matrix(ws * (y - pnorm(index)) / dnorm(index), nrow = T, ncol = N)
 lpsi   <- (T/(T - 1)) * as.vector(matrix(rbind(rep(0, N), psi[-T, ]), ncol = 1))
 l2psi  <- (T/(T - 2)) * as.vector(matrix(rbind(rep(0, N), rep(0, N), 
                                                psi[-c((T - 1):T), ]), ncol = 1))
 l3psi  <- (T/(T - 3)) * as.vector(matrix(rbind(rep(0, N), rep(0, N), 
                                                rep(0, N), psi[-c((T - 2):T), ]), ncol = 1))
 l4psi  <- (T/(T - 4)) * as.vector(matrix(rbind(rep(0, N), rep(0, N), rep(0, N), 
                                                rep(0, N), psi[-c((T - 3):T), ]), ncol = 1))
 B      <- (1/2) * apply((t(resx * ws * (index - 2 * lpsi - 2 * L2 * l2psi - 2 * L3 * l3psi - 2 * L4 * l4psi)) %*% FE) / kronecker(ws %*% FE, rep(1,ncol(x))), 1, mean)
 D      <- (1/2) * apply((t(resx * ws * index) %*% TE) / kronecker(ws %*% TE, rep(1,ncol(x))), 1, mean)
 W      <- (t(resx) %*% (ws * resx))/(N*T)
 return(solve(W) %*% (B/T + D/N))
}


# Version of speedglm function that returns linear predictors and design matrix;

spglm <- function(formula, data, family = gaussian(), weights = NULL, 
                   start = NULL, etastart = NULL, mustart = NULL, offset = NULL, 
                   maxit = 25, k = 2, sparse = NULL, set.default = list(), x = FALSE, ...) {
  call <- match.call()
  tf <- terms(formula)
  M <- model.frame(tf, data)
  y <- M[[1]]
  X <- model.matrix(tf, M)
  offset <- model.offset(M)
  intercept <- attributes(tf)$intercept
  set <- list(sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
              row.chunk = NULL, tol.solve = .Machine$double.eps, acc = 1e-08, 
              tol.values = 1e-07, tol.vectors = 1e-07, method = "eigen")
  nmsC <- names(set)
  set[(namc <- names(set.default))] <- set.default
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in set.default: ", 
            paste(noNms, collapse = ", "))
  rval <- speedglm.wfit(y = y, X = X, family = family, weights = weights, 
                        start = start, etastart = etastart, mustart = mustart, 
                        offset = offset, intercept = intercept, row.chunk = set$row.chunk, 
                        maxit = maxit, k = k, acc = set$acc, sparselim = set$sparselim, 
                        camp = set$camp, eigendec = set$eigendec, tol.solve = set$tol.solve, 
                        sparse = sparse, tol.values = set$tol.values, tol.vectors = set$tol.vectors, 
                        method = set$method)
  rval$terms <- tf
  rval$call <- call
  rval$linear.predictors <- as.vector(X %*% rval$coefficients) 
  if (x) 
    rval$x <- X
  class(rval) <- "speedglm"
  if ((rval$iter == maxit) & (!rval$convergence)) 
    warning("Maximum number of iterations reached without convergence")
  rval
}



