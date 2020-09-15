################################################################################
### Equivalent of Matlab tests in R (for TPs) - Controlled simulations case 
################################################################################

## Note: this script assumes that the working directory is the one with this file

source("./utils_contour.R")
library(R.matlab)
library(DiceDesign)
library(hetGP)
library(parallel)

set.seed(42)

d <- 2

k <- 20
I <- 20
budget <- 150
m0 <- 50
x1 <- x2 <- seq(0, 1, length.out = m0)
xt <- as.matrix(expand.grid(x1, x2))
xtt <- xt
w <- rep(c(0.4, 0.6), times = c(0.8, 0.2)*nrow(xtt))
w <- w * nrow(xtt) / sum(w)
plotF <- FALSE
control <- list(multi.start = 20, maxit = 100, tol_dist = 1e-4)
lower <- rep(0.18, d)
# Note: upper is increased as iteration increases (from 0.5 to 5)
uppers <- matrix(seq(0.5, 8, length.out = budget - I), nrow = budget - I,ncol = d) 
known <- list(nu = 3)
noiseControl <- list(sigma2_bounds = c(0.3^2, 2^2), g_bounds = c(0.3^2, 2^2))


# Test cases

ftype <- c("gp", "tp")

cases <- list(c('normal', 'small'),
              c('normal', 'large'),
              c('t_constdf', 'small'),
              c('t_constdf', 'large'))

crits <- c("crit_MCU", 
           "crit_cSUR",
           "crit_ICU",
           "crit_tMSE")

casescrits <- expand.grid(ftype, cases, crits)
mleFit <- match.fun("mleHomTP")

# initialize variables
m <- nrow(xt)
x_seq <- array(dim = c(budget, d, k))
y_seq <- matrix(NA, budget, k)
er <- ee <- bias <- metric <- matrix(NA, budget - I, k)
Ef <- Varf <- array(dim = c(nrow(xt), budget - I, k)) # just for plotting?
sigma2s <- noisevars <- nus <- matrix(NA, budget - I, k) # store process variance
thetas <- array(NA, dim = c(budget - I, d, k)) # lengthscales
timings <- numeric(k)

up_it <- unique(c(1, 2, 4, 8, 16, 32, 40, 64, 80, 90, 128, 130, 190,
                  250, 350, 450, budget - I, seq(1, budget-I, by = 100)))

n0 <- 40 # numver of design used to defined the GP/TP test functions

#' @param X matrix of designs (one per row)
#' @param ftype either \code{gp} or \code{tp}
#' @param case vector with noise type (normal or student) and noise variance (small or large)
#' @param fixedData list with elements X0, Y0, Kmatinv
#' @param irun which element of the data to use for consistency with Matlab
fun_wrap <- function(X, case, fixedData, irun){
  if(is.null(nrow(X))){
    stop("Maybe a bug here!") # could be nrow = 1 instead
    X <- matrix(X, ncol = 1)
  }
  
  noiseStructure <- case[1]  
  noisevar <- case[2]
  
  sd0 <- switch(noisevar,
                "small" = 0.1,
                "large" = 1
  )
  
  f <- cov_gen(X, fixedData$X0[,,irun], theta = c(2*0.2^2, 2*0.2^2), type = "Gaussian") %*% (fixedData$Kmatinv[,,irun] %*% fixedData$f[,irun])
  
  res <- switch(noiseStructure,
                "normal" = drop(rnorm(nrow(X), mean = f, sd = sd0)),
                "t_constdf" = drop(f + sd0 * rt(n = nrow(X), df = 3)),
                "null" = f
  )
  return(res)
}


run_cont <- function(ii){
  
  for(i in 1:k)
  {
    print(i)
    t0 <- Sys.time()
    ftype <- casescrits[ii, 1]
    case <- unlist(casescrits[ii, 2])
    crit <- as.character(casescrits[ii, 3])
    
    if(ftype == "gp"){
      fixedData <- readMat("./gp.mat")
    }
    if(ftype == "tp"){
      fixedData <- readMat("./tp.mat")
    }
    
    ft <- cov_gen(xtt, fixedData$X0[,,i], theta = c(2*0.2^2, 2*0.2^2), type = "Gaussian") %*% (fixedData$Kmatinv[,,i] %*% fixedData$f[,i])
    
    
    X <- maximinESE_LHS(lhsDesign(n = I, dimension = d, seed = i)$design)$design
    Y <- fun_wrap(X = X, case = case, fixedData = fixedData, irun = i)
    
    if(d == 1 && plotF) plot(X, Y, xlim = c(0, 1), ylim = c(-2, 2))
    if(d == 2 && plotF) plot(X, xlim = c(0,1), ylim = c(0, 1))
    
    model <- mleFit(X = X, Z = Y, lower = lower, upper = uppers[1,], known = known, noiseControl = noiseControl)
    for(j in 1:(budget - I)){
      if(j %% 50 == 0) print(j)
      if(crit == "crit_ICU"){
        preds <- predict(model, x = xtt)
        kxprime <- cov_gen(X1 = model$X0, X2 = xtt, theta = model$theta, type = model$covtype)
        sol <- crit_optim(model, crit = crit, control = control, Xref = xtt, w = w,
                          preds = preds, kxprime = kxprime)
      }
      if(crit == "crit_tMSE"){
        sol <- crit_optim(model, crit = crit, control = control, seps = 0)
      }
      
      if(crit == "crit_cSUR"){
        sol <- crit_optim(model, crit = crit, control = control)
      }
      
      if(crit == "crit_MCU"){
        preds <- predict(model, x = xtt)
        gamma <- diff(quantile(preds$mean, probs = c(0.25,0.75))) / (3 * mean(sqrt(preds$sd2)))
        sol <- crit_optim(model, crit = crit, control = control, gamma = gamma)
      }
      
      metric[j, i] <- sol$value
      
      Xnew <- sol$par
      X <- rbind(X, Xnew); Ynew <- fun_wrap(X = Xnew, case = case, fixedData = fixedData, irun = i)
      Y <- c(Y, Ynew)
      if(j %in% up_it) maxit_up <- 1 else maxit_up <- 0
      model <- update(model, Xnew=Xnew, Znew=Ynew, ginit = model$g, maxit = maxit_up)
      
      if(d == 1 && plotF) points(Xnew, Ynew)
      if(d == 2 && plotF) points(Xnew)
      
      ## periodically restart model fit
      if(j %in% up_it){ 
        mod2 <- mleFit(X=list(X0=model$X0, Z0=model$Z0, mult=model$mult), Z=model$Z,
                       lower = lower, upper = uppers[j,], known = known, noiseControl = noiseControl)
        model <- mod2
      } 
      
      # hyperparameter values
      
      thetas[j, ,i] <- model$theta
      noisevars[j, i] <- model$g
      if(class(model) %in% c("homTP", "hetTP")){
        nus[j, i] <- model$nu
        sigma2s[j, i] <- model$sigma2
      } 
      
      
      # Performance evaluation
      perfs <- gp_perf(model, xtt, ft)
      er[j, i] <- perfs$er
      ee[j, i] <- perfs$ee
      bias[j, i] <- perfs$bias
      preds <- predict(model, xt)
      Ef[, j, i] <- preds$mean
      Varf[, j, i] <- preds$sd2
      
    }
    
    tf <- Sys.time()
    
    if(d == 1 && plotF){
      preds <- predict(model, matrix(xt, ncol = 1))
      ## Display mean predictive surface
      lines(xt, preds$mean, col = 'red', lwd = 2)
      ## Display 95% confidence intervals
      lines(xt, qnorm(0.05, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
      lines(xt, qnorm(0.95, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
      ## Display 95% prediction intervals
      lines(xt, qnorm(0.05, preds$mean, sqrt(preds$sd2 + preds$nugs)), 
            col = 3, lty = 2)
      lines(xt, qnorm(0.95, preds$mean, sqrt(preds$sd2 + preds$nugs)), 
            col = 3, lty = 2)
      lines(xt, fun(xt), col = 'blue')
      ids <- which(pnorm(-abs(preds$mean)/sqrt(preds$sd2)) > 0.05)
      points(xt[ids], rep(min(fun(xt)),length(ids)), pch = 20)
    }
    
    # Save run data
    x_seq[,,i] <- X
    y_seq[,i] <- Y
    timings[i] <- difftime(tf, t0, units = "secs")
  }
  
  writeMat(paste0(class(model), "_", crit, "_", ftype, "_", case[1], "_", case[2], "_", d, "d_simu.mat"),
           type = paste0(class(model), "_", crit),
           noise_case = paste0(case[1], "_", case[2]),
           x_seq = x_seq, y_seq = y_seq, ee = ee, er = er,
           bias = bias, Ef = Ef, Varf = Varf, timings = timings,
           sigma2s = sigma2s, thetas = thetas, nus = nus, noisevars = noisevars)
  i
  
}
# }

# run_cont(1)
mclapply(1:nrow(casescrits), run_cont, mc.cores = min(16, detectCores()))


## To exploit results
if(FALSE){
  crit <- "ICU"
  dim <- 2
  ddir <- paste0("./")
  ftype <- "gp"
  res_ts <- readMat(paste0(ddir, "homTP_crit_", crit, "_", ftype, "_t_constdf_small_", dim, "d_simu.mat"))
  res_tl <- readMat(paste0(ddir, "homTP_crit_", crit, "_", ftype, "_t_constdf_large_", dim, "d_simu.mat"))
  res_ns <- readMat(paste0(ddir, "homTP_crit_", crit, "_", ftype, "_normal_small_", dim, "d_simu.mat"))
  res_nl <- readMat(paste0(ddir, "homTP_crit_", crit, "_", ftype, "_normal_large_", dim, "d_simu.mat"))
  
  cat(signif(mean(res_ns$er[nrow(res_ns$er),]), 3) * 100, signif(sd(res_ns$er[nrow(res_ns$er),]), 3) * 100, "\n")
  cat(signif(mean(res_nl$er[nrow(res_nl$er),]), 3) * 100, signif(sd(res_nl$er[nrow(res_nl$er),]), 3) * 100, "\n")
  cat(signif(mean(res_ts$er[nrow(res_ts$er),]), 3) * 100, signif(sd(res_ts$er[nrow(res_ts$er),]), 3) * 100, "\n")
  cat(signif(mean(res_tl$er[nrow(res_tl$er),]), 3) * 100, signif(sd(res_tl$er[nrow(res_tl$er),]), 3) * 100, "\n")
}
