################################################################################
### Equivalent of Matlab tests in R for TPs
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  print("No argument provided, use d = 1 by default")
  d <- 1 # 1, 2, 6
}else{
  d <- as.numeric(args[1]) # dimension: 1 or 2 or 6
}


source("../examples/contour/utils_contour.R")
library(R.matlab)
library(DiceDesign)
library(hetGP)
library(parallel)


set.seed(42)

if(d == 1){
  k <- 100
  I <- 10
  budget <- 100
  m0 <- 1000
  xt <- xtt <- seq(0, 1, length.out = m0)
  fun <- fun1d
  ft <- fun(xtt)
  xt <- xtt <- matrix(xtt, ncol = 1)
  w <- rep(1, nrow(xtt))
  plotF <- FALSE
  control <- list(multi.start = 20, maxit = 100, tol_dist = 1e-4)
  lower <- rep(0.18, d)
  # Note: upper is increased as iteration increases (from 0.5 to 5)
  uppers <- matrix(seq(0.5, 2, length.out = budget - I), nrow = budget - I,ncol = d) 
  known <- list(nu = 3)
  noiseControl <- list(sigma2_bounds = c(0.3^2, 2^2), g_bounds = c(0.3^2, 2^2))
}

if(d == 2){
  k <- 20
  I <- 20
  budget <- 150
  m0 <- 50
  x1 <- x2 <- seq(0, 1, length.out = m0)
  xt <- as.matrix(expand.grid(x1, x2))
  xtt <- readMat("../examples/contour/testPoints2D.mat")$xtt
  fun <- branin3
  ft <- fun(xtt)
  w <- rep(c(0.4, 0.6), times = c(0.8, 0.2)*nrow(xtt))
  w <- w * nrow(xtt) / sum(w)
  plotF <- FALSE
  control <- list(multi.start = 20, maxit = 100, tol_dist = 1e-4)
  lower <- rep(0.18, d)
  # Note: upper is increased as iteration increases (from 0.5 to 5)
  uppers <- matrix(seq(0.5, 8, length.out = budget - I), nrow = budget - I,ncol = d) 
  known <- list(nu = 3)
  noiseControl <- list(sigma2_bounds = c(0.3^2, 2^2), g_bounds = c(0.3^2, 2^2))
}

if(d == 6){
  k <- 20
  I <- 60
  budget <- 1000
  xt <- xtt <- readMat("../examples/contour/testPoints6D.mat")$xtt
  fun <- hartman6m
  ft <- apply(xtt, 1, fun)
  w <- rep(c(0.4, 0.6), times = c(0.8, 0.2)*nrow(xtt))
  w <- w * nrow(xtt) / sum(w)
  # control <- list(multi.start = 20, maxit = 20)
  # control <- list(multi.start = 20, maxit = 25, method = "pso")
  control <- list(multi.start = 20, maxit = 10, method = "DE", tol_dist = 1e-4)
  lower <- rep(0.18, d)
  # Note: upper is increased as iteration increases (from 0.5 to 5)
  uppers <- matrix(seq(0.5, 8, length.out = budget - I), nrow = budget - I,ncol = d)
  known <- list(nu = 3)
  noiseControl <- list(sigma2_bounds = c(0.3^2, 2^2), g_bounds = c(0.3^2, 2^2))
  library(DEoptim)
  library(pso)
}


# Test cases
cases <- list(c('t_constdf', 'small'),
              c('t_constdf', 'large'),
              c('mixed', 'mixed'),
              c('t_heterodf', 'hetero'))
crits <- c("crit_MCU", 
           "crit_cSUR",
           "crit_ICU",
           "crit_tMSE")

casescrits <- expand.grid(cases, crits)
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

run_cont <- function(ii){
  
  for(i in 1:k)
  {
    print(i)
    t0 <- Sys.time()
    case <- unlist(casescrits[ii, 1])
    crit <- as.character(casescrits[ii, 2])
    X <- maximinESE_LHS(lhsDesign(n = I, dimension = d, seed = i)$design)$design
    Y <- genFun(X, fun = fun, noiseStructure = case[1], noisevar = case[2])
    
    if(d == 1 && plotF) plot(X, Y, xlim = c(0, 1), ylim = c(-2, 2))
    if(d == 2 && plotF) plot(X, xlim = c(0,1), ylim = c(0, 1))
    
    model <- mleFit(X = X, Z = Y, lower = lower, upper = uppers[1,], known = known, noiseControl = noiseControl)
    for(j in 1:(budget - I)){
      if(j %% 50 == 0) print(j)
      if(crit == "crit_ICU"){
        # t0 <- Sys.time()
        preds <- predict(model, x = xtt)
        kxprime <- cov_gen(X1 = model$X0, X2 = xtt, theta = model$theta, type = model$covtype)
        sol <- crit_optim(model, crit = crit, control = control, Xref = xtt, w = w,
                          preds = preds, kxprime = kxprime)
        # print(Sys.time() - t0)
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
      X <- rbind(X, Xnew); Ynew <- genFun(Xnew, fun = fun, noiseStructure = case[1], case[2])
      Y <- c(Y, Ynew)
      if(i %in% up_it) maxit_up <- 1 else maxit_up <- 0
      model <- update(model, Xnew=Xnew, Znew=Ynew, ginit = model$g, maxit = maxit_up)
      
      if(d == 1 && plotF) points(Xnew, Ynew)
      if(d == 2 && plotF) points(Xnew)
      
      ## periodically restart model fit
      if(j %in% up_it){ 
        mod2 <- mleFit(X=list(X0=model$X0, Z0=model$Z0, mult=model$mult), Z=model$Z,
                       lower = lower, upper = uppers[j,], known = known, noiseControl = noiseControl)
        # if(mod2$ll > model$ll) 
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
  
  writeMat(paste0(class(model), "_", crit, "_", case[1], "_", case[2], "_", d, "d_rev1.mat"),
           type = paste0(class(model), "_", crit),
           noise_case = paste0(case[1], "_", case[2]),
           x_seq = x_seq, y_seq = y_seq, ee = ee, er = er,
           bias = bias, Ef = Ef, Varf = Varf, timings = timings,
           sigma2s = sigma2s, thetas = thetas, nus = nus, noisevars = noisevars)
  i
  
}
# }

# run_cont(1)
# run_cont(4)
# mclapply(9:12, run_cont, mc.cores = min(16, detectCores() - 1))
mclapply(1:nrow(casescrits), run_cont, mc.cores = min(16, detectCores()))

### Some graphs
if(FALSE){
  data_res <- readMat("../examples/contour/New_Results_1d/homTP_crit_MCU_t_heterodf_hetero_1d_rev1.mat")
  matplot(data_res$ee, type = 'l')
  lines(rowMeans(data_res$ee), lwd = 2)
  matplot(data_res$er, type = 'l')
  matplot(data_res$bias, type = 'l')
  plot(data_res$x.seq, data_res$y.seq)
  plot(data_res$Ef[,90,2])
  boxplot(data_res$ee[90,])
  
  dat1 <- readMat("../examples/contour/Results_2d/homTP_crit_MCU_t_constdf_large_2d.mat")
  dat2 <- readMat("../examples/contour/Results_2d/homTP_crit_MCU_t_constdf_large_2d_rev1.mat")
  par(mfrow = c(3,2))
  plot(rowMeans(dat1$ee))
  lines(rowMeans(dat2$ee))
  boxplot(dat1$ee[130,], dat2$ee[130,])
  plot(rowMeans(dat1$er))
  lines(rowMeans(dat2$er))
  boxplot(dat1$er[130,], dat2$er[130,])
  plot(rowMeans(dat1$bias))
  lines(rowMeans(dat2$bias))
  abline(h = 0)
  boxplot(dat1$bias[130,], dat2$bias[130,])
  par(mfrow = c(1,1))
}

### For Table 3
if(FALSE){
  crit <- "ICU"
  dim <- 6
  ddir <- paste0("../examples/contour/New_Results4_", dim, "d/")
  res_ts <- readMat(paste0(ddir, "homTP_crit_", crit, "_t_constdf_small_", dim, "d_rev1.mat"))
  res_tl <- readMat(paste0(ddir, "homTP_crit_", crit, "_t_constdf_large_", dim, "d_rev1.mat"))
  res_mm <- readMat(paste0(ddir, "homTP_crit_", crit, "_mixed_mixed_", dim, "d_rev1.mat"))
  res_hh <- readMat(paste0(ddir, "homTP_crit_", crit, "_t_heterodf_hetero_", dim, "d_rev1.mat"))
  
  cat(signif(mean(res_ts$er[nrow(res_ts$er),]), 3) * 100, signif(sd(res_ts$er[nrow(res_ts$er),]), 3) * 100, "\n")
  cat(signif(mean(res_tl$er[nrow(res_tl$er),]), 3) * 100, signif(sd(res_tl$er[nrow(res_tl$er),]), 3) * 100, "\n")
  cat(signif(mean(res_mm$er[nrow(res_mm$er),]), 3) * 100, signif(sd(res_mm$er[nrow(res_mm$er),]), 3) * 100, "\n")
  cat(signif(mean(res_hh$er[nrow(res_hh$er),]), 3) * 100, signif(sd(res_hh$er[nrow(res_hh$er),]), 3) * 100, "\n")
}

## Comparing results
if(FALSE){
  crit <- "ICU"
  dim <- 6
  ddir <- paste0("../examples/contour/Results_", dim, "d/")
  res_ts <- readMat(paste0(ddir, "homTP_crit_", crit, "_t_constdf_small_", dim, "d_rev1.mat"))
  res_tl <- readMat(paste0(ddir, "homTP_crit_", crit, "_t_constdf_large_", dim, "d_rev1.mat"))
  res_mm <- readMat(paste0(ddir, "homTP_crit_", crit, "_mixed_mixed_", dim, "d_rev1.mat"))
  res_hh <- readMat(paste0(ddir, "homTP_crit_", crit, "_t_heterodf_hetero_", dim, "d_rev1.mat"))
  
  ddir1 <- paste0("../examples/contour/New_Results_", dim, "d/")
  res_ts1 <- readMat(paste0(ddir1, "homTP_crit_", crit, "_t_constdf_small_", dim, "d_rev1.mat"))
  res_tl1 <- readMat(paste0(ddir1, "homTP_crit_", crit, "_t_constdf_large_", dim, "d_rev1.mat"))
  res_mm1 <- readMat(paste0(ddir1, "homTP_crit_", crit, "_mixed_mixed_", dim, "d_rev1.mat"))
  res_hh1 <- readMat(paste0(ddir1, "homTP_crit_", crit, "_t_heterodf_hetero_", dim, "d_rev1.mat"))
  
  ddir2 <- paste0("../examples/contour/New_Results4_", dim, "d/")
  res_ts2 <- readMat(paste0(ddir2, "homTP_crit_", crit, "_t_constdf_small_", dim, "d_rev1.mat"))
  res_tl2 <- readMat(paste0(ddir2, "homTP_crit_", crit, "_t_constdf_large_", dim, "d_rev1.mat"))
  res_mm2 <- readMat(paste0(ddir2, "homTP_crit_", crit, "_mixed_mixed_", dim, "d_rev1.mat"))
  res_hh2 <- readMat(paste0(ddir2, "homTP_crit_", crit, "_t_heterodf_hetero_", dim, "d_rev1.mat"))
  
  # ddir3 <- paste0("../examples/contour/New_Results3_", dim, "d/")
  # res_ts3 <- readMat(paste0(ddir3, "homTP_crit_", crit, "_t_constdf_small_", dim, "d_rev1.mat"))
  # res_tl3 <- readMat(paste0(ddir3, "homTP_crit_", crit, "_t_constdf_large_", dim, "d_rev1.mat"))
  # res_mm3 <- readMat(paste0(ddir3, "homTP_crit_", crit, "_mixed_mixed_", dim, "d_rev1.mat"))
  # res_hh3 <- readMat(paste0(ddir3, "homTP_crit_", crit, "_t_heterodf_hetero_", dim, "d_rev1.mat"))
  
  
  par(mfrow = c(2, 1))
  lr <- nrow(res_ts$er) #
  boxplot(res_ts$ee[lr,], res_ts1$ee[lr,], res_ts2$ee[lr,], #res_ts3$ee[lr,],
          res_tl$ee[lr,], res_tl1$ee[lr,], res_tl2$ee[lr,], #res_tl3$ee[lr,],
          res_mm$ee[lr,], res_mm1$ee[lr,], res_mm2$ee[lr,], #res_mm3$ee[lr,],
          res_hh$ee[lr,], res_hh1$ee[lr,], res_hh2$ee[lr,]) #res_hh3$ee[lr,])
  title(paste(crit, " d=", dim))
  boxplot(res_ts$er[lr,], res_ts1$er[lr,], res_ts2$er[lr,], #res_ts3$er[lr,],
          res_tl$er[lr,], res_tl1$er[lr,], res_tl2$er[lr,], #res_tl3$er[lr,],
          res_mm$er[lr,], res_mm1$er[lr,], res_mm2$er[lr,], #res_mm3$er[lr,],
          res_hh$er[lr,], res_hh1$er[lr,], res_hh2$er[lr,]) #res_hh3$er[lr,])
  par(mfrow = c(1, 1))
  
}

if(FALSE){
  boxplot(res_ts$ee[90,], res_tl$ee[90,], res_mm$ee[90,], res_hh$ee[90,])
  boxplot(res_ts$er[90,], res_tl$er[90,], res_mm$er[90,], res_hh$er[90,])
}


# ################################################################################
# if(FALSE){
#   X <- matrix(seq(0,1, length.out = 20))
#   Y <- genFun(X, fun, "small", "tconstdf")
#   
#   plot(X, Y, ylim = c(-1, 1))
#   model <- mleHetGP(X = X, Z = Y, lower = 0.001, upper = 10)
#   
#   
#   for(i in 1:180){
#     ## find the next point and update 
#     # crit_values <- crit_tMSE(x = Xgrid, model = model, seps = 0.01)
#     crit_values <- crit_cSUR(x = Xgrid, model = model)
#     
#     Xnew <- Xgrid[which.max(crit_values),,drop = FALSE]
#     X <- c(X, Xnew); Ynew <- genFun(Xnew, fun, "small", "tconstdf"); Y <- c(Y, Ynew)
#     model <- update(model, Xnew=Xnew, Znew=Ynew, ginit=model$g*1.01)
#     
#     points(Xnew, Ynew)
#     ## periodically restart model fit
#     if(i %% 25 == 0){ 
#       mod2 <- mleHetTP(X=list(X0=model$X0, Z0=model$Z0, mult=model$mult), Z=model$Z, init = list(nu = 6),
#                        lower=0.0001, upper=1)
#       if(mod2$ll > model$ll) model <- mod2
#     } 
#   }
#   
#   preds <- predict(model, Xgrid)
#   ## Display mean predictive surface
#   lines(xgrid, preds$mean, col = 'red', lwd = 2)
#   ## Display 95% confidence intervals
#   lines(xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#   lines(xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#   ## Display 95% prediction intervals
#   lines(xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2 + preds$nugs)), 
#         col = 3, lty = 2)
#   lines(xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2 + preds$nugs)), 
#         col = 3, lty = 2)
#   lines(xgrid, fun(Xgrid), col = 'blue')
#   lines(xgrid, fun(Xgrid) + 0.1*qt(p = 0.95, df = 3), col = 'blue')
#   
#   ids <- which(pnorm(-abs(preds$mean)/sqrt(preds$sd2)) > 0.05)
#   points(Xgrid[ids], rep(-1,length(ids)), pch = 20)
#   
#   plot(X)
#   abline(h = 0.75, col = 'red')
# }


