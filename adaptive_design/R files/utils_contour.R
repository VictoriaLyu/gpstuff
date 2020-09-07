##' 1D test function
##' @param x test location(s)
fun1d <- function(x){
  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
  x^2 - 0.75^2
}

##' Modified Branin function
branin3 <- function(x){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  x1 <- x[,1] * 15
  x2 <- x[,2] * 15 - 5
  ((x1 - 5.1/(4 * pi^2) * (x2^2) + 5/pi * x2 - 20)^2 + (10 - 10/(8 * pi)) * cos(x2) - 181.47)/178
  
}

##' Modified Hartman6 function
hartman6m <- function(x){
  # works only on a vector
  hartman6_uni <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    x3 <- x[3]
    x4 <- x[4]
    x5 <- x[5]
    x6 <- x[6]
    a <- matrix(c(8, 3, 10, 3.50, 1.7, 6,
                  0.5, 8, 10, 1, 6, 9,
                  3, 3.5, 1.7, 8, 10, 6,
                  10, 6, 0.5, 8, 1, 9), 
                6, 4)
    p <- matrix(c(0.1312, 0.2329, 0.2348, 0.4047, 0.1696, 0.4135, 
                  0.1451, 0.8828, 0.5569, 0.8307, 0.3522, 0.8732, 0.0124, 
                  0.3736, 0.2883, 0.5743, 0.8283, 0.1004, 0.3047, 0.1091, 
                  0.5886, 0.9991, 0.665, 0.0381), 6, 4, byrow = TRUE)
    c <- c(0.2, 0.22, 0.28, 0.3)
    d <- matrix(0, 1, 4)
    for (i in seq(1, 4)) {
      d[i] <- sum(a[, i] * (x - p[, i])^2)
    }
    f <- -sum(c * exp(-d))
    return((f + 0.1)/0.1)
  }
  
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  return(apply(x, 1, hartman6_uni))
}

##' @noRd
##' @examples 
##' x <- seq(0,1, by = 0.001)
##' plot(x, am_put_oracle(x))
am_put_oracle <- function(x){
  # Initialize the parameters 
  bdr1 <- c(35.11630, 35.38308, 35.63883, 35.89936, 36.18202, 36.51018, 36.90874, 37.39437, 37.95088)
  K <- 40
  boundaries <- bdr1
  sigma <- 0.2
  r <- 0.06
  dt <- 0.04
  
  x <- 10 * x + 30
  len <- length(x)
  M <- length(boundaries)
  stopPayoff <- pmax(0, K-x)
  curX <- x*exp(rnorm(len)*sigma*sqrt(dt) + (r- sigma^2/2)*dt)
  payoff <- stopPayoff
  contNdx <- seq(1, len, by = 1)
  i <- 1
  while (i <=M && length(contNdx) > 0){
    payoff[contNdx]  <- exp(-i*dt*r) * pmax(0, K- curX[contNdx]);
    contNdx <- contNdx[( curX[contNdx] > boundaries[i] | payoff[contNdx] == 0)]
    if (any(contNdx != 0))
      curX[contNdx] <- curX[contNdx]*exp(rnorm(length(contNdx))*sigma * sqrt(dt) + (r- sigma^2/2)*dt)
    i <- i + 1
  }
  payoff[contNdx]  <- exp(-i*dt*r) * pmax(0, K- curX[contNdx]);
  y <- payoff-stopPayoff
  
}

##' @noRd
##' @examples 
##' ### Synthetic 1d case
##' set.seed(42)
##' print(genFun(0.75, fun1d, "tconstdf", "small"))
##' xgrid <- seq(0,1, length.out = 201)
##' Xgrid <- t(t(xgrid))
##' par(mfrow = c(3,2))
##' plot(xgrid, genFun(Xgrid, fun1d, "normal", "small"))
##' lines(xgrid, fun1d(Xgrid), col = 'red')
##' plot(xgrid, genFun(Xgrid, fun1d, "t_constdf", "middle"))
##' lines(xgrid, fun1d(Xgrid), col = 'red')
##' plot(xgrid, genFun(Xgrid, fun1d, "t_heterodf", "large"))
##' lines(xgrid, fun1d(Xgrid), col = 'red')
##' plot(xgrid, genFun(Xgrid, fun1d, "tconstdf", "hetero"))
##' lines(xgrid, fun1d(Xgrid), col = 'red')
##' plot(xgrid, genFun(Xgrid, fun1d, "mixed", "mixed"))
##' lines(xgrid, fun1d(Xgrid), col = 'red')
##' plot(xgrid, genFun(Xgrid, fun1d, "t_heterodf", "hetero"))
##' lines(xgrid, fun1d(Xgrid), col = 'red')
##' par(mfrow = c(1, 1))
##' @note r is supposed to be 1.
genFun <- function(X, fun, noiseStructure, noisevar){
  if(is.null(nrow(X))) X <- matrix(X, ncol = 1)
  f <- fun(X)
  
  sd <- switch(noisevar,
               "small" = 0.1,
               "middle" = 0.2,
               "large" = 0.5,
               "hetero" = 0.4 * (4 * X[,1] + 1),
               "mixed" = c(0.5, 1)
  )
  
  gmsamp <- function(f, sd) sample(rnorm(2, mean = f, sd=sd), 1)
  
  res <- switch(noiseStructure,
                "normal" = drop(rnorm(nrow(X), mean = f, sd = sd)),
                "t_constdf" = drop(f + sd * rt(n = nrow(X), df = 3)),
                "t_heterodf" = drop(f + sd * rt(n = nrow(X), df = 6 - 4 * X[,1])),
                "mixed" = drop(sapply(f, gmsamp, sd = sd)),
                "null" = f
  )
  return(res)
}

##' Performance evaluation for contour study
##' @noRd
##' @param model for prediction
##' @param xt design locations
##' @param ft corresponding true value
##' @param pcr critical probability pcr in eq. 5.2
##' @param lambda
gp_perf <- function(model, xt, ft, pcr = 0.4, lambda = 0.8){
  d <- ncol(xt)
  m <- nrow(xt)
  nt1 <- lambda*m
  preds <- predict(model, xt)
  Ef <- preds$mean
  Varf <- preds$sd2
  
  if(class(model) %in% c("homTP", "hetTP")){
    lee <- pt(-abs(Ef)/sqrt(Varf), df = model$nu + length(model$Z)) 
  }else{
    lee <- pnorm(-abs(Ef)/sqrt(Varf))
  } 
  
  if(d == 1){
    ee <- sum(lee)/m
    er <-  sum(Ef * ft < 0)/m
    bias <- (sum(Ef > 0 & ft < 0) - sum(Ef < 0 & ft > 0))/m
  }else{
    ee <- pcr * sum(lee[1:nt1])/nt1 + (1 - pcr) * sum(lee[-(1:nt1)])/(m - nt1)
    er <- pcr * sum(Ef[1:nt1] * ft[1:nt1] < 0)/nt1 + (1 - pcr) * sum(Ef[-(1:nt1)] * ft[-(1:nt1)] < 0)/(m - nt1)
    bias <- pcr * ((sum(Ef[1:nt1] > 0 & ft[1:nt1] < 0) - sum(Ef[1:nt1] < 0 & ft[1:nt1] > 0))/nt1) -
      (1 - pcr) * (sum(Ef[-(1:nt1)] > 0 & ft[-(1:nt1)] < 0) - sum(Ef[-(1:nt1)] < 0 & ft[-(1:nt1)] > 0))/(m - nt1)
  }
  
  return(list(ee = ee, er = er, bias = bias))
}


##' Multi-start optimization of contour criteria
crit_optim <- function(model, crit, 
                       control = list(method = "L-BFGS-B", tol_dist = 1e-4,
                                      multi.start = 20, maxit = 100), ...){
  if(is.null(control$method)) control$method <- "L-BFGS-B"
  if(is.null(control$tol_dist)) control$tol_dist <- 1e-4
  if(control$method == "L-BFGS-B"){
    # Xstart <- lhsDesign(control$multi.start, ncol(model$X0), seed = sample(1:2^15, 1))$design
    Xstart <- matrix(runif(control$multi.start * ncol(model$X0)), nrow = control$multi.start)
    
    res <- list(par = NA, value = -Inf, new = NA)
    
    for(i in 1:nrow(Xstart)){
      out <- optim(Xstart[i,, drop = FALSE], fn = match.fun(crit), method = "L-BFGS-B", lower = rep(0, d), upper = rep(1, d),
                   model = model, control = list(maxit = control$maxit, fnscale = -1, factr = 1e-3), ... = ...)
      if(out$value > res$value)
        res <- list(par = out$par, value = out$value, new = TRUE, id = NULL)
    }
  }else{
    if(control$method == "DE"){
    crit2 <- function(x, model, ...){
      return(-match.fun(crit)(x = x, model = model, ...))
    }
    res <- DEoptim(fn = crit2, lower = rep(0, d), upper = rep(1, d), model = model, ... = ...,
                   control = DEoptim.control(trace = FALSE, itermax = control$maxit, reltol = 1e-3))
    res <- list(par = matrix(res$optim$bestmem, nrow = 1), value = res$optim$bestval)
    }
    if(control$method == "pso"){
      res <- psoptim(par = rep(NA, ncol(model$X)), fn = match.fun(crit), lower = rep(0, d), upper = rep(1, d),
                     model = model, ... = ..., control = list(fnscale = -1, maxit = control$maxit, s = control$multi.start))
      res$par <- matrix(res$par, nrow = 1)
    }
  }
  
  ## Check if new design is not to close to existing design
  dists <- sqrt(hetGP:::distance_cpp(res$par, model$X0))
  if(min(dists) < control$tol_dist){
    res <- list(par = model$X0[which.min(dists),,drop = FALSE],
                value = match.fun(crit)(x = model$X0[which.min(dists),, drop = F], model = model, ... = ...),
                new = FALSE, id = which.min(dists))
  }
  
  ## Ensure that not more than 10 x dimension of the problem is allocated to a single design, otherwise explore
  max_mult <- 10*ncol(model$X0)
  if(max(model$mult) > max_mult && all(res$par == model$X0[which.max(model$mult),])){
    newX <- matrix(runif(ncol(model$X0)), nrow = 1)
    res <- list(par = newX,
                value = match.fun(crit)(x = newX, model = model, ... = ...),
                new = TRUE)
  }
  
  return(res)
}

