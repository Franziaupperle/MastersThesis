
### Estimation with Double selection

double.mcfun = function(no_runs_mc) {
  
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  para = double.parallel(
    no_runs_mc = no_runs_mc,
    pf = pf,
    Rd2 = Rd2,
    Ry2 = Ry2
  )
  
  double.parallel = function(no_runs_mc, pf, Rd2, Ry2) {
    
    alpha.est = NULL
    alpha.se = NULL
    bias = c()
    mse = c()
    ci = c()
    ci.laenge = c()
    model.size = c()
    
    cat("start parallelisierung \n \n")
    
     # initialize cores
    numCores <- detectCores()
    cl <- makeCluster(numCores[1] - 1)
    registerDoParallel(cl)
    
    ci.coverfun = ci.coverfun
    generate.data = generate.data
    
    temp <- foreach (
      i = 1:no_runs_mc,
      .combine = 'rbind',
      .inorder = TRUE,
      .multicombine = TRUE,
      .maxcombine = 1001
    ) %dopar% {
      
      library(MASS)
      library(hdm)
      
      set.seed(i + 12345)
      
      data = generate.data(pf = pf,
                           Rd2 = Rd2,
                           Ry2 = Ry2)
      x = data$x
      y = data$y
      d = data$d
      alpha0 = pf$alpha0
      
      double.ATE = rlassoEffect(
        x = x,
        y = y,
        d = d,
        method = "double selection",
        I3 = NULL
      )
      
      alpha.est = double.ATE$alpha
      alpha.se = double.ATE$se[1]
      model.size = length(double.ATE$coefficients.reg) - 1
      
      bias = alpha.est - alpha0
      mse = (alpha.est - alpha0) ^ 2 + alpha.se ^ 2
      confi.int = ci.coverfun(alpha.est = alpha.est,
                              alpha.se = alpha.se,
                              alpha0 = alpha0)
      ci = confi.int$truefalse
      ci.laenge = confi.int$ci.upper - confi.int$ci.lower
      
      result = cbind(alpha.est, alpha.se, bias, mse, ci, ci.laenge, model.size)
      
      return(result)
      
      stopCluster(cl)
    }
    
    return(temp)
    
  }
  
  result = list(
    Result = para,
    Bias = round(sqrt(n) * mean(para[, 3]), 2),
    SE.Bias = round(sqrt(n) * std.error(para[, 3]), 2),
    MSE = round(n * mean(para[, 4]), 2),
    SE.MSE = round(n * std.error(para[, 4]), 2),
    Cover = round(mean(para[, 5]), 2),
    SE.Cover = round(std.error(para[, 5]), 2),
    Length = round(mean(para[, 6]), 2),
    SE.Length = round(std.error(para[, 6]), 2),
    model.size = para[, 7],
    beta = beta,
    gamma = gamma,
    type = type,
    Ry2 = Ry2,
    Rd2 = Rd2
  )
  
  return(result)
  
}




