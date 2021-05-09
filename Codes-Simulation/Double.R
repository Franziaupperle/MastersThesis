
### Estimation with Double selection

# no_runs_mc: number of MC iterations
double.mcfun = function(no_runs_mc) {
  
  # generate relevant parameters
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  # start parallel computing
  para = double.parallel(
    no_runs_mc = no_runs_mc,
    pf = pf,
    Rd2 = Rd2,
    Ry2 = Ry2
  )
  
  # function to execute simulation iterations parallel
  
  # no_runs_mc: number of MC-iterations
  # pf: selected parameters
  # Rd2: coefficient of determination in the treatment equation
  # Ry2: coefficient of determination in the main equation
  double.parallel = function(no_runs_mc, pf, Rd2, Ry2) {
    
    # declare variables/vectors
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
    
    # introduce relevant functions to server
    ci.coverfun = ci.coverfun
    generate.data = generate.data
    
     # distribute MC-iterations
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
      
       # get data
      data = generate.data(pf = pf,
                           Rd2 = Rd2,
                           Ry2 = Ry2)
      x = data$x
      y = data$y
      d = data$d
      alpha0 = pf$alpha0
      
      # perform Double-Selection with function rlasso
      double.ATE = rlassoEffect(
        x = x,
        y = y,
        d = d,
        method = "double selection",
        I3 = NULL
      )
      
      
      # get results
      alpha.est = double.ATE$alpha
      alpha.se = double.ATE$se[1]
      model.size = length(double.ATE$coefficients.reg) - 1
      
      # calculate relevant performance measures
      bias = alpha.est - alpha0
      mse = (alpha.est - alpha0) ^ 2 + alpha.se ^ 2
      confi.int = ci.coverfun(alpha.est = alpha.est,
                              alpha.se = alpha.se,
                              alpha0 = alpha0)
      ci = confi.int$truefalse
      ci.laenge = confi.int$ci.upper - confi.int$ci.lower
      
      # bind results of i-th MC-iteration in a vector
      result = cbind(alpha.est, alpha.se, bias, mse, ci, ci.laenge, model.size)
      
      return(result)
      
      stopCluster(cl)
    }
    # return result of i-th MC iteration
    return(temp)
    
  }
  
  # return final results of Double-Stability with the selected parameter setting in the beginning
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




