
# # PODS - Projection onto Double Selection

pods.mcfun = function(no_runs_mc) {
  
  # generate relevant parameters
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  # start parallel computing
  para = pods.parallel(
    no_runs_mc = no_runs_mc,
    pf = pf,
    Rd2 = Rd2,
    Ry2 = Ry2
  )
  
  # function to execute simulation iterations parallel
  
  pods.parallel = function(no_runs_mc, pf, Rd2, Ry2) {
    
    # initialize variables
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
    IF.varestbiascorr = IF.varestbiascorr
    generate.data = generate.data
    Projection = Projection
    ci.coverfun = ci.coverfun
    
     # distribute MC-iterations
    temp <- foreach(
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
      
      proj = Projection(
        alpha0 = alpha0,
        y = y,
        d = d,
        x = x
      )
      
      alpha.est = proj$alpha
      alpha.se = proj$se
      
      bias = alpha.est - alpha0
      mse = (alpha.est - alpha0) ^ 2 + alpha.se ^ 2
      confi = ci.coverfun(alpha.est, alpha.se, alpha0)
      ci = confi$truefalse
      ci.laenge = (confi$ci.upper - confi$ci.lower)
      model.size = proj$model.size
      
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
