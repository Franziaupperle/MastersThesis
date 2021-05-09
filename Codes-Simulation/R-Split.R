
### R - Split

r.split_mcfun = function(no_runs_mc, B, dfmin) {
  
  # generate relevant parameters
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  # start parallel computing
  para = parallel_r_split(
    B = B,
    pf = pf,
    dfmin = dfmin,
    Rd2 = Rd2,
    Ry2 = Ry2
  )
  
  # function to execute simulation iterations parallel
  
  parallel_r_split = function(B, pf, dfmin, Rd2, Ry2) {
    
    cat("start parallelisierung \n \n")
    
     # initialize cores
    numCores <- detectCores()
    cl <- makeCluster(numCores[1] - 1)
    registerDoParallel(cl)
    
     # introduce relevant functions to server
    generate.data = generate.data
    HolpAlasso.set = HolpAlasso.set
    glmnet = glmnet
    cv.glmnet = cv.glmnet
    ci.coverfun = ci.coverfun
    IF.varestbiascorr = IF.varestbiascorr
    Projection.set = Projection.set
    
    temp <-
      foreach (
        i = 1:no_runs_mc,
        .combine = 'rbind',
        .inorder = TRUE,
        .multicombine = TRUE,
        .maxcombine = 1001
      ) %dopar% {
        
        library(glmnet)
        library(MASS)
        
        alpha.est = c()
        model.size = c()
        
        set.seed(i + 12345)
        
        data = generate.data(pf = pf,
                             Rd2 = Rd2,
                             Ry2 = Ry2)
        x = data$x
        y = data$y
        d = data$d
        alpha0 = pf$alpha0
        
        Y.count = matrix(data = NA,
                         nrow = B,
                         ncol = pf$n)
        z = cbind(d, x)
       
        
        for (b in 1:B) {
          index.subsam = sample(seq(1, pf$n), pf$split.size)
          y.s = y[index.subsam]
          z.s = z[index.subsam,]
          index.test = setdiff(seq(1, pf$n), index.subsam)
          y.t = y[index.test]
          z.t = z[index.test, ]
          Y.count[b, ] = tabulate(index.subsam, pf$n)
          
          M.hat = HolpAlasso.set(
            y = y.s,
            x = z.s,
            m0 = 1,
            size.holp = pf$size.holp,
            default = TRUE,
            dfmax = length(y.t) - 6,
            dfmin = dfmin
          )
          
          alpha.est[b] = lm(y.t ~ z.t[, M.hat])$coef[2]
          model.size[b] = length(M.hat) - 1
          
        }
        
        alpha.smoothed = mean(alpha.est)
        alpha.se = IF.varestbiascorr(Y.count, alpha.est, pf$n, pf$split.size)
        confi.int = ci.coverfun(
          alpha.est = alpha.smoothed,
          alpha.se = alpha.se,
          alpha0 = pf$alpha0
        )
        
        bias = alpha.smoothed - pf$alpha0
        mse = (alpha.smoothed - pf$alpha0) ^ 2 + alpha.se ^ 2
        ci =  confi.int$truefalse
        laenge = confi.int$ci.upper - confi.int$ci.lower
        modelsize = mean(model.size)
        min.modelsize = min(model.size)
        max.modelsize = max(model.size)
        
        result = cbind(alpha.smoothed,
                       alpha.se,
                       bias,
                       mse,
                       ci,
                       laenge,
                       modelsize,
                       min.modelsize,
                       max.modelsize)
        
        return(result)
        
        stopCluster(cl)
        
      }
    
    return(temp)
  }
  
  print(para)
  
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
    beta = beta,
    gamma = gamma,
    type = type,
    Ry2 = Ry2,
    Rd2 = Rd2,
    dfmin = dfmin
  )
  
  return(result)
  
}
