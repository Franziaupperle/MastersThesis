
# # PODS - Split; Projection-On-Double-Selection Split

pods.split_mcfun = function(no_runs_mc, B, dfmin) {
  
  # generate relevant parameters
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  # start parallel computing
  para = parallel_pods_split(
    B = B,
    pf = pf,
    dfmin = dfmin,
    Rd2 = Rd2,
    Ry2 = Ry2
  )
  
  # function to execute simulation iterations parallel
  # B: number of subsampling steps
  # pf: generated parameters
  # Rd2: coefficient of determination in the treatment equation
  # Ry2: coefficient of determination in the main equation
  parallel_pods_split = function(B, pf, dfmin, Rd2, Ry2) {

    cat("start parallelisierung \n \n")

     # initialize cores
    numCores <- detectCores()
    cl <- makeCluster(numCores[1] - 1)
    registerDoParallel(cl)

     # introduce relevant functions to server
    HolpAlasso.set = HolpAlasso.set
    Projection.set = Projection.set
    glmnet = glmnet
    cv.glmnet = cv.glmnet
    ci.coverfun = ci.coverfun
    IF.varestbiascorr = IF.varestbiascorr
    generate.data = generate.data

     # distribute MC-iterations
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
      
        # declare variables
        alpha.est = NULL
        alpha.se = NULL
        bias = NULL
        mse = NULL
        ci = NULL
        laenge = NULL
        model.size = c()
        
        set.seed(i + 12345)
        
         # get data
        data = generate.data(pf = pf,
                             Rd2 = Rd2,
                             Ry2 = Ry2)
        x = data$x
        y = data$y
        d = data$d
        split.size = pf$split.size # number of observations used for model selection
        size.holp = pf$size.holp
        alpha0 = pf$alpha0
        
        # declare matrix used to store the selection indices of the sample
        Y.count = matrix(data = NA,
                         nrow = B,
                         ncol = pf$n)
        z = cbind(d, x)
        
         # start repeated subsampling, model selection treatment estimation
        for (b in 1:B) {
          # randomly selection subsample
          index.subsam = sample(seq(1, pf$n), split.size)
          y.s = y[index.subsam]
          z.s = z[index.subsam, ]
          d.s = z.s[, 1]
          x.s = z.s[, -1]
          index.test = setdiff(seq(1, pf$n), index.subsam)
          y.t = y[index.test]
          z.t = z[index.test, ]
          d.t = z.t[, 1]
          x.t = z.t[, -1]
          Y.count[b, ] = tabulate(index.subsam, pf$n)
          
          # perform model selection with PODS
          M.hat = Projection.set(
            y = y.s,
            d = d.s,
            x = x.s,
            size.holp = size.holp,
            dfmax = length(y.t)/2 - 3,
            dfmin = dfmin/2,
            HolpAlasso.set.fun = HolpAlasso.set
          )
          
          # refit the model
          alpha.est[b] = lm(y.t ~ d.t + x.t[, M.hat])$coef[2]
          model.size[b] = length(M.hat)
          
        }
        
        # get smoothed estimator and its standard eroor
        alpha.smoothed = mean(alpha.est)
        alpha.se = IF.varestbiascorr(Y.count, alpha.est, pf$n, split.size)
        
        # calculate relevant performance measures
        confi.int = ci.coverfun(alpha.est = alpha.smoothed,
                                alpha.se = alpha.se,
                                alpha0 = alpha0)
        
        bias = alpha.smoothed - alpha0
        mse = (alpha.smoothed - alpha0) ^ 2 + alpha.se ^ 2
        ci = confi.int$truefalse
        laenge = confi.int$ci.upper - confi.int$ci.lower
        modelsize = mean(model.size)
        min.modelsize = min(model.size)
        max.modelsize = max(model.size)
        
        # bind results of i-th MC-iteration in a vector
        result = cbind(alpha.smoothed, alpha.se, bias, mse, ci, laenge, modelsize, min.modelsize, max.modelsize)
        
        return(result)
        
        stopCluster(cl)
        
      }
    # return result of i-th MC iteration
    return(temp)
    
  }
  
  print(para)
  # return final results of PODS-Split with the selected parameter setting from the beginning
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
    dfmin = dfmin,
    tresh = tresh
  )
  
  return(result)
  
}
