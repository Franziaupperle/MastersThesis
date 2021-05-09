
## METHOD FOR PEFORMING DOUBLE-STABILITY

# no_runs_mc: number of MC iterations
# B: number of subsampling steps
# tresh: selected value for the treshold
db.stab_mcfun = function(no_runs_mc, B, tresh) {
  
  # generate relevant parameters
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  # start parallel computing
  para = parallel_db_stab(
    B = B,
    pf = pf,
    tresh = tresh,
    Rd2 = Rd2,
    Ry2 = Ry2
  )
  
  # function to execute simulation iterations parallel
  
  # B: number of subsampling steps
  # tresh: selected value for the treshold
  # Rd2: coefficient of determination in the treatment equation
  # Ry2: coefficient of determination in the main equation
  parallel_db_stab = function(B, pf, tresh, Rd2, Ry2) {
    alpha.est = NULL
    alpha.se = NULL
    ci = NULL
    laenge = NULL
    
    cat("start parallelisierung \n \n")
    
    # initialize cores
    numCores <- detectCores()
    cl <- makeCluster(numCores[1] - 1)
    registerDoParallel(cl)
    
    # introduce relevant functions to server
    HolpAlasso.set = HolpAlasso.set
    glmnet = glmnet
    cv.glmnet = cv.glmnet
    ci.coverfun = ci.coverfun
    IF.varestbiascorr = IF.varestbiascorr
    generate.data = generate.data
    db.stab.set = db.stab.set
    
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
        library(sandwich)
        library(lmtest)
        
        # initialize vecotrs
        M.yx.all = c()
        M.dx.all = c()
        
        set.seed(i + 12345)
        
        # get data
        data = generate.data(pf = pf,
                             Rd2 = Rd2,
                             Ry2 = Ry2)
        x = data$x
        y = data$y
        d = data$d
        
        z = cbind(d, x)
        
        # start repeated subsampling and model selection 
        for (b in 1:B) {
          # randomly selection subsample
          index.subsam = sample(seq(1, pf$n), 0.6 * pf$n)
          y.s = y[index.subsam]
          z.s = z[index.subsam, ]
          d.s = z.s[, 1]
          x.s = z.s[, -1]
          
          # perform double-selection with LASSO cross-validation
          M.yx = db.stab.set(x = x.s, y = y.s)
          M.dx = db.stab.set(x = x.s, y = d.s)
          
          # story selected coefficients
          M.yx.all = append(M.yx.all, M.yx)
          M.dx.all = append(M.dx.all, M.dx)
          
        }
        
        # get variables that have a relative inclusion frequency of at least 'tresh' in the LASSO selections of the main equation
        M.Freq = as.data.frame(table(M.yx.all) / B)
        M.hat = M.Freq[M.Freq$Freq > tresh, ]
        M.yx.hat = M.hat$M.yx.all
        
        # get variables that have a relative inclusion frequency of at least 'tresh' in the LASSO selections of the treatment equation
        M.Freq = as.data.frame(table(M.dx.all) / B)
        M.hat = M.Freq[M.Freq$Freq > tresh, ]
        M.dx.hat = M.hat$M.dx.all
        
        # unite the selected variables 
        M.hat = unique(sort(c(M.yx.hat, M.dx.hat)))
        # refit the model with the selected variables
        linearmodel = lm(y ~ d + x[, M.hat]) # Post-lasso step
        alpha.est = linearmodel$coef[2]
        # get robust standard errors
        linearmodel_HC3 = coeftest(linearmodel, vcov = vcovHC(linearmodel, type =
                                                                "HC3"))
        alpha.se = linearmodel_HC3[2, 2]
        
        # calculate relevant performance measures
        bias = alpha.est - pf$alpha0
        mse = (alpha.est - pf$alpha0) ^ 2 + alpha.se ^ 2
        confi_int = ci.coverfun(
          alpha.est = alpha.est,
          alpha.se = alpha.se,
          alpha0 = pf$alpha0
        )
        ci =  confi_int$truefalse
        laenge = confi_int$ci.upper - confi_int$ci.lower
        model.size = length(M.hat)
        
        # bind results of i-th MC-iteration in a vector
        result = cbind(alpha.est, alpha.se, bias, mse, ci, laenge, model.size)
        
        return(result)
        
        stopCluster(cl)
      }
    # return result of i-th MC iteration
    return(temp)
  }
  
  print(para)
  
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
