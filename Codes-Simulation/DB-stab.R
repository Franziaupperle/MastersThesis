
## Stability Double Selection

db.stab_mcfun = function(no_runs_mc, B, tresh) {
  
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  para = parallel_db_stab(
    B = B,
    pf = pf,
    tresh = tresh,
    Rd2 = Rd2,
    Ry2 = Ry2
  )
  
  parallel_db_stab = function(B, pf, tresh, Rd2, Ry2) {
    alpha.est = NULL
    alpha.se = NULL
    ci = NULL
    laenge = NULL
    
    cat("start parallelisierung \n \n")
    
    numCores <- detectCores()
    
    cl <- makeCluster(numCores[1] - 1)
    registerDoParallel(cl)
    
    HolpAlasso.set = HolpAlasso.set
    glmnet = glmnet
    cv.glmnet = cv.glmnet
    ci.coverfun = ci.coverfun
    IF.varestbiascorr = IF.varestbiascorr
    generate.data = generate.data
    db.stab.set = db.stab.set
    
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
        
        M.yx.all = c()
        M.dx.all = c()
        
        set.seed(i + 12345)
        
        data = generate.data(pf = pf,
                             Rd2 = Rd2,
                             Ry2 = Ry2)
        x = data$x
        y = data$y
        d = data$d
        
        z = cbind(d, x)
        
        for (b in 1:B) {
          index.subsam = sample(seq(1, pf$n), 0.6 * pf$n)
          y.s = y[index.subsam]
          z.s = z[index.subsam, ]
          d.s = z.s[, 1]
          x.s = z.s[, -1]
          
          M.yx = db.stab.set(x = x.s, y = y.s)
          M.dx = db.stab.set(x = x.s, y = d.s)
          
          # DOUBLE SELECTION
          # M.yx = HolpAlasso.set(
          #   y = y.s,
          #   x = x.s,
          #   m0 = 1,
          #   size.holp = pf$p,
          #   default1 = FALSE,
          #   default2 = TRUE,
          #   dfmax = pf$p,
          #   dfmin = 0
          # )
          
          
          # M.dx = HolpAlasso.set(
          #   y = d.s,
          #   x = x.s,
          #   m0 = 1,
          #   size.holp = pf$p, 
          #   default1 = FALSE,
          #   default2 = TRUE,
          #   dfmax = pf$p,
          #   dfmin = 0
          # )
          
          
          M.yx.all = append(M.yx.all, M.yx)
          M.dx.all = append(M.dx.all, M.dx)
          
        }
        
        M.Freq = as.data.frame(table(M.yx.all) / B)
        M.hat = M.Freq[M.Freq$Freq > tresh, ]
        M.yx.hat = M.hat$M.yx.all
        
        M.Freq = as.data.frame(table(M.dx.all) / B)
        M.hat = M.Freq[M.Freq$Freq > tresh, ]
        M.dx.hat = M.hat$M.dx.all
        
        M.hat = unique(sort(c(M.yx.hat, M.dx.hat))) # M.hat does not include the treatment
        
        linearmodel = lm(y ~ d + x[, M.hat]) # Post-lasso step
        alpha.est = linearmodel$coef[2]
        linearmodel_HC3 = coeftest(linearmodel, vcov = vcovHC(linearmodel, type =
                                                                "HC3"))
        alpha.se = linearmodel_HC3[2, 2]
        
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
        
        result = cbind(alpha.est, alpha.se, bias, mse, ci, laenge, model.size)
        
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
    dfmin = dfmin,
    tresh = tresh
  )
  
  return(result)
  
}
