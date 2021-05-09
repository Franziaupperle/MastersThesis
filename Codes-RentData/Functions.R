
# R-SPLIT

Split.smooth <- function(y, d, x, B) {
  
  z = cbind(d,x)
  n = length(y)
  p = length(x[1,])
  
  alpha.est = NULL
  alpha.se = NULL
  M.all = c()
  model.size = c()
  
  Y.count = matrix(data = NA, nrow = B, ncol = n)
  C.count = matrix(data = NA, nrow = B, ncol = p)
  
  for (b in 1:B) {
    
    set.seed(b+12345)
    index.subsam = sample(seq(1,n), 0.7*n) 
    y.s = y[index.subsam]
    z.s = z[index.subsam,]
    x.s = z.s[,-1]
    index.test = setdiff( seq(1,n), index.subsam )
    y.t = y[index.test]
    z.t = z[index.test,]
    x.t = z.t[,-1]
    d.t = z.t[,1]
    Y.count[b,] = tabulate(index.subsam,n)
    
    M.hat = select.model(y = y.s, x = x.s)
    
    alpha.est[b] = lm(y.t ~ d.t + x.t[,M.hat])$coef[2]
    
    model.size[b] = length(M.hat)
    
    cat("Iteration number: ")
    print(b)
    
    M.all = append(M.all, M.hat)
    C.count[b,] = tabulate(M.hat, p)
    
  }
  
    C.rel = colMeans(C.count)

    M.hat09 = which(C.rel > 0.9)
    M.hat08 = which(C.rel > 0.8)
    M.hat07 = which(C.rel > 0.7)
    M.hat06 = which(C.rel > 0.6)
  
    M.Freq = as.data.frame(table(M.all)/B)

  alpha.smoothed = mean(alpha.est)
  alpha.se = IF.varestbiascorr(Y.count, alpha.est, n, 0.7*n)
  
  result = list(alpha.smoothed = alpha.smoothed, 
                alpha.est = alpha.est,
                alpha.se = alpha.se, 
                Covariates06 = colnames(x[,c(M.hat06)]),
                Covariates07 = colnames(x[,c(M.hat07)]),
                Covariates08 = colnames(x[,c(M.hat08)]),
                Covariates09 = colnames(x[,c(M.hat09)]),
                Freq = table(M.all),
                Mean.model.size = mean(model.size),
                model.size = model.size,
                C.count = C.count)
  
  return(result)
  
}


# PODS-SPLIT

PODS.Split <- function(y, d, x, B) {
  
  n = length(y)
  p = length(x[1,])
  
  Y.count = matrix(data = NA, nrow = B, ncol = n)
  C.count = matrix(data = NA, nrow = B, ncol = p)
  z = cbind(d,x)
  model.size = c()
  alpha.est = c()
  M.all = c()
  
  for(b in 1:B) {
    
    set.seed(b+12345)
    index.subsam = sample(seq(1,n), 0.7*n)
    y.s = y[index.subsam]
    z.s = z[index.subsam,]
    d.s = z.s[,1]
    x.s = z.s[,-1]
    index.test = setdiff( seq(1, n), index.subsam )
    y.t = y[index.test]
    z.t = z[index.test,]
    d.t = z.t[,1]
    x.t = z.t[,-1]
    Y.count[b,] = tabulate(index.subsam,n)
    
    M.hat = Projection.set(y = y.s, d = d.s, x = x.s)
    
    alpha.est[b] = lm(y.t ~ d.t + x.t[,M.hat])$coefficients[2]
    model.size[b] = length(M.hat)
    
    cat("Iteration no: ")
    print(b)
    
    M.all = append(M.all, M.hat)
    C.count[b,] = tabulate(M.hat, p)
    
  }

  C.rel = colMeans(C.count)
  
  M.hat09 = which(C.rel > 0.9)
  M.hat08 = which(C.rel > 0.8)
  M.hat07 = which(C.rel > 0.7)
  M.hat06 = which(C.rel > 0.6)
  
  alpha.smoothed = mean(alpha.est)
  alpha.se = IF.varestbiascorr(Y.count, alpha.est, n, 0.7*n)
  
  result = list(alpha.smoothed = alpha.smoothed, 
                alpha.est = alpha.est,
                alpha.se = alpha.se, 
                Covariates06 = colnames(x)[c(M.hat06)],
                Covariates07 = colnames(x)[c(M.hat07)],
                Covariates08 = colnames(x)[c(M.hat08)],
                Covariates09 = colnames(x)[c(M.hat09)],
                Mean.model.size = mean(model.size),
                model.size = model.size,
                C.count = C.count)
  
  return(result)
  
}

# DOUBLE STABILITY

Double.Stability <- function(y, d, x, B) {
  
  alpha.est = NULL
  alpha.se = NULL
  
  library(MASS)
  library(sandwich) 
  library(lmtest)
  
  n = length(y)
  p = length(x[1,])       
  M.yx.all = c()
  M.dx.all = c()
  z = cbind(d,x)
  Cy.count = matrix(data = NA, nrow = B, ncol = p)
  Cd.count = matrix(data = NA, nrow = B, ncol = p)
  
  for(b in 1:B) {
    
    set.seed(b+12345)
    index.subsam = sample(seq(1,n), 0.7*n)
    y.s = y[index.subsam]
    z.s = z[index.subsam,]
    d.s = z.s[,1] 
    x.s = z.s[,-1]
    
    # DOUBLE SELECTION
    M.dx = select.model(y = d.s, x = x.s)
    M.yx = select.model(y = y.s, x = x.s)
    
    M.yx.all = append(M.yx.all, M.yx)
    M.dx.all = append(M.dx.all, M.dx)
    cat("Iteration no: ")
    print(b)
    Cy.count[b,] = tabulate(M.yx, p)
    Cd.count[b,] = tabulate(M.dx, p)
    
  }
  
  Cy.rel = colMeans(Cy.count)
  Cd.rel = colMeans(Cd.count)

  Cy095 = which(Cy.rel > 0.95)
  Cd095 = which(Cd.rel > 0.95)
  M.hat095 = unique(sort(c(Cy095,Cd095)))
  Cy09 = which(Cy.rel > 0.9)
  Cd09 = which(Cd.rel > 0.9)
  M.hat09 = unique(sort(c(Cy09,Cd09)))
  Cy08 = which(Cy.rel > 0.8)
  Cd08 = which(Cd.rel > 0.8)
  M.hat08 = unique(sort(c(Cy08,Cd08)))
  Cy07 = which(Cy.rel > 0.7)
  Cd07 = which(Cd.rel > 0.7)
  M.hat07 = unique(sort(c(Cy07,Cd07)))
  Cy06 = which(Cy.rel > 0.6)
  Cd06 = which(Cd.rel > 0.6)
  M.hat06 = unique(sort(c(Cy06,Cd06)))
  
  
  lm06 = lm(y ~ d + x[,M.hat06]) # Post-lasso OLS step
  alpha.est06 = lm06$coef[2]
  lm_HC106 = coeftest(lm06, vcov = vcovHC(lm06, type="HC1"))
  alpha.se06 = lm_HC106[2,2]


  lm07 = lm(y ~ d + x[,M.hat07]) # Post-lasso OLS step
  alpha.est07 = lm07$coef[2]
  lm_HC107 = coeftest(lm07, vcov = vcovHC(lm07, type="HC1"))
  alpha.se07 = lm_HC107[2,2]


  lm08 = lm(y ~ d + x[,M.hat08]) # Post-lasso OLS step
  alpha.est08 = lm08$coef[2]
  lm_HC108 = coeftest(lm08, vcov = vcovHC(lm08, type="HC1"))
  alpha.se08 = lm_HC108[2,2]

  
  lm09 = lm(y ~ d + x[,M.hat09]) # Post-lasso OLS step
  alpha.est09 = lm09$coef[2]
  lm_HC109 = coeftest(lm09, vcov = vcovHC(lm09, type="HC1"))
  alpha.se09 = lm_HC109[2,2]
  
  lm095 = lm(y ~ d + x[,M.hat095]) # Post-lasso OLS step
  alpha.est095 = lm095$coef[2]
  lm_HC1095 = coeftest(lm095, vcov = vcovHC(lm095, type="HC1"))
  alpha.se095 = lm_HC1095[2,2]
  
  result = list(lm06 = summary(lm06),
                lmhc1_06 = lm_HC106,
                alpha.est06 = alpha.est06,
                alpha.se06 = alpha.se06,
                M.hat06 = colnames(x)[c(M.hat06)],
                lm07 = summary(lm07),
                lmhc1_07 = lm_HC107,
                alpha.est07 = alpha.est07,
                alpha.se07 = alpha.se07,
                M.hat07 = colnames(x)[c(M.hat07)],
                lm08 = summary(lm08),
                lmhc1_08 = lm_HC108,
                alpha.est08 = alpha.est08,
                alpha.se08 = alpha.se08,
                M.hat08 = colnames(x)[c(M.hat08)],
                lm09 = summary(lm09),
                lmhc1_09 = lm_HC109,
                alpha.est09 = alpha.est09,
                alpha.se09 = alpha.se09,
                M.hat09 = colnames(x)[c(M.hat09)],
                lm095 = summary(lm095),
                lmhc1_095 = lm_HC1095,
                alpha.est095 = alpha.est095,
                alpha.se095 = alpha.se095,
                M.hat095 = colnames(x)[c(M.hat095)],
                Cd.count = Cd.count,
                Cy.count = Cy.count,
                M.hat095 = M.hat095,
                M.hat09 = M.hat09,
                M.hat08 = M.hat08,
                M.hat07 = M.hat07,
                M.hat06 = M.hat06)
  
  return(result)
  
}

# Function for model selection in R-Split, PODS-Split and Double-Stability

select.model = function(y, x){
  
  p = length(x[1,])
  M = seq(1,p,1)
  
  # Perform model selection with either Rigorous Lasso oder Cross-Validation

  # index.lasso = rlasso(y ~ x)$index
  # 
  # M.hat = sort(M[index.lasso], decreasing = FALSE)
  # 
  # if(sum(index.lasso) == 0){
  #   M.hat = NULL
  # }


  fit.lasso = cv.glmnet(x = x, y = y,
                        standardize = TRUE, intercept = TRUE)
  lambda = fit.lasso$lambda.min
  fit.lasso0 = glmnet(x = x, y = y, standardize = TRUE, intercept = TRUE,
                      lambda = lambda)
  theta.lasso = as.vector(coef(fit.lasso0))
  index.lasso = theta.lasso[-1]!=0

  M.hat = sort(M[index.lasso], decreasing = FALSE)

  if(sum(index.lasso) == 0){
    M.hat = NULL
  }
  print(lambda)

  return(M.hat) 
}

# Function to calculate standard error in R-Split and PODS-Split 

IF.varestbiascorr = function(Y.count, alpha.est, n, split.size) {
  
  Brep = length(alpha.est)
  n2 = n-split.size
  
  cov.vec = t(sweep(Y.count, 2, apply(Y.count, 2, mean))) %*% (alpha.est - mean(alpha.est))/Brep
  var.est = (sum(cov.vec^2))
  r.corr = (n/split.size)^2 * ((n-1)/n)^2
  var.est0 = var.est * r.corr
  
  var.bias = sum( (alpha.est - mean(alpha.est))^2 )/Brep^2 * n  * n2/(n-n2)
  
  return( sqrt(abs(var.est0 - var.bias))) 
  
}

# Function to select a model and perform OLS in PODS

Projection = function(y, d, x) { 
  
  p = length(x[1,])
  n = length(y)
  
  M = seq(1,p)
  fit.dx = rlasso(d ~ x)$index # Step 1: select a set of variables M_hat_D which are associated with D
  
  if(sum(fit.dx) == 0){ # index of selected variables (logical vector)
    fit.dx = 1
  }
  
  M.dx = M[fit.dx] # FRISCH WAUGH THEOREM - variables, that are associated with D are kicked out
  
  x.res = lm(x[,-M.dx]~ d + x[,M.dx])$res # Step 2: to remove the components associated with D, we project (X, Y) onto a space which is orthogonal
  y.res = lm(y ~ d + x[,M.dx])$res # to the space spanned by D and X_M_hat_D
  
  fit.yx = rlasso(y.res ~ x.res)$index  # Step 3: selecting additional variables that are expected to have low correlation with D and so the
  
  temp0 = M[-M.dx]; # over-fitting bias can be controlled; select a model fit.yx for regressing y.res on x.res
  if(sum(fit.yx) == 0){ # select a model M_hat_Y_star for regressing Y on D and X_M_hat_star
    fit.yx = 1
  }
  
  M.yx = temp0[fit.yx] 
  
  M.hat = sort(c(M.dx, M.yx)) # unite relevant covariates 
  
  reg.dx = lm(d ~ x[,M.hat])
  reg.yx = lm(y ~ d + x[,M.hat]) # Step 4: Regress Y on D and M.hat to get alpha.hat, which is the estimated coefficient of D
  
  nu = reg.dx$residuals
  eps = reg.yx$residuals*sqrt(n/(n - length(M.hat) - 1))
  
  alpha.hat = reg.yx$coefficients[2]
  alpha.se = sqrt(1/n * 1/mean(nu^2) * mean(eps^2 * nu^2) * 1/mean(nu^2))
  
  result = list(alpha = alpha.hat,
                se = alpha.se, 
                model.size = length(M.hat), 
                M.hat = M.hat,
                M.dx = M.dx, 
                M.yx = M.yx)
  
  return(result)
}

