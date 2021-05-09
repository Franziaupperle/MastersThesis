
# R-SPLIT

# y: dependent variable
# d: treatment variable
# x: covariates
# B: number of subsampling iterations
Split.smooth <- function(y, d, x, B) {
  
  z = cbind(d,x)
  n = length(y)
  p = length(x[1,])
  
  # declare variables for the results
  alpha.est = NULL
  alpha.se = NULL
  M.all = c()
  model.size = c()
  
  # declare matrix used to store the selection indices of the sample
  Y.count = matrix(data = NA, nrow = B, ncol = n)
  # declare matrix used to store the selection indices of the control variables
  C.count = matrix(data = NA, nrow = B, ncol = p)
  
  for (b in 1:B) {
    
    set.seed(b+12345)
    # start repeated subsampling, model selection and treatment estimation
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
    
    # select model with LASSO Cross-Validation or Rigorous LASSO, which one is chosen is selected in function select.model()
    M.hat = select.model(y = y.s, x = x.s)
    # refit the model
    alpha.est[b] = lm(y.t ~ d.t + x.t[,M.hat])$coef[2]
    
    model.size[b] = length(M.hat)
    
    cat("Iteration number: ")
    print(b)
    
    M.all = append(M.all, M.hat)
    C.count[b,] = tabulate(M.hat, p)
    
  }
  
    # get relative inclusion frequencies of each variable over all subsampling steps
    C.rel = colMeans(C.count)

    M.hat09 = which(C.rel > 0.9)
    M.hat08 = which(C.rel > 0.8)
    M.hat07 = which(C.rel > 0.7)
    M.hat06 = which(C.rel > 0.6)
  
    M.Freq = as.data.frame(table(M.all)/B)

  # get the smoothed estimator and its standard error
  alpha.smoothed = mean(alpha.est)
  alpha.se = IF.varestbiascorr(Y.count, alpha.est, n, 0.7*n)
  
  # return results
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

# y: dependent variable
# d: treatment variable
# x: covariates
# B: number of subsampling iterations
PODS.Split <- function(y, d, x, B) {
  
  n = length(y)
  p = length(x[1,])
  
  # declare matrix used to store the selection indices of the sample
  Y.count = matrix(data = NA, nrow = B, ncol = n)
  # declare matrix used to store the selection indices of the control variables
  C.count = matrix(data = NA, nrow = B, ncol = p)
  z = cbind(d,x)
  
  # declare variables for the results
  model.size = c()
  alpha.est = c()
  M.all = c()
  
  for(b in 1:B) {
    
    set.seed(b+12345)
    # start repeated subsampling, model selection and treatment estimation
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
    
    # select the model
    M.hat = Projection.set(y = y.s, d = d.s, x = x.s)
    # refit the model
    alpha.est[b] = lm(y.t ~ d.t + x.t[,M.hat])$coefficients[2]
    model.size[b] = length(M.hat)
    
    cat("Iteration no: ")
    print(b)
    
    M.all = append(M.all, M.hat)
    C.count[b,] = tabulate(M.hat, p)
    
  }

  # get relative inclusion frequencies of each variable over all subsampling steps
  C.rel = colMeans(C.count)
  
  M.hat09 = which(C.rel > 0.9)
  M.hat08 = which(C.rel > 0.8)
  M.hat07 = which(C.rel > 0.7)
  M.hat06 = which(C.rel > 0.6)
  
   # get the smoothed estimator and its standard error
  alpha.smoothed = mean(alpha.est)
  alpha.se = IF.varestbiascorr(Y.count, alpha.est, n, 0.7*n)
  
  # return results
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

# y: dependent variable
# d: treatment variable
# x: covariates
# B: number of subsampling iterations
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
  # declare matrix used to store the selection indices of the control variables associated with the rent price
  Cy.count = matrix(data = NA, nrow = B, ncol = p)
  # declare matrix used to store the selection indices of the control variables associated with the rental brake
  Cd.count = matrix(data = NA, nrow = B, ncol = p)
  
  for(b in 1:B) {
    
    set.seed(b+12345)
    # start repeated subsampling, model selection and treatment estimation
    index.subsam = sample(seq(1,n), 0.7*n)
    y.s = y[index.subsam]
    z.s = z[index.subsam,]
    d.s = z.s[,1] 
    x.s = z.s[,-1]
    
    # DOUBLE SELECTION with LASSO Cross-Validation or Rigorous LASSO, which one is chosen is selected in function select.model()
    M.dx = select.model(y = d.s, x = x.s)
    M.yx = select.model(y = y.s, x = x.s)
    
    M.yx.all = append(M.yx.all, M.yx)
    M.dx.all = append(M.dx.all, M.dx)
    cat("Iteration no: ")
    print(b)
    Cy.count[b,] = tabulate(M.yx, p)
    Cd.count[b,] = tabulate(M.dx, p)
    
  }
  
  # get relative inclusion frequencies of the variable selection from regression on the rent price 
  Cy.rel = colMeans(Cy.count)
  # get relative inclusion frequencies of the variable selection from regression on the rental brake
  Cd.rel = colMeans(Cd.count)

  # get indices of variables with a reltive inclusion frequency higher 0.95
  Cy095 = which(Cy.rel > 0.95)
  Cd095 = which(Cd.rel > 0.95)
  M.hat095 = unique(sort(c(Cy095,Cd095)))
  # get indices of variables with a reltive inclusion frequency higher 0.9
  Cy09 = which(Cy.rel > 0.9)
  Cd09 = which(Cd.rel > 0.9)
  M.hat09 = unique(sort(c(Cy09,Cd09)))
  # get indices of variables with a reltive inclusion frequency higher 0.8
  Cy08 = which(Cy.rel > 0.8)
  Cd08 = which(Cd.rel > 0.8)
  M.hat08 = unique(sort(c(Cy08,Cd08)))
  # get indices of variables with a reltive inclusion frequency higher 0.7
  Cy07 = which(Cy.rel > 0.7)
  Cd07 = which(Cd.rel > 0.7)
  M.hat07 = unique(sort(c(Cy07,Cd07)))
  # get indices of variables with a reltive inclusion frequency higher 0.6
  Cy06 = which(Cy.rel > 0.6)
  Cd06 = which(Cd.rel > 0.6)
  M.hat06 = unique(sort(c(Cy06,Cd06)))
  
  # refit the model with indices of variables with a reltive inclusion frequency higher 0.6
  lm06 = lm(y ~ d + x[,M.hat06]) # Post-lasso OLS step
  alpha.est06 = lm06$coef[2]
  lm_HC106 = coeftest(lm06, vcov = vcovHC(lm06, type="HC1"))
  alpha.se06 = lm_HC106[2,2]

  # refit the model with indices of variables with a reltive inclusion frequency higher 0.7
  lm07 = lm(y ~ d + x[,M.hat07]) # Post-lasso OLS step
  alpha.est07 = lm07$coef[2]
  lm_HC107 = coeftest(lm07, vcov = vcovHC(lm07, type="HC1"))
  alpha.se07 = lm_HC107[2,2]

  # refit the model with indices of variables with a reltive inclusion frequency higher 0.8
  lm08 = lm(y ~ d + x[,M.hat08]) # Post-lasso OLS step
  alpha.est08 = lm08$coef[2]
  lm_HC108 = coeftest(lm08, vcov = vcovHC(lm08, type="HC1"))
  alpha.se08 = lm_HC108[2,2]

  # refit the model with indices of variables with a reltive inclusion frequency higher 0.9
  lm09 = lm(y ~ d + x[,M.hat09]) # Post-lasso OLS step
  alpha.est09 = lm09$coef[2]
  lm_HC109 = coeftest(lm09, vcov = vcovHC(lm09, type="HC1"))
  alpha.se09 = lm_HC109[2,2]
  
  # REFIT THE MODEL WITH indices of variables with a reltive inclusion frequency higher 0.95
  lm095 = lm(y ~ d + x[,M.hat095]) # Post-lasso OLS step
  alpha.est095 = lm095$coef[2]
  lm_HC1095 = coeftest(lm095, vcov = vcovHC(lm095, type="HC1"))
  alpha.se095 = lm_HC1095[2,2]
  
  # return results
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
  
  # Perform model selection with either Rigorous Lasso (line 320 - 326) or Cross-Validation (line 329 - 342)
  # one LASSO method always has to be commented out 

  
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
# used as in the simualtion studies 

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
# used as in the simulation studies

Projection = function(y, d, x) { 
  
  p = length(x[1,])
  n = length(y)
  
  M = seq(1,p)
  fit.dx = rlasso(d ~ x)$index 
  
  if(sum(fit.dx) == 0){ 
    fit.dx = 1
  }
  
  M.dx = M[fit.dx] 
  
  x.res = lm(x[,-M.dx]~ d + x[,M.dx])$res 
  y.res = lm(y ~ d + x[,M.dx])$res 
  
  fit.yx = rlasso(y.res ~ x.res)$index  
  
  temp0 = M[-M.dx]; 
  if(sum(fit.yx) == 0){ 
    fit.yx = 1
  }
  
  M.yx = temp0[fit.yx] 
  
  M.hat = sort(c(M.dx, M.yx))
  
  reg.dx = lm(d ~ x[,M.hat])
  reg.yx = lm(y ~ d + x[,M.hat]) 
  
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

