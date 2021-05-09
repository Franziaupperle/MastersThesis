
# THIS FILE CONTAINS THE METHODS USED FOR CONDUCTING THE ORACLE ESTIMATOR AND THE REMAINING 5 METHODS


# This method calculates upper and lower bounds of the 95%-CI 
# It further proves, if the true value of alpha lies within the 95%-CI
# This is also known as the coverage probability

# alpha.est: estimated value of alpha
# alpha.se: estimated standard error of alpha
# alpha0: true value of alpha
ci.coverfun = function(alpha.est, alpha.se, alpha0){
  ci.lower = alpha.est - 1.96*alpha.se
  ci.upper = alpha.est + 1.96*alpha.se
  truefalse = ifelse(alpha0 >= ci.lower & alpha0 <= ci.upper, 1, 0)
  # returns, if the true alpha lies within the 95%-CI, upper bound of 95%-CI and lower bound of 95%CI
  result = list(truefalse = truefalse, ci.upper = ci.upper, ci.lower = ci.lower)
  return(result)
}


# This function performs the Oracle estimation

# alpha0: true value of alpha
# y: dependent variable, 
# d: treatment
# x: regressor matrix, 
# m_x0: index vector representing the true model
oracle = function(alpha0, y, d, x, m_x0){
  
  fit = lm(y ~ d + x[,m_x0])
  alpha.est = fit$coefficients[2]
  alpha.se = summary(fit)$coefficients[2,2]
  ci = ci.coverfun(alpha.est, alpha.se, alpha0)$truefalse
  # returns estimated value of alpha, its standard error, information about the 95%-CI and the whole model
  result = list(alpha = alpha.est, se = alpha.se,
                ci = ci, fit = fit)
  return(result)
}


# This function performs Projection-on-Double-Selection

# alpha0: true value of alpha, 
# y: dependent variable, 
# d: treatment, 
# x: regressor matrix
Projection = function(alpha0, y, d, x) {
  
  p = length(x[1,])
  n = length(y)
  
  M = seq(1,p)
  # select a set of variables M_hat_D which are associated with D
  fit.dx = rlasso(d ~ x)$index 
  # index of selected variables (logical vector)
  if(sum(fit.dx) == 0){
    fit.dx = 1
  }
  
  # FRISCH WAUGH THEOREM (line 65 - 66)
  M.dx = M[fit.dx] 
  
  # to remove the components associated with D, (X, Y) is projected onto a space which is orthogonal to the space spanned by D and X_M_hat_D
  x.res = lm(x[,-M.dx]~ d + x[,M.dx])$res 
  y.res = lm(y ~ d + x[,M.dx])$res 
  
  # selecting additional variables that are expected to have low correlation with D and so the
  # over-fitting bias can be controlled. Thus, select a model fit.yx for regressing y.res on x.res
  fit.yx = rlasso(y.res ~ x.res)$index 
  temp0 = M[-M.dx]; 
  if(sum(fit.yx) == 0){ 
    fit.yx = 1
  }
  
  M.yx = temp0[fit.yx] 
  
  # unite relevant covariates from the two selection steps 
  M.hat = c(M.dx, M.yx)
  
  # Auxiliary regression needed for estimating the standard error of the estimated alpha
  reg.dx = lm(d ~ x[,M.hat])
  # Step 4: Regress Y on D and M.hat to get alpha.hat, which is the estimated coefficient of D
  reg.yx = lm(y ~ d + x[,M.hat]) 
  
  # Standard error of alpha.HAT is calculated (line 85 - 86)
  nu = reg.dx$residuals
  eps = reg.yx$residuals/sqrt(n/(n - length(M.hat) - 1))
  
  # get estimator for alpha.hat
  alpha.hat = reg.yx$coefficients[2]
  # get its standard error
  alpha.se = sqrt(1/n * 1/mean(nu^2) * mean(eps^2 * nu^2) * 1/mean(nu^2))
  
  # returns estimated value of alpha, its standard error and the model size of the selected variables
  result = list(alpha = alpha.hat, se = alpha.se,
                model.size = length(M.hat))
  
  return(result)
  
}


# This function calculates the standard errors for the estimated alpha in R-Split and PODS-Split

# Y.count: Matrix with dimension (B x n) with 0-1 entries. 1 indicates iteration j with j = 1,...,B contains observation k with k = 1,..., n
# alpha.est: vector of dimension (B x 1) that contains estimated values of alpha for each iteration
# n: sample size
# split.size: size of the subsample that was used for model selection
# for the calculation, see Wang et al (2019) equation (5)
IF.varestbiascorr = function(Y.count, alpha.est, n, split.size){
  Brep = length(alpha.est)
  n2 = n-split.size
  
  cov.vec = t(sweep(Y.count, 2, apply(Y.count, 2, mean))) %*% (alpha.est - mean(alpha.est))/Brep
  var.est = (sum(cov.vec^2))
  r.corr = (n/split.size)^2 * ((n-1)/n)^2
  var.est0 = var.est * r.corr
  
  var.bias = sum( (alpha.est - mean(alpha.est))^2 )/Brep^2 * n  * n2/(n-n2)
  # returns standard error for estimated alpha
  return( sqrt(abs(var.est0 - var.bias))) 
  
}


# This function selects the model for R-Split and both stages of PODS-Split 
# It contains first variable screening

# y: dependent variable, 
# x: regressor matrix,
# m0: 1 for R-SPlit, 0 for PODS-Split
# size.holp: number of variables kept in the HOLP-screening step
# default: TRUE for R.Split, FALSE for PODS-Split
# dfmax: upper bound of the model size
# dfmin: lower bound of the model size
HolpAlasso.set = function(y, x, m0, size.holp, default, dfmax, dfmin){
  
  p = length(x[1,])
  n = length(y)
  
  index.p = seq(1,p,1)
  
  # prepare data for HOLP estimation
  X = scale(x)
  Y = y - mean(y)
  # HOLP ridge estimator
  OLS = t(X) %*% solve(X %*% t(X) + diag(n) * 1, Y) 
  ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
  # just keep the size.holp highest variables
  index.Holp = ranking$ix[1:size.holp] 
  
  if(default){
    # delete duplicates, m0=1 is used here to ensure that 
    # the treatment effect d (which is the first column of z) is exactly once in M.hat for R-Split
    index.Holp = unique(c(m0,index.Holp)) 
 
  }
  
  # build new matrix out of old matrix x with the as relevant identified indices from screening
  x.screen = x[,index.Holp] 
  # define penalty weights vector for the adaptive LASSO
  penality.weight = 1/abs(OLS[index.Holp,1]) 
  # treatment gets penatly weight 0
  if(default){
    penality.weight[match(m0,index.Holp)] = 0 

  }
  
  # perform LASSO with cross-validation for calculating lambda
  fit.lasso = cv.glmnet(x = x.screen, y = y, penalty.factor = penality.weight, 
                        standardize = TRUE, intercept = TRUE, pmax = dfmax) 
  
  # nzero: number of non-zero coefficients at each LAMBDA 
  # get lambda.dfmin: index of the lambda-vector where model size is at least of size 'dfmin'
  lambda.dfmin = as.numeric(which(fit.lasso$nzero >= dfmin)[1]) 
  # get lambda from just discovered index
  lambda.dfmin = fit.lasso$lambda[lambda.dfmin] 
  # get value of lambda that delivers minimum cvm
  # note: this value could be larger than lambda.dfmin and thus delivering a smaller (maybe too small model)
  lambda.lasso = fit.lasso$lambda.min 
   # get the lambda that delivers the minimum cvm with regard to the minimum for selected model size
  lambda = min(lambda.lasso, lambda.dfmin)
  
  if(is.na(lambda)){
    lambda = min(fit.lasso$lambda) 
  }
  
  # rebuilding the model using glmnet() function
  # glmnet needs to be used again to refit the model with the newly selected lambda
  fit.lasso0 = glmnet(x = x.screen, y = y, standardize = TRUE, intercept = TRUE, 
                       lambda = lambda, penalty.factor = penality.weight)  
  # get nonzero coef of refitted model
  theta.lasso = as.vector(coef(fit.lasso0)) 
  # delete the estimated intercept of the LASSO regression with glmnet
  index.lasso = theta.lasso[-1]!=0 
  
  # get the final selection indices
  M.hat = sort(index.Holp[index.lasso], decreasing = FALSE)
  
  # if no variable was selected, R-Split keeps at least the treatment
  if(default){
  if(sum(index.lasso) == 0){
    M.hat = 1
  }
  }
  
  return(M.hat)
}


# This function selects the variables of Double-Stability and thus is applied twice per iteration

# y: dependent variable, 
# x: regressor matrix
db.stab.set = function(x, y){
  
  p = length(x[1,])
  index.Holp = seq(1,p,1)
  
  # perform LASSO with cross-validation for calculating lambda
  fit.lasso = cv.glmnet(x = x, y = y, 
                        standardize = TRUE, intercept = TRUE) 
  
  # get value of lambda that delivers minimum cvm
  lambda = fit.lasso$lambda.min 
  
  # rebuilding the model using glmnet() function
  fit.lasso0 = glmnet(x = x, y = y, standardize = TRUE, intercept = TRUE, 
                      lambda = lambda) 
  # get nonzero coef of refitted model
  theta.lasso = as.vector(coef(fit.lasso0)) 
  # delete the estimated intercept of the LASSO regression with glmnet
  index.lasso = theta.lasso[-1]!=0 
  
  # get the final selection indices
  M.hat = sort(index.Holp[index.lasso], decreasing = FALSE)
  
  # return NULL, if no variable was selected
  if(sum(index.lasso) == 0){
    M.hat = NULL
  }
  
  return(M.hat) 
}

# This function processes model selection for PODS-Split

# y: dependent variable, 
# d: treatment variable,
# x: regressor matrix,
# size.holp: number of variables kept in the HOLP-screening step,
# dfmax: upper bound of the model size
# dfmin: lower bound of the model size
# HolpAlasso.set.fun: make function HolpAlasso.set available for parallelization
Projection.set = function(y, d, x, size.holp, dfmax, dfmin, HolpAlasso.set.fun) {
  
  p = length(x[1,])
  n = length(y)
  
  # Perform HOLP + adaptive LASSO model selection for variables associated with d
  MD.hat = HolpAlasso.set.fun(y = d, x = x, m0 = 1, size.holp = size.holp, default = FALSE, 
                          dfmax = dfmax, dfmin = dfmin)
  
  # include the projection step as in PODS
  x.res = lm(x[,-MD.hat] ~ d + x[,MD.hat])$res
  y.res = lm(y ~ d + x[,MD.hat])$res
  
  index.p = seq(1,p)
  
  # Perform HOLP + adaptive LASSO model selection for variables associated with y
  MX.hat = HolpAlasso.set.fun(y = y.res, x = x.res, m0 = 1, size.holp = size.holp, default = FALSE,
                          dfmax = dfmax, dfmin = dfmin)
  
  temp0 = index.p[-MD.hat]
  MX.hat = temp0[MX.hat] 
  
  # get indices for final selected variables
  M.hat = unique(sort(c(MD.hat,MX.hat)))
  
  return(M.hat)
  
}












