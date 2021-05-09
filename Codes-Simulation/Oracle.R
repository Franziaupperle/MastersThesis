
library(caret)
library(lattice)
library(ggplot2)

### Estimation with oracle model selection
# Oracle estimator only used when beta is (moderately) sparse

oracle.mcfun = function(no_runs_mc) {
  
  pf = generate.parameter(
    n = n,
    p = p,
    type = type,
    beta = beta,
    gamma = gamma
  )
  
  bias = c()
  mse = c()
  ci = c()
  ci.laenge = c()
  oracle0 = list(NULL, NULL, NULL, NULL)
  
  for (i in 1:no_runs_mc) {
    set.seed(i + 12345)
    
     # get data
    data = generate.data(pf = pf, Rd2 = Rd2, Ry2 = Ry2)
    x = data$x
    y = data$y
    d = data$d
    m_x0 = pf$m_x0
    alpha0 = pf$alpha0
    
    oracle = oracle(
      alpha0 = alpha0,
      y = y,
      d = d,
      x = x,
      m_x0 = m_x0
    )
    
    oracle0 = Map(c, oracle0, oracle)
    
    bias[i] = oracle$alpha - alpha0
    mse[i] = (oracle$alpha - alpha0) ^ 2 + oracle$se ^ 2
    ci[i] = oracle$ci
    confi.int = ci.coverfun(
      alpha.est = oracle$alpha,
      alpha.se = oracle$se,
      alpha0 = alpha0
    )
    ci.laenge[i] = confi.int$ci.upper - confi.int$ci.lower
    
    cat("Iteration number: ", i, "\n")
    
  }
  
  result = list(
    oracle = oracle0,
    Bias = round(sqrt(n) * mean(bias), 2),
    SE.Bias = round(sqrt(n) * std.error(bias), 2),
    MSE = round(n * mean(mse), 2),
    SE.MSE = round(n * std.error(mse), 2),
    Cover = round(mean(ci), 2),
    SE.Cover = round(std.error(ci), 2),
    Length = round(mean(ci.laenge), 2),
    SE.Length = round(std.error(ci.laenge), 2),
    beta = beta,
    gamma = gamma,
    type = type,
    Ry2 = Ry2,
    Rd2 = Rd2
  )
  
  return(result)
  
}

