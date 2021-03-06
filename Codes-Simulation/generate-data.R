
# FUNCTION TO GENERATE ALL RELEVANT DATA 

# pf: selected parameters
# Rd2: coefficient of determination in the treatment equation
# Ry2: coefficient of determination in the main equation
generate.data = function(pf, Rd2, Ry2) {
  
  # get parameter
  Sigma = pf$Sigma
  gamma = pf$gamma
  beta = pf$beta
  m_x0 = pf$m_x0
  split.size = pf$split.size
  size.holp = pf$size.holp
  alpha0 = pf$alpha0
  p = pf$p
  n = pf$n
  
  # generate covariates and error terms
  x = mvrnorm(n = n,
              mu = rep(0, p),
              Sigma = Sigma)
  nu = rnorm(n, 0, 1)
  w = rnorm(n, 0, 1)
  
  # calculate disturbance terms to get desired Rd2 and Ry2
  
  vec.d = seq(0.1, 1.6, by = 0.1)
  Rd.2 = NULL
  vec.y = seq(0.1, 1.6, by = 0.1)
  Ry.2 = NULL
  
  for (j in 1:length(vec.d)) {
    d = 0.5  + x %*% gamma * vec.d[j] + nu
    Rd.2[j] = 1 - var(nu) / var(d)
  }
  
  cd = vec.d[which.min(abs(Rd.2 - Rd2))]
  # calculate the final treatment d
  d = 0.5 + x %*% gamma * cd + nu
  
  for (k in 1:length(vec.y)) {
    y = 1  + alpha0 * d + x %*% beta * vec.y[k] + w
    Ry.2[k] = 1 - (var(w) + alpha0 ^ 2 * var(nu)) / var(y)
  }
  
  cy =  vec.y[which.min(abs(Ry.2 - Ry2))]
  # calculate the final dependend variable
  y = 1 + alpha0 * d + x %*% beta * cy + w
  
  print(cd)
  print(cy)
  
  # return the generated data 
  result = list(x = x,
                d = d,
                y = y)
  
  return(result)
  
}
  
