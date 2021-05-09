

## Function to generate the main parameters as requested in "Results.R"

generate.parameter = function(n, p, type, beta, gamma) {
  
  ## Parameters
  
  ## treatment effect
  alpha0 = 1.5
  ## number of variable kept after a screening step - HOLP
  size.holp = 300
  ## sample size for model selection
  split.size = 0.7 * n
  
  Sigma = generate.sigma(p, type)
  
  if (beta == "sparse") {
    beta = c(1, 1, 1, 1, rep(0, p - 4))
    m_x0 = c(1, 2, 3, 4) # m_x0 is for Oracle estimation, that knows the true model with m_x0 representing the true model
  } else if (beta == "moderately sparse") {
    beta = c(rep(c(5, 1), each = 10), rep(0, p - 20))
    m_x0 = c(seq(1:20))
  } else if (beta == "approximately sparse") {
    beta = c(1 / (seq(p) ^ 2))
    m_x0 = NULL # m_x0 not necessary due to the non existent oracle estimator for this setting
  } else if (beta == "dense") {
    beta = c(1 / (sqrt(seq(1, p))))
    m_x0 = NULL # m_x0 not necessary due to the non existent oracle estimator for this setting
  }
  
  if (gamma == "sparse") {
    gamma = c(0, 0, 0, 0, 1, 1, 1, 1, rep(0, p - 8))
  } else if (gamma == "dense") {
    gamma = c(1 / (sqrt(seq(1, p))))
  }
  
  result = list(
    Sigma = Sigma,
    gamma = gamma,
    beta = beta,
    m_x0 = m_x0,
    split.size = split.size,
    size.holp = size.holp,
    alpha0 = alpha0,
    p = p,
    n = n
  )
  
  return(result)
  
} 


