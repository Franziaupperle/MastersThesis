
generate.sigma = function(p, type) {
  
  #Independent
  if (type == "ind") {
    Sigma = matrix(data = 0,
                   nrow = p,
                   ncol = p)
    diag(Sigma) = 1
  }
  
  #toeplitz matrix
  else if (type == "toeplitz9") {
    Sigma = matrix(data = 0,
                   nrow = p,
                   ncol = p)
    p.vec = 1:p
    Sigma = 0.9 ^ (toeplitz(p.vec) - 1)
  }
  
  #equal correlation matrix 1
  else if (type == "equalcorr9") {
    Sigma = matrix(data = 0,
                   nrow = p,
                   ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] = 0.9
      }
    }
    diag(Sigma) = 1
  }
  
  #equal correlation matrix 2
  else if (type == "equalcorr3") {
    Sigma = matrix(data = 0,
                   nrow = p,
                   ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] = 0.3
      }
    }
    diag(Sigma) = 1
  }
  
  return(Sigma)
}
