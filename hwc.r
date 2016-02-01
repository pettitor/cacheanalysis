hwca <- function(l,b,d,rho,C) {
  
  a <- matrix(data = NA, nrow = length(l), ncol = length(rho))
  
  for (m in 1:dim(a)[1]) {
    for (k in 1:dim(a)[2]) {
      a[m,k] <- l[m]/rho[k]*b[m]*d[m]/length(rho)
    }
  }
  return(a)
}

hwc <- function(l,b,d,rho,C) {
  
  a <- matrix(data = NA, nrow = length(l), ncol = length(rho))
  
  for (m in 1:dim(a)[1]) {
    for (k in 1:dim(a)[2]) {
      a[m,k] <- l[m]*b[m]*d[m]/rho[k]/length(rho)
    }
  }
  
  X <- matrix(data = rep(0,length(l)*length(rho)), nrow = length(l), ncol = length(rho))
  
  for (m in 1:dim(X)[1]) {
    for (k in 1:dim(X)[2]) {
      if (sum(X[,k]) >=  C[k] || sum(X[m,]) >= sum(a[m,])) {
        X[m,k] <- 0
      } else {
        X[m,k] <- 1
      }
    }
  }
  return(X)
}