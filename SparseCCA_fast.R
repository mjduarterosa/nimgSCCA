SparseCCA_fast <- function (x, y, v, penaltyx, penaltyz, niter,
                            upos, uneg, vpos, vneg) 
{
  vold <- rnorm(length(v))
  u <- rnorm(ncol(x))
  for (i in 1:niter) {
    if (sum(is.na(u)) > 0 || sum(is.na(v)) > 0) {
      v <- rep(0, length(v))
      vold <- v
    }
    if (sum(abs(vold - v)) > 1e-06) {
      unew <- rep(NA, ncol(x))
      
      argu <- matrix(y %*% v, nrow = 1) %*% x
      if (upos) 
        argu <- pmax(argu, 0)
      if (uneg) 
        argu <- pmin(argu, 0)
      lamu <- myBinarySearch(argu, penaltyx * sqrt(ncol(x)))
      su <- soft(argu, lamu)
      u <- matrix(su/l2n(su), ncol = 1)
      
      vnew <- rep(NA, ncol(y))
      
      vold <- v
      argv <- matrix(x %*% u, nrow = 1) %*% y
      if (vpos) 
        argv <- pmax(argv, 0)
      if (vneg) 
        argv <- pmin(argv, 0)
      lamv <- myBinarySearch(argv, penaltyz * sqrt(ncol(y)))
      sv <- soft(argv, lamv)
      v <- matrix(sv/l2n(sv), ncol = 1)
      
    }
  }
  
  d <- sum((x %*% u) * (y %*% v))
  if (sum(is.na(u)) > 0 || sum(is.na(v)) > 0) {
    u <- matrix(rep(0, ncol(x)), ncol = 1)
    v <- matrix(rep(0, ncol(y)), ncol = 1)
    d <- 0
  }
  return(list(u = u, v = v, d = d))
}
