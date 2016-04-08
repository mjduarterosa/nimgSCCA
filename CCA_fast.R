CCA_fast <- function (x, z, penaltyx = NULL, penaltyz = NULL, K = 1, niter = 15, 
          v = NULL, standardize = TRUE, upos = FALSE, 
          uneg = FALSE, vpos = FALSE, vneg = FALSE) 
{
  
  if (standardize) {
    sdx <- apply(x, 2, sd)
    sdz <- apply(z, 2, sd)
    if (min(sdx) == 0) 
      stop("Cannot standardize because some of the columns of x have std. dev. 0")
    if (min(sdz) == 0) 
      stop("Cannot standardize because some of the columns of z have std. dev. 0")
    x <- scale(x, TRUE, sdx)
    z <- scale(z, TRUE, sdz)
  }
  
  v <- checkVs(v, x, z, K)
  
  out <- CCAalgorithm(x = x, z = z, v = v, 
                      penaltyx = penaltyx, penaltyz = penaltyz, K = K, niter = niter, 
                      upos = upos, 
                      uneg = uneg, vpos = vpos, vneg = vneg)
  out$penaltyx <- penaltyx
  out$penaltyz <- penaltyz
  out$K <- K
  out$niter <- niter
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  out$v.init <- v
  out$cors <- numeric(K)
  for (k in 1:K) {
    if (sum(out$u[, k] != 0) > 0 && sum(out$v[, k] != 0) > 
          0) 
      out$cors[k] <- cor(x %*% out$u[, k], z %*% out$v[, 
                                                       k])
  }
  class(out) <- "CCA_fast"
  return(out)
}