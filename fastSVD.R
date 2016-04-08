fastSVD <- function (x, z) 
{
  xx = x %*% t(x)
  xx2 = MSqrt(xx)
  y = t(z) %*% xx2
  a = try(svd(y), silent = TRUE)
  iter <- 1
  if (class(a) == "try-error" && iter < 10) {
    a = try(svd(y), silent = TRUE)
    iter <- iter + 1
  }
  if (iter == 10) 
    stop("too many tries.")
  v = a$u
  d = a$d
  zz = z %*% t(z)
  zz2 = MSqrt(zz)
  y = t(x) %*% zz2
  a = try(svd(y), silent = TRUE)
  iter <- 1
  if (class(a) == "try-error" && iter < 10) {
    a = try(svd(y), silent = TRUE)
    iter <- iter + 1
  }
  if (iter == 10) 
    stop("too many tries.")
  u = a$u
  return(list(u = u, v = v, d = d))
}