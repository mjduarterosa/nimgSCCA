MSqrt <- function (x) 
{
  eigenx <- eigen(x)
  return(eigenx$vectors %*% diag(sqrt(pmax(0, eigenx$values))) %*% 
           t(eigenx$vectors))
}