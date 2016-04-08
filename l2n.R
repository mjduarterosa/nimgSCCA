l2n <- function (vec) 
{
  a <- sqrt(sum(vec^2))
  if (a == 0) 
    a <- 0.05
  return(a)
}