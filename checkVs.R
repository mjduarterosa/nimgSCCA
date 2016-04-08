checkVs <- function (v, x, z, K) 
{
  if (!is.null(v) && !is.matrix(v)) 
    v <- matrix(v, nrow = ncol(z))
  if (!is.null(v) && ncol(v) < K) 
    v <- NULL
  if (!is.null(v) && ncol(v) > K) 
    v <- matrix(v[, 1:K], ncol = K)
  if (is.null(v) && ncol(z) > nrow(z) && ncol(x) > nrow(x)) {
    v <- try(matrix(fastSVD(x, z)$v[, 1:K], ncol = K), silent = TRUE)
    attempt <- 1
    while (class(v) == "try-error" && attempt < 10) {
      v <- try(matrix(fastSVD(x, z)$v[, 1:K], ncol = K), 
               silent = TRUE)
      attempt <- attempt + 1
    }
    if (attempt == 10) 
      stop("Problem computing SVD.")
  }
  else if (is.null(v) && (ncol(z) <= nrow(z) || ncol(x) <= 
                            nrow(x))) {
    attempt <- 1
    v <- try(matrix(svd(t(x) %*% z)$v[, 1:K], ncol = K), 
             silent = TRUE)
    while (class(v) == "try-error" && attempt < 10) {
      v <- try(matrix(svd(t(x) %*% z)$v[, 1:K], ncol = K), 
               silent = TRUE)
      attempt <- attempt + 1
    }
    if (attempt == 10) 
      stop("Problem computing SVD.")
  }
  return(v)
}