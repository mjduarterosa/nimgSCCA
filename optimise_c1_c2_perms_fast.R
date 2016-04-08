optimise_c1_c2_perms_fast <- function (x, z, penaltyxs = NULL, penaltyzs = NULL, niter = 3, 
                                  v = NULL, nperms = 25, standardize = TRUE, upos = FALSE, uneg = FALSE, 
                                  vpos = FALSE, vneg = FALSE) 
{
  u <- NULL
  out <- scca_permute_fast(x = x, z = z, penaltyxs = penaltyxs, penaltyzs = penaltyzs, 
                      niter = niter, v = v, nperms = nperms, 
                      standardize = standardize, 
                      upos = upos, uneg = uneg, vpos = vpos, vneg = vneg)
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  class(out) <- "scca_permute_fast"
  return(out)
}