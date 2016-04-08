optimise_c1_c2_perms <- function (x, z, typex = c("standard", "ordered"), typez = c("standard", 
                                                                                    "ordered"), penaltyxs = NULL, penaltyzs = NULL, niter = 3, 
                                  v = NULL, trace = TRUE, nperms = 25, standardize = TRUE, 
                                  chromx = NULL, chromz = NULL, upos = FALSE, uneg = FALSE, 
                                  vpos = FALSE, vneg = FALSE, outcome = NULL, y = NULL, cens = NULL) 
{
  u <- NULL
  typex <- match.arg(typex)
  typez <- match.arg(typez)
  call <- match.call()
  if (is.null(penaltyxs)) 
    penaltyxs <- seq(0.1, 0.7, len = 10)
  if (is.null(penaltyzs)) 
    penaltyzs <- seq(0.1, 0.7, len = 10)
  out <- scca_permute(x = x, z = z, typex = typex, 
                            typez = typez, penaltyxs = penaltyxs, penaltyzs = penaltyzs, 
                            niter = niter, v = v, trace = trace, nperms = nperms, 
                            standardize = standardize, chromx = chromx, chromz = chromz, 
                            upos = upos, uneg = uneg, vpos = vpos, vneg = vneg, 
                            outcome = outcome, y = y, cens = cens)
  out$call <- call
  out$upos <- upos
  out$uneg <- uneg
  out$vpos <- vpos
  out$vneg <- vneg
  class(out) <- "CCA.permute"
  return(out)
}