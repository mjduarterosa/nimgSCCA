scca_permute_fast_one <- function (x, z, penaltyxs, penaltyzs, niter, v, nperms, standardize, upos, uneg, vpos, 
                               vneg)
{
  if (standardize) {
    x <- scale(x, TRUE, TRUE)
    z <- scale(z, TRUE, TRUE)
  }
  v <- checkVs(v, x, z, 1)
  
  # Initialise Arrays
  # # # # # # # # # #
  npenaltyx <- length(penaltyxs)
  npenaltyz <- length(penaltyzs)
  ccperms   = nnonzerous.perms = nnonzerovs.perms = stabperms = array(NA,dim=c(npenaltyx,nperms))
  ccs       = nnonzerous = nnonzerovs = stab = array(0,dim=c(npenaltyx))
  
  for (i in 1:nperms) {
    
    print(sprintf('Permutation:>> %d out of %d',i,nperms))
    
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
    
    # # # # # # #
    #  MY CODE! #
    # # # # # # #
    
    for (j in 1:npenaltyx) {      
      
      if (i == 1) {
        
        out <- CCA_fast(x, z, 
                        penaltyx = penaltyxs[j], penaltyz = penaltyzs[j], 
                        niter = niter, 
                        v = v, upos = upos, uneg = uneg, 
                        vpos = vpos, vneg = vneg, standardize = FALSE)    
        
    
        
        
        nnonzerous[j] <- sum(out$u != 0)
        nnonzerovs[j] <- sum(out$v != 0)
        if (mean(out$u == 0) != 1 && mean(out$v == 0) != 
              1) {
          ccs[j] <- cor(x %*% out$u, z %*% out$v)
        }
        else {
          ccs[j] <- 0
        }
      }
      out <- CCA_fast(x[sampx, ], z[sampz, ], penaltyx = penaltyxs[j], penaltyz = penaltyzs[j], 
                      niter = niter, 
                      v = v, upos = upos, uneg = uneg, 
                      vpos = vpos, vneg = vneg, standardize = FALSE)
      

      
      nnonzerous.perms[j,i] <- sum(out$u != 0)
      nnonzerovs.perms[j,i] <- sum(out$v != 0)
      if (mean(out$u == 0) != 1 && mean(out$v == 0) != 
            1) {
        ccperms[j,i] <- cor(x[sampx, ] %*% out$u, z[sampz, 
                                                      ] %*% out$v)
      }
      else {
        ccperms[j,i] <- 0
      }  
    }    
  }
  
  cc.norm     <- fishertran(ccs)
  ccperm.norm <- fishertran(ccperms)
  
  zstats <- (cc.norm - rowMeans(ccperm.norm))/(apply(ccperm.norm, 1, sd) + 0.05)
  pvals  <- apply(sweep(ccperms, 1, ccs, "-") >= 0, 1, mean)
  

  
  bpenaltyx <- penaltyxs[which.max(zstats)]
  bpenaltyz <- penaltyzs[which.max(zstats)]
  
  
  results <- list(zstats = zstats, penaltyxs = penaltyxs, penaltyzs = penaltyzs, 
                  bestpenaltyx = bpenaltyx, bestpenaltyz = bpenaltyz, 
                  cors = ccs, corperms = ccperms, ft.cors = cc.norm, ft.corperms = rowMeans(ccperm.norm), 
                  nnonzerous = nnonzerous, nnonzerovs = nnonzerovs, nnonzerous.perm = rowMeans(nnonzerous.perms), 
                  nnonzerovs.perm = rowMeans(nnonzerovs.perms), 
                  v.init = v, pvals = pvals, nperms = nperms, 
                  pvalbestz = pvals[which.max(zstats)])
  return(results)
}