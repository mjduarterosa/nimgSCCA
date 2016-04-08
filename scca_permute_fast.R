scca_permute_fast <- function (x, z, penaltyxs, penaltyzs, niter, v, nperms, standardize, upos, uneg, vpos, 
                          vneg)
{
  if (standardize) {
    x <- scale(x, TRUE, TRUE)
    z <- scale(z, TRUE, TRUE)
  }
  v <- checkVs(v, x, z, 1)
  
  #ccperms = nnonzerous.perms = nnonzerovs.perms = matrix(NA, npenaltyx, nperms)
  #ccs = nnonzerous = nnonzerovs = numeric(npenaltyx)
  
  # Initialise Arrays
  # # # # # # # # # #
  npenaltyx <- length(penaltyxs)
  npenaltyz <- length(penaltyzs)
  ccperms = nnonzerous.perms = nnonzerovs.perms = stabperms = array(NA,dim=c(npenaltyx,npenaltyz,nperms))
  ccs = nnonzerous = nnonzerovs = stab = array(0,dim=c(npenaltyx,npenaltyz))
  
  penaltyx.matrix <- matrix(rep(penaltyxs, npenaltyx), nrow=npenaltyx, byrow=F)
  penaltyz.matrix <- matrix(rep(penaltyzs, npenaltyz), nrow=npenaltyz, byrow=T)
  
  for (i in 1:nperms) {
    
    print(sprintf('Permutation:>> %d out of %d',i,nperms))
    
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
    
    # # # # # # #
    #  MY CODE! #
    # # # # # # #
    
    for (j in 1:npenaltyx) {

      for (k in 1:npenaltyz) {
      
        if (i == 1) {
          
          out <- CCA_fast(x, z, 
                     penaltyx = penaltyxs[j], penaltyz = penaltyzs[k], 
                     niter = niter, 
                     v = v, upos = upos, uneg = uneg, 
                     vpos = vpos, vneg = vneg, standardize = FALSE)
          

          
          
          nnonzerous[j,k] <- sum(out$u != 0)
          nnonzerovs[j,k] <- sum(out$v != 0)
          if (mean(out$u == 0) != 1 && mean(out$v == 0) != 
                1) {
            ccs[j,k] <- cor(x %*% out$u, z %*% out$v)
          }
          else {
            ccs[j,k] <- 0
          }
        }
        out <- CCA_fast(x[sampx, ], z[sampz, ], penaltyx = penaltyxs[j], penaltyz = penaltyzs[k], 
                   niter = niter, 
                   v = v, upos = upos, uneg = uneg, 
                   vpos = vpos, vneg = vneg, standardize = FALSE)
        

        
        nnonzerous.perms[j,k,i] <- sum(out$u != 0)
        nnonzerovs.perms[j,k,i] <- sum(out$v != 0)
        if (mean(out$u == 0) != 1 && mean(out$v == 0) != 
              1) {
          ccperms[j,k,i] <- cor(x[sampx, ] %*% out$u, z[sampz, 
                                                        ] %*% out$v)
        }
        else {
          ccperms[j,k,i] <- 0
        }
      }  
    }    
  }
  
  cc.norm     <- fishertran(ccs)
  ccperm.norm <- fishertran(ccperms)
  

  
  zstats <- (cc.norm - apply(ccperm.norm, c(1,2), mean))/(apply(ccperm.norm, c(1,2), sd) + 0.05)
  pvals  <- apply(sweep(ccperms, c(1,2), ccs, "-") >= 0, c(1,2), mean)
  

  
  bpenaltyx <- penaltyx.matrix[zstats == max(zstats)]
  bpenaltyz <- penaltyz.matrix[zstats == max(zstats)]
  
  bpenaltyx <- bpenaltyx[1]
  bpenaltyz <- bpenaltyz[1]
  
  results <- list(zstats = zstats, penaltyxs = penaltyxs, penaltyzs = penaltyzs, 
                  bestpenaltyx = bpenaltyx, bestpenaltyz = bpenaltyz, 
                  cors = ccs, corperms = ccperms, ft.cors = cc.norm, ft.corperms = rowMeans(ccperm.norm), 
                  nnonzerous = nnonzerous, nnonzerovs = nnonzerovs, nnonzerous.perm = rowMeans(nnonzerous.perms), 
                  nnonzerovs.perm = rowMeans(nnonzerovs.perms), 
                  v.init = v, pvals = pvals, nperms = nperms, 
                  pvalbestz = pvals[which.max(zstats)])
  return(results)
}