# Clean workspace
rm(list=ls())

source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/l2n.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/soft.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/myBinarySearch.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/SparseCCA_fast.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/CCAalgorithm.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/CCA_fast.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/fishertran.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/MSqrt.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/fastSVD.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/checkVs.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/get_stability_fast.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/scca_permute_fast.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/scca_permute_fast_one.R')
source('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/optimise_c1_c2_perms_fast.R')

# Create data
u     <- matrix(c(rep(1,10),rep(0,10)),ncol=1)
w     <- matrix(c(rep(0,10),rep(1,10)),ncol=1)
v1    <- matrix(c(rep(1,100),rep(0,900)),ncol=1)
v2    <- matrix(c(rep(0,50),rep(1,50),rep(0,200),rep(1,100),rep(0,600)),ncol=1)
v3    <- matrix(c(rep(0,900),rep(1,100)),ncol=1)
v4    <- matrix(c(rep(0,850),rep(1,100),rep(0,50)),ncol=1)
x     <- u%*%t(v1)+w%*%t(v3) + matrix(rnorm(u%*%t(v1)+w%*%t(v3),20*1000),ncol=1000)
y     <- u%*%t(v2)+w%*%t(v4) + matrix(rnorm(u%*%t(v2)+w%*%t(v4),20*1000),ncol=1000)

# Information
Kv      <- 10
typevar <- "standard"

penalx <- c(seq(0.3,0.9,by=0.1))
penalz <- c(seq(0.3,0.9,by=0.1))

# Find best paramters
perm.out <- optimise_c1_c2_perms_fast(x,y,penaltyxs=penalx,penaltyzs=penalz,nperms=1000,standardize=TRUE,upos=TRUE,vpos=TRUE)
print(sprintf('Best penalty x:>> %g',perm.out$bestpenaltyx))
print(sprintf('Best penalty z:>> %g',perm.out$bestpenaltyz))

# Run SCCA
outperm <- CCA_fast(x,y,K=Kv,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init, standardize=TRUE,upos=TRUE,vpos=TRUE)
print(sprintf('Correlation:>> %g',outperm$cors))

# FPR and FNR
# nzx <- which(v1 != 0)
# zx  <- which(v1 == 0)
# pu  <- which(outperm$u!=0)
# nu  <- which(outperm$u==0)
# tpu <- length(intersect(pu,nzx))
# tnu <- length(intersect(nu,zx))
# fpu <- length(intersect(pu,zx))
# fnu <- length(intersect(nu,nzx))
# fpru <- fpu/(fpu+tnu)
# fnru <- fnu/(fnu+tpu)
# 
# nzy <- which(v2 != 0)
# zy  <- which(v2 == 0)
# pv  <- which(outperm$v!=0)
# nv  <- which(outperm$v==0)
# tpv <- length(intersect(pv,nzy))
# tnv <- length(intersect(nv,zy))
# fpv <- length(intersect(pv,zy))
# fnv <- length(intersect(nv,nzy))
# fprv <- fpv/(fpv+tnv)
# fnrv <- fnv/(fnv+tpv)

# Save output
# u = outperm$u
# v = outperm$v
# t <- u[,1]
# w <- v[,1]
# write.csv(t,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv1_uperm_grey050.csv')
# write.csv(w,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv1_vperm_grey050.csv')

filename <- "/home/k1340324/Data_Analysis/GAPF/SCCA_rev/results_toydata_Kv10_permpos_grey050.RData"
save.image(file=filename)
