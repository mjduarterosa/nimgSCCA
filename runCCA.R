# Load PMA package
library(PMA)

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

# Import data
x.data    <- as.matrix(read.table('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_HAL_nomean_grey050.txt'))
y.data    <- as.matrix(read.table('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_ARI_nomean_grey050.txt'))
z.data    <- as.matrix(read.table('/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_PLA_nomean_grey050.txt'))

x.data    <- scale(x.data,center=TRUE,scale=FALSE)
y.data    <- scale(y.data,center=TRUE,scale=FALSE)
z.data    <- scale(z.data,center=TRUE,scale=FALSE)

id        <- c(1:8,11:20)
x         <- t(x.data);
x         <- x[id,]
y         <- t(y.data);
z         <- t(z.data);
z         <- z[id,]

x   <- x-z
y   <- y-z

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

# Save output
u = outperm$u
v = outperm$v
t <- u[,1]
w <- v[,1]
write.csv(t,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv10c6_uweights.csv')
write.csv(w,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv10c6_vweights.csv')
t <- u[,2]
w <- v[,2]
write.csv(t,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv10c2_uweights.csv')
write.csv(w,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv10c2_vweights.csv')
t <- u[,3]
w <- v[,3]
write.csv(t,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv10c3_uweights.csv')
write.csv(w,file='/home/k1340324/Data_Analysis/GAPF/SCCA_rev/GAPF_Kv10c3_vweights.csv')

filename <- "/home/k1340324/Data_Analysis/GAPF/SCCA_rev/results_GAPF_Kv10_permpos_grey050.RData"
save.image(file=filename)
