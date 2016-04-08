rm(list=ls())
library(PMA)
library(psych)

#load("/Users/mariarosa/Documents/Work/scca/rev_results/results_toydata_Kv10_permpos_grey050.RData")

load("/Users/mariarosa/Documents/Work/scca/rev_results/results_GAPF_Kv10_permpos_grey050.RData")

nperm <- 1000
sizex <- dim(x)[1]
pval  <- 0

cor_cor <- matrix(0,nperm,Kv)
cor_smp <- matrix(0,nperm,Kv)

out <- CCA(x,y,typex="standard",typez="standard",K=Kv,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init,upos=TRUE,vpos=TRUE,standardize=TRUE)

for (p in 1:nperm)

{

print(sprintf("Permutation..........................>>>>: %d out of %d", p, nperm))

smpx     <- sample(sizex)
smpy     <- sample(sizex)

outperm <- CCA(x[smpx,],y,typex="standard",typez="standard",K=Kv,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init,upos=TRUE,vpos=TRUE,standardize=TRUE, trace = FALSE)

for (i in 1:Kv)
{

cor_smp[p,i] <- cor(scale(x[smpx,]) %*% outperm$u[,i], scale(y) %*% outperm$v[,i])
cor_cor[p,i] <- cor(scale(x) %*% out$u[,i], scale(y) %*% out$v[,i])

}

}

p_vals <- colSums(cor_smp>=cor_cor)/nperm
f_cor  <- fisherz(cor_cor[1,])
f_smp  <- fisherz(cor_smp)
c_mean <- colMeans(f_smp)
c_sd   <- apply(f_smp, 2, sd)

z_stat <- (f_cor - c_mean)/c_sd






