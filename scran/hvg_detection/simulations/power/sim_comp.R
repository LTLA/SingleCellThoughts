# This script simulates HVGs involving different abundances.

library(scran)
source("../functions.R")
ngenes <- 1000
ncells <- 500
nspikes <- 500 # For stability, for the time being.

id <- rep(c("high", "low"), each=ngenes/2)
mus <- rbind(matrix(rep(c(20, 10), c(10, ncells-10)), nrow=ngenes/2, ncol=ncells, byrow=TRUE),
             matrix(rep(c(2, 0.1), c(10, ncells-10)), nrow=ngenes/2, ncol=ncells, byrow=TRUE))
disp <- computeDisp(mus)

sim.genes <- matrix(rnbinom(ngenes*ncells, mu=mus, size=1/disp), nrow=ngenes, ncol=ncells)
sim.spikes <- sampleCounts(ngenes=nspikes, nsamples=ncells, scaling=1)
counts <- rbind(sim.genes, sim.spikes)
is.spike <- seq_len(nspikes) + ngenes
                        
my.env <- detectHVGs(counts, is.spike, chosen=seq_len(ngenes))
pdf("relative_power.pdf")
boxplot(split(log10(my.env$ref$p.value[-is.spike]), id), ylab="-log p-value", Main="Variance of logs")
boxplot(split(pmax(-15, log10(my.env$out.raw$p.value[-is.spike])), id), ylab="-log p-value", main='Brennecke')
boxplot(split(pmax(-15, log10(my.env$out.imp$p.value[-is.spike])), id), ylab="-log p-value", main='Improved CV2')
dev.off()
