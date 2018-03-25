# This script simulates HVGs involving different abundances.

library(scran)
source("../functions.R")
ngenes <- 1000
ncells <- 500
nspikes <- 500 # For stability, for the time being.

id <- rep(c("high", "low"), each=ngenes/2)
collected.l <- collected.r <- collected.i <- vector("list", 20)
for (it in seq_len(20)) { 
    mus <- rbind(matrix(rep(c(20, 10), c(10, ncells-10)), nrow=ngenes/2, ncol=ncells, byrow=TRUE),
                 matrix(rep(c(2, 0.1), c(10, ncells-10)), nrow=ngenes/2, ncol=ncells, byrow=TRUE))
    disp <- computeDisp(mus)
    
    sim.genes <- matrix(rnbinom(ngenes*ncells, mu=mus, size=1/disp), nrow=ngenes, ncol=ncells)
    sim.spikes <- sampleCounts(ngenes=nspikes, nsamples=ncells, scaling=1)
    counts <- rbind(sim.genes, sim.spikes)
    is.spike <- seq_len(nspikes) + ngenes
    
    my.env <- detectHVGs(counts, is.spike, chosen=seq_len(ngenes))
    collected.l[[it]] <- split(log10(my.env$ref$p.value[-is.spike]), id)
    collected.r[[it]] <- split(pmax(-15, log10(my.env$out.raw$p.value[-is.spike])), id)
    collected.i[[it]] <- split(pmax(-15, log10(my.env$out.imp$p.value[-is.spike])), id)
}

reorder <- c(seq(1, 40, by=2), seq(2, 40, by=2))
color <- rep(c("salmon", "aquamarine"), each=20)

pdf("relative_power.pdf")
boxplot(unlist(collected.l, recursive=FALSE)[reorder], col=color, ylab="-log p-value", main="Variance of logs")
boxplot(unlist(collected.r, recursive=FALSE)[reorder], col=color, ylab="-log p-value", main="Brennecke")
boxplot(unlist(collected.i, recursive=FALSE)[reorder], col=color, ylab="-log p-value", main="Improved CV2")
dev.off()
